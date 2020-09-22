# Return initialization parameters for the fit
.initializer = function(X, K, tail, init, pareto.shape = list(min.val = 0.01, max.val = 5))
{
  if(is.list(init)) {

    stopifnot(c('beta', 'shape', 'scale') %in% names(init))
    stopifnot(c('a', 'b') %in% rownames(init$beta))

    return(init)
  }

  # Sampler for Beta parameters with random mean.
  # Reject samples if a/ b are negative
  bsamp = function(m = runif(1)) {
    repeat{
      v = runif(n = 1, min = 0.0001, max = 0.01)
      p = c(mobster:::.estBetaParams(m, v), mean=m, var=v)
      if(all(p > 0)) return(p)
    }
  }

  if(is.character(init) && init == 'random')
    beta = sapply(runif(K), bsamp)

  if(is.character(init) && init == 'even')
  {
    stop('Not yet implemented')
    #
    # units = 1/Kbeta
    # init$mean = units * 0:(Kbeta-1) + units/2
    #
    # for(i in 1:Kbeta) {
    #   lower = (i-1) * units
    #   upper = i * units
    #   sample.var = var(X[X > lower & X < upper])
    #   init$var[i] = ifelse(is.na(sample.var), 0.01, sample.var)
    #
    #   params = estBetaParams(init$mean[i], init$var[i])
    #   a[i] = unlist(params$alpha)
    #   b[i] = unlist(params$alpha)
    #
    # }
  }

  if(is.character(init) & init == 'peaks')
  {
    # Compute KDE
    h = hist(X, breaks = seq(0, 1, 0.01), plot = FALSE)

    # Detect peaks
    peaks = mobster:::.find_peaks(h$density, 1)
    x.peaks = (peaks * 0.01)

    # store only peaks above 0.1
    peaks = peaks[x.peaks > 0.1]
    x.peaks = x.peaks[x.peaks > 0.1]
    peakValues = h$density[peaks]

    # Cluster their x-coordinates
    clus = cnt = NULL
    if(length(x.peaks) > K)
    {
      clus = kmeans(x.peaks, K, nstart = 100)
      cnt = as.vector(clus$centers)
    }
    else{
      clus = kmeans(x.peaks, length(x.peaks) - 1, nstart = 100)
      clus = as.vector(clus$centers)
      cnt = c(
        clus,
        runif(K - length(clus))
      )
    }

    # Get Beta peaks from means
    beta = sapply(cnt, bsamp)
  }

  # Prepare tibbles for Betas
  Betas = NULL

  for(i in 1:K)
    Betas =  dplyr::bind_rows(
      Betas,
      tibble::tribble(
        ~cluster, ~type, ~fit.value, ~init.value,
        paste0('C',i), "a", NA, unlist(beta['a', i]) %>% as.vector,
        paste0('C',i), "b", NA, unlist(beta['a', i]) %>% as.vector,
        paste0('C',i), "Mean", NA, unlist(beta['mean', i]) %>% as.vector,
        paste0('C',i), "Variance", NA, unlist(beta['var', i]) %>% as.vector
      )
    )

  # Tail, sample and make a tibble
  Pareto = NULL

  if(tail)
  {
    shape = runif(1, min = pareto.shape$min.val, max = pareto.shape$max.val)
    scale = min(X) - 1e-9
    mv = mobster:::.MeanVarPareto(shape, scale)

    Pareto = tibble::tribble(
      ~cluster, ~type, ~fit.value, ~init.value,
      'Tail', "Shape", NA, shape,
      'Tail', "Scale", NA, scale,
      'Tail', "Mean", NA, mv$mean %>% as.vector,
      'Tail', "Variance", NA, mv$var %>% as.vector
    )
  }


  # Mixing proportions -- uniform distribution
  K = K + as.numeric(tail)

  Clusters = dplyr::bind_rows(Pareto, Betas)

  pi.Pareto = tibble::tribble(
    ~cluster, ~type, ~fit.value, ~init.value,
    'Tail', "Mixing proportion", NA, ifelse(tail, 1/K, 0)
  )

  pi.Betas = NULL
  for(x in unique(Betas$cluster))
    pi.Betas = dplyr::bind_rows(pi.Betas,
                                tibble::tribble(
                                  ~cluster, ~type, ~fit.value, ~init.value,
                                  x, "Mixing proportion", NA, 1/K
                                ))

  dplyr::bind_rows(Pareto, Betas, pi.Pareto, pi.Betas)
}

# Decide if the fit stops or not
.stoppingCriterion = function(i, prevNLL, NLL, prevpi, pi, fit.type, epsilon, isDebug, K)
{
  # Compute NLL difference after ith iteration
  NLL.Diff  <- prevNLL - NLL

  # Fit by MLE
  if(fit.type == 'MLE')
  {
    # Correctness: we know by EM properties that the LL should always increase
    if (NLL.Diff < 0)
    {
      stop("Negative log likelihood increases: ", NLL.Diff, ".\nSomething is wrong, aborting!")
    }

    # Check for convergence
    stopping = NLL.Diff < epsilon
  }

  # Fit by MM
  if(fit.type == 'MM')
  {
    # With MM in general we are not interested in the logLik
    # and we only look for variations in the actual mixing
    # proportions. When the fit is however single-cluster, 
    # we need to check also the logLik because the proportions
    # never change by definition (do the abs beacuse Jensen's ineq. is invalid)
    
    if(K > 1)
    {
      # Compute pi's difference after ith iteration
      pi.Diff  <- abs(prevpi - pi)
  
      # Check for convergence -- all pi's have changed less than epsilon
      stopping = all(pi.Diff < epsilon)
    }
    else
      stopping = abs(NLL.Diff) < epsilon
  }

  # Some printing
  if (isDebug)
  {
    pi = round(pi, 3)
    NLL = round(NLL, 2)
    NLL.Diff = format(NLL.Diff, scientific = T)
    pi.Diff = max(prevpi - pi)

    cat("\t[ step ]", i, "\t[ NLL ]", sprintf("%.13s", NLL), 'delta', NLL.Diff,
        "[ pi ]", paste(pi, collapse = ', '), 'delta', pi.Diff, '\r')
    # " \r")
  }

  return(stopping)
}

# Peaks detection function
.find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

# Use the log sum exp trick for having numeric stability
.log_sum_exp <- function(x) {
  # Computes log(sum(exp(x))
  offset <- max(x)
  s <- log(sum(exp(x - offset))) + offset
  i <- which(!is.finite(s))
  if (length(i) > 0) { s[i] <- offset }
  return(s)
}

# Switch from mean/ var to a and b for Beta
.estBetaParams <- function(mu, var)
{
  if(var == 0) var = 1e-9
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(a = alpha, b = beta))
}

# Switch from a and b to mean/ var for Beta
.MeanVarBeta = function(a, b) {
  m = a/(a+b)
  return(list(mean = m, var = (m*(1-m))/(a+b+1)))
}

# Switch from shape and scale to mean/ var for Pareto
.MeanVarPareto = function(shape, scale) {
  meanPareto = ifelse(shape > 1, (shape * scale) / (shape - 1), Inf)
  varPareto = ifelse(shape > 2, (shape * scale**2) / ((shape - 1)**2 * (shape - 2)), Inf)
  
  return(list(mean = meanPareto, var = varPareto))
}

# Extract Beta parameters
.params_Beta = function(x, init = FALSE)
{
  if(init) x$Clusters$fit.value = x$Clusters$init.value
  
  x$Clusters %>%
    dplyr::filter(cluster != 'Tail', type == 'a' | type == 'b') %>%
    dplyr::select(-init.value) %>%
    tidyr::spread(key = type, value = fit.value)
}

# Extract Pareto parameters
.params_Pareto = function(x, init = FALSE)
{
  if(init) x$Clusters$fit.value = x$Clusters$init.value
  
  suppressWarnings(
    x$Clusters %>%
      dplyr::filter(cluster == 'Tail', type == 'Shape' | type == 'Scale') %>%
      dplyr::select(-init.value) %>%
      tidyr::spread(key = type, value = fit.value)
  )
}

# Extract mixing proportions parameters
.params_Pi = function(x, init = FALSE)
{
  if(init) x$Clusters$fit.value = x$Clusters$init.value
  
  v = x$Clusters %>%
    dplyr::filter(type == 'Mixing proportion') %>%
    dplyr::select(-init.value) %>%
    tidyr::spread(key = type, value = fit.value)
  
  # pi = v$`Mixing proportion`
  # names(pi) = v$cluster
  
  pi = pio:::nmfy(v$cluster, v$`Mixing proportion`)
  
  # ord.pi = c(pi['Tail'], pi[names(pi) != 'Tail'])
  ord.pi = c(
    pi['Tail'],
    pi[sort(names(pi)[names(pi) != 'Tail'], decreasing = FALSE)]
    )
  
  ord.pi
}

# Set Beta parameters
.set_params_Beta = function(fit, a, b)
{
  names.BetaC = names(a)
  if(any(is.null(names.BetaC))) stop("params -- named vector required?")
  
  # Dangerous for ordering of Beta clusters..
  # fit$Clusters[
  #   fit$Clusters$cluster %in% names.BetaC & fit$Clusters$type == 'a',
  #   'fit.value'
  #   ] =  a
  
  # fit$Clusters[
  #   fit$Clusters$cluster %in% names.BetaC & fit$Clusters$type == 'b',
  #   'fit.value'
  #   ] =  b
  
  
  fit$Clusters = fit$Clusters %>%
    dplyr::mutate(
      fit.value = 
        ifelse(
          (cluster %in% names.BetaC) & (type == 'a'), 
          a[cluster],
          fit.value
        )
    )
  
  fit$Clusters = fit$Clusters %>%
    dplyr::mutate(
      fit.value = 
        ifelse(
          (cluster %in% names.BetaC) & (type == 'b'), 
          b[cluster],
          fit.value
        )
    )
  
  
  
  for (s in names.BetaC)
  {
    mv = mobster:::.MeanVarBeta(a[s], b[s])
    
    # fit$Clusters[
    #   fit$Clusters$cluster == s & fit$Clusters$type == 'Mean',
    #   'fit.value'
    #   ] =  mv$mean
    
    # fit$Clusters[
    #   fit$Clusters$cluster == s & fit$Clusters$type == 'Variance',
    #   'fit.value'
    #   ] =  mv$var
    
    fit$Clusters = fit$Clusters %>%
      dplyr::mutate(
        fit.value = 
          ifelse(
            (cluster == s) & (type == 'Mean'), 
            mv$mean,
            fit.value
          )
      )
    
    fit$Clusters = fit$Clusters %>%
      dplyr::mutate(
        fit.value = 
          ifelse(
            (cluster %in% s) & (type == 'Variance'), 
            mv$var,
            fit.value
          )
      )
  }
  
  fit
}

# Set Pareto parameters
.set_params_Pareto = function(fit, shape, scale)
{
  if(!fit$fit.tail) return(fit)
  
  mv = mobster:::.MeanVarPareto(shape, scale)
  
  fit$Clusters[
    fit$Clusters$cluster == 'Tail' & fit$Clusters$type == 'Shape',
    'fit.value'
    ] = shape
  
  fit$Clusters[
    fit$Clusters$cluster == 'Tail' & fit$Clusters$type == 'Scale',
    'fit.value'
    ] =  scale
  
  fit$Clusters[
    fit$Clusters$cluster == 'Tail' & fit$Clusters$type == 'Mean',
    'fit.value'
    ] =  mv$mean
  
  fit$Clusters[
    fit$Clusters$cluster == 'Tail' & fit$Clusters$type == 'Variance',
    'fit.value'
    ] =  mv$var
  
  fit
}

# Set mixing proportions parameters
.set_params_Pi = function(fit, pi)
{
  # Dangerous, not robust to ordering permutation..
  # fit$Clusters[
  #   fit$Clusters$type == 'Mixing proportion',
  #   'fit.value'
  #   ] =  pi
  
  fit$Clusters = fit$Clusters %>%
    dplyr::mutate(
      fit.value = ifelse(
        type == 'Mixing proportion', 
        pi[cluster],
        fit.value
      )
    )
  
  fit
}

# Compute SSE of a fit
.compute_fit_sqerr = function(x, binning = 1e-2)
{
  densities = mobster:::template_density(
    x,
    x.axis = seq(binning, 1-binning, by = binning), # Restricted for numerical errors
    binwidth = binning,
    reduce = TRUE)
  
  densities = tibble::as_tibble(densities)
  densities = densities %>% group_by(x) %>% summarise(y = sum(y), cluster = 'f(x)')
  
  # Empirical density
  empirical = hist(x$data$VAF, breaks = seq(0, 1, binning), plot = FALSE)$density
  empirical = empirical[-length(empirical)]
  empirical = empirical * binning # adjust for binwidth
  
  # Error
  error = densities
  error$cluster = 'e(x)'
  error$y = (densities$y - empirical)^2
  error$cum.y = cumsum(error$y)
  
  error
}

.get_clusters_labels = function(x)
{
  unique(x$data$cluster)
}


# Very fast analysis, raw results but good to get a grasp of the data
# - max 2 clones (1 subclone);
# - random start
# - 2 samples
# - mild epsilon/ maxIter
# - not parallel
template_parameters_fast_setup = function()
{
  return(
    list(
      K = 1:2,
      samples = 2,
      init = 'random',
      tail = c(TRUE, FALSE),
      epsilon = 1e-6,
      maxIter = 100,
      fit.type = 'MM',
      seed = 12345,
      model.selection = 'reICL',
      trace = FALSE,
      parallel = FALSE,
      pi_cutoff = 0.02,
      N_cutoff = 10
    )
  )
}

auto_setup = function(x)
{
  if(x == "FAST")  return(template_parameters_fast_setup())
  
  stop("Auto setup unknown: use \"FAST\".")
}

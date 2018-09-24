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
      p = c(.estBetaParams(m, v), mean=m, var=v)
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

  if(is.character(init) && init == 'peaks')
  {
    # Compute KDE
    h = hist(X, breaks = seq(0, 1, 0.01), plot = FALSE)
  
    # Detect peaks
    peaks = .find_peaks(h$density, 1)
    x.peaks = (peaks * 0.01)
    
    # store only peaks above 0.1
    peaks = peaks[x.peaks > 0.1]
    x.peaks = x.peaks[x.peaks > 0.1]
    peakValues = h$density[peaks]
    
    # Cluster their x-coordinates
    clus = cnt = NULL
    if(length(x.peaks) >= K)
    {
      clus = kmeans(x.peaks, K, nstart = 100)
      cnt = as.vector(clus$centers)
    }
    else{
      clus = kmeans(x.peaks, length(x.peaks), nstart = 100)
      cnt = c(
        as.vector(clus$centers),
        runif(K - lenght(x.peaks))
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
        paste0('C',i), "a", NA, unlist(beta['a', i]),
        paste0('C',i), "b", NA, unlist(beta['a', i]),
        paste0('C',i), "Mean", NA, unlist(beta['mean', i]),
        paste0('C',i), "Variance", NA, unlist(beta['var', i])
      )
    )
  
  # Tail, sample and make a tibble
  Pareto = NULL
  
  if(tail)
  {
    shape = runif(1, min = pareto.shape$min.val, max = pareto.shape$max.val)
    scale = min(X) - 1e-9
    mv = .MeanVarPareto(shape, scale)
  
    Pareto = tibble::tribble(
        ~cluster, ~type, ~fit.value, ~init.value,
        'Tail', "Shape", NA, shape,
        'Tail', "Scale", NA, scale,
        'Tail', "Mean", NA, mv$mean,
        'Tail', "Variance", NA, mv$var
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

.stoppingCriterion = function(i, prevNLL, NLL, prevpi, pi, fit.type, epsilon, isDebug)
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
    # Compute pi's difference after ith iteration
    pi.Diff  <- abs(prevpi - pi)
    
    # Check for convergence -- all pi's have changed less than epsilon
    stopping = all(pi.Diff < epsilon)
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



#' Title
#'
#' @param x
#' @param peak.range
#' @param m.peak
#'
#' @return
#' @export
#'
#' @examples
guess_purity_from_peaks = function(x, peak.range = c(0.2, 0.5), m.peak = 3)
{
  pio::pioTit(paste("Guessing purity via peak detection in", peak.range[1], '--', peak.range[2], 'for diploid SNVs'))
  
  # Detect peaks
  x = x[x>peak.range[1]]
  x = x[x<peak.range[2]]
  h = hist(x, breaks = seq(0, 1, 0.01), plot = FALSE)
  
  peaks = dbpmm:::.find_peaks(h$density, m.peak)
  x.peaks = (peaks * 0.01)
  
  # diploid tumour and normal
  guess = 2 * x.peaks
  
  pio::pioDisp(guess)
  
  guess
}

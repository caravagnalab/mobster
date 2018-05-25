.initializer = function(X, K, init, pareto.shape = list(min.val = 0.01, max.val = 5))
{
  if(is.list(init)) {

    stopifnot(c('beta', 'shape', 'scale') %in% names(init))
    stopifnot(c('a', 'b') %in% rownames(init$beta))

    return(init)
  }

  if(is.character(init) && init == 'random')
  {
    # Sampler for Beta parameters -- reject samples where the mean/variance lead to negative
    # parameters (a and b are required to be strictly positive)
    bsamp = function(fake) {
      repeat{
        m = runif(1)
        v = runif(n = 1, min = 0.001, max = 0.25)
        p = c(.estBetaParams(m, v), mean=m, var=v)
        if(all(p > 0)) return(p)
      }
    }

    # K Beta + 1 Pareto
    beta = sapply(1:K, bsamp)
    shape = runif(1, min = pareto.shape$min.val, max = pareto.shape$max.val)
    scale = min(X) - 1e-9

    return(list(beta = beta, shape = shape, scale = scale))
  }

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
    # h = hist(X, breaks = seq(0, 1, 0.01), freq = F, plot = FALSE)
    h = hist(X, breaks = seq(0, 1, 0.01), plot = FALSE)
    # maxD = max(h$density)

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
      # print(x.peaks)
      clus = kmeans(x.peaks, length(x.peaks), nstart = 100)
      cnt = c(
        as.vector(clus$centers),
        runif(K - lenght(x.peaks))
      )
    }

    # Colors to spot them
    # col = RColorBrewer::brewer.pal(8, 'Set2')

    # Print the peaks, coloured by cluster assignment
    # df = data.frame(peaks = x.peaks, peakValues, c = clus$cluster)
    # lapply(split(df, f = df$c), function(w)
    #   points(w$peaks, w$peakValues, col = col[w$c], pch = 19))

    # A line for the mean
    # abline(v = as.vector(clus$centers), lty = 2, col = col)

    # Rectangles
    # col = ggplot2::alpha(col, 0.5)
    # for(i in 1:K){
    #   d = split(df, f = df$c)[[i]]
    #   x.left = min(d$peaks)
    #   rect(
    #     min(d$peaks), 0,
    #     max(d$peaks),
    #     maxD,
    #     col = col[i], density = 23, border = NA)
    # }

    # Random Beta
    # Sampler for Beta parameters -- reject samples where the variance leads to negative
    # parameters (a and b are required to be strictly positive)
    bsamp = function(m) {
      repeat{
        v = runif(n = 1, min = 0.001, max = 0.25)
        p = c(dbpmm:::.estBetaParams(m, v), mean=m, var=v)
        if(all(p > 0)) return(p)
      }
    }

    beta = sapply(cnt, bsamp)

    # domain = seq(0, 1, 0.01)
    # for(i in 1:K) {
    #   par = bsamp(as.vector(clus$centers)[i])
    #   lines(domain, dbeta(domain, par$alpha, par$beta), col = col[i])
    # }

    shape = runif(1, min = pareto.shape$min.val, max = pareto.shape$max.val)
    scale = min(X) - 1e-9
    # lines(domain, sads::dpareto(domain, shape, scale), col = 'red')

    return(list(beta = beta, shape = shape, scale = scale))
  }

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



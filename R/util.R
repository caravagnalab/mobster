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
  return(params = list(alpha = alpha, beta = beta))
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

# Setup parallel for multi-core executions
.setup_parallel = function(cores.ratio)
{
  require(parallel)
  require(doParallel)
  require(crayon)

  # set the number of cores to be used in the parallelization
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1) cores = 1

  cat(cyan('\n\tRegistering to use', cores, 'cores out of', detectCores(), 'via \"parallel\" ...'))
  cl = makeCluster(cores)
  registerDoParallel(cl)
  cat(bgGreen(" OK.\n"))

  return(cl)
}

# Stop parallel for multi-core executions
.stop_parallel =  function(cl)
{
  cat(cyan('\n\tStopping \"parallel\" clusters: '))

  stopCluster(cl)
  cat(bgGreen("OK.\n"))
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


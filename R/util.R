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


.params_Beta = function(x, init = FALSE)
{
  if(init) x$Clusters$fit.value = x$Clusters$init.value
    
  x$Clusters %>%
    filter(cluster != 'Tail', type == 'a' | type == 'b') %>%
    select(-init.value) %>%
    spread(key = type, value = fit.value)   
}

.params_Pareto = function(x, init = FALSE)
{
  if(init) x$Clusters$fit.value = x$Clusters$init.value
  
  x$Clusters %>%
    filter(cluster == 'Tail', type == 'Shape' | type == 'Scale') %>%
    select(-init.value) %>%
    spread(key = type, value = fit.value)   
}

.params_Pi = function(x, init = FALSE)
{
  if(init) x$Clusters$fit.value = x$Clusters$init.value
  
  v = x$Clusters %>%
    filter(type == 'Mixing proportion') %>%
    select(-init.value) %>%
    spread(key = type, value = fit.value)   
  
  pi = v$`Mixing proportion`
  names(pi) = v$cluster
  
  ord.pi = c(pi['Tail'], pi[names(pi) != 'Tail'])
  
  ord.pi
}

.set_params_Beta = function(fit, a, b)
{
  names.BetaC = names(a)
  if(any(is.null(names.BetaC))) stop("params -- named vector required?")
  
  fit$Clusters[
  fit$Clusters$cluster %in% names.BetaC & fit$Clusters$type == 'a',
  'fit.value'
  ] =  a

  fit$Clusters[
    fit$Clusters$cluster %in% names.BetaC & fit$Clusters$type == 'b',
    'fit.value'
    ] =  b
  
  for (s in names.BetaC)
  {
    mv = .MeanVarBeta(a[s], b[s])
    
    fit$Clusters[
      fit$Clusters$cluster == s & fit$Clusters$type == 'Mean',
      'fit.value'
      ] =  mv$mean
    
    fit$Clusters[
      fit$Clusters$cluster == s & fit$Clusters$type == 'Variance',
      'fit.value'
      ] =  mv$var
  }
  
  fit
}

.set_params_Pareto = function(fit, shape, scale)
{
  if(!fit$fit.tail) return(fit)

  mv = .MeanVarPareto(shape, scale)
  
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

.set_params_Pi = function(fit, pi)
{
  fit$Clusters[
    fit$Clusters$type == 'Mixing proportion',
    'fit.value'
    ] =  pi
  
  fit
}

# 
# .extract.DP.fit = function(mfit, with.summary = FALSE)
# {
#   fit = NULL
# 
#   # return the posterior density
#   # fit$densp.m = mfit$densp.m
# 
#   # then the order statistics as well
#   rsave = mfit$save.state$randsave
#   maxK = max(apply(rsave, 1, function(x)length(unique(x))))
#   K_fit = round(mfit$coefficients['ncluster'])
#   alpha_fit = mfit$coefficients['alpha']
# 
#   srsave = spropsave = matrix(NA, nrow(rsave), maxK)
#   for(i in 1:nrow(srsave))
#   {
#     v = sort(unique(rsave[i, ]), decreasing = T)
#     srsave[i, 1:length(v)] = v
#     spropsave[i, 1:length(v)] = (table(rsave[i, ])/ncol(rsave))[as.character(v)]
#   }
# 
#   p = matrixStats::colMedians(srsave, na.rm = T)[1:K_fit]
#   pi = matrixStats::colMedians(spropsave, na.rm = T)[1:K_fit]
# 
#   fit$theta_k = paste(p, collapse = ':')
#   fit$pi_k = paste(pi/sum(pi), collapse = ':')
#   # fit$alpha = alpha_fit
#   # fit$prior = mfit$prior
#   fit$K_fit = K_fit
# 
#   if(with.summary) fit$summary = summary(mfit)
#   data.frame(fit)
# }
# 
# 
# .extract.Bmix.fit = function(bmix.fit)
# {
#   fit = NULL
# 
#   tk = ttk = NULL
#   if(bmix.fit$K[1] > 0)
#   {
#    tk = paste(bmix.fit$B.params, collapse = ':')
#   }
# 
#   if(bmix.fit$K[2] > 0)
#   {
#     ttk = paste(bmix.fit$BB.params[1, ], collapse = ':')
#   }
# 
#   fit$theta_k = paste(tk, ttk, sep = '')
#   fit$pi_k = paste(bmix.fit$pi, collapse = ':')
#   fit$K_fit = sum(bmix.fit$K)
# 
#   data.frame(fit)
# }
# 
# .extract.vbdbmm.fit = function(mfit, pi.cutoff = 0.01, with.summary = FALSE)
# {
#   fit = NULL
# 
#   which.selection = mfit$pi_k > pi.cutoff
#   # which.selection = TRUE
# 
#   fit$theta_k = paste(mfit$theta_k[which.selection], collapse = ':')
#   fit$pi_k = paste(mfit$pi_k[which.selection], collapse = ':')
#   fit$K_fit = length(mfit$pi_k[which.selection])
# 
#   if(with.summary) fit$summary = summary(mfit)
# 
#   cat('Extracted mixtures components with minimum pi', pi.cutoff, '\n')
#   data.frame(fit)
# }



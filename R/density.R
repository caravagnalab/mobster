

#' Density function for a Dirichlet Beta-Pareto Mixture Model
#'
#' @param x an object of class dbpmm
#' @param data data to compute the (log)-likelihood, if this is null the data stored inside the object is used.
#' @param components a vector specifying which components should be used to compute the density (this is used in the implementation of the
#' fitting algorithms for this mixture -- most users would not use this option).
#' @param init if TRUE, it uses the initial values for model fit, otherwise final.
#' @param log log-likelihood or just likelihood.
#'
#' @return The density for the data. Notice, this is NOT the negative LL returned during fit (which is this value with different sign).
#' @export
#' @import sads
#'
#' @examples Will make some..
ddbpmm = function(x, data = NULL, components = 1:x$K, 
                  a = NULL, b = NULL, pi = NULL, scale = NULL, shape = NULL,
                  log = TRUE) {

  stopifnot(length(components) > 0)

  # Get parameters and data if these are not passed explicitly
  if(all(is.null(data))) data = x$data$VAF

  if(is.null(scale)) scale = .params_Pareto(x)$Scale
  if(is.null(shape)) shape = .params_Pareto(x)$Shape

  if(all(is.null(a))) a = .params_Beta(x)$a
  if(all(is.null(b))) b = .params_Beta(x)$b
  
  if(all(is.null(pi))) pi = .params_Pi(x)
  
  # Mask what we do not compute
  mask = (1:x$K) %in% components

  log_pi = log(pi)
  
  logLik = 0
  
  # Tail
  if(mask[1])
  {
    if(x$fit.tail)
      logLik = log_pi[1] + sads::dpareto(data, 
              shape = shape, 
              scale = scale, log = TRUE) * mask[1]
    else logLik = rep(-Inf, length(data)) # Otherwise -Inf is the min achievable
  }
  
  log_pi = log_pi[2:x$K]
  
  for(k in 1:length(log_pi)) 
    if(mask[k + 1])
      logLik = logLik +
        log_pi[k] + dbeta(data,
                         shape1 = a[k], shape2 = b[k], log = TRUE) 
      
  if(!log) logLik = exp(logLik)

  logLik
}


template_density = function(x, x.axis = seq(0, 1, 0.01), init = FALSE, binwidth = 0.01, reduce = FALSE)
{
  labels = names(.params_Pi(x))
  
  values = lapply(
    1:x$K, # Component wise
    function(w)
      data.frame(
        cluster = labels[w],
        x = x.axis,
        y = ddbpmm(x, data = x.axis, components = w, 
                   # init = init, 
                   log = FALSE)))
  
  names(values) = labels
  
  # Scale density wrt binwidth
  values = lapply(values, function(w) {w$y = w$y * binwidth; w}) 
  
  # Reduce if required
  if(reduce) values = Reduce(rbind, values)
  values
}



# UNUSED FUNCTION: numerical nlogLik for a Pareto component (the actual solution is analytical)
# a, b parameters of the mixture, k the component ID, and x data
# Creates a functional that depends only on a and b
# .NLLParetoMix = function(x, z_nk, pi, scale)
# {
#   f = function(shape) # NLL
#   {
#     # z[x,k] * { log pi[k] + log Pareto(x | shape, scale) }
#     -1 * z_nk[, 1] %*% (log(pi[1]) + sads::dpareto(x, shape = shape, scale = scale, log = TRUE))
#   }
# 
#   return(f)
# }



#' Generate a random sample from a DBPMM model.
#'
#' @param x an object with computed from MOBSTER
#' @param n number of samples required
#' @param a if passed, `a` for the Beta components
#' @param b if passed, `a` for the Beta components
#' @param pi if passed, `pi` for the mixing proportions
#' @param shape if passed, `shape` of the tail
#' @param scale if passed, `scale` of the tail
#' @param tail.cutoff reject Pareto Type I samples above this value (ensure that will be returned
#' the number of samples required anyway).
#'
#' @return n samples from the mixture
#' @export
#'
#' @examples
#' a = b = 50
#' names(a) = names(b) = "C1"
#' 
#' pi = c('Tail' = .6, 'C1' = .4)
#' 
#' v = rdbpmm(x = NULL, n = 1000, a = a, b = b, pi = pi, shape = 2, scale = 0.05)
#' hist(v, breaks = seq(0, 1, 0.01))
rdbpmm = function(x, 
                  a = NULL,
                  b = NULL,
                  pi = NULL,
                  shape = NULL,
                  scale = NULL,
                  n = 1,
                  tail.cutoff = 1) 
{
  # Get parameters if these are not passed explicitly
  if(all(is.null(scale)) | all(is.null(shape))) 
  {
    scale = .params_Pareto(x)$Scale
    shape = .params_Pareto(x)$Shape
  }
  
  if(all(is.null(a)) | all(is.null(b))) 
  {
    Betas = .params_Beta(x)
    a = Betas$a
    b = Betas$b
    names(a) = names(b) = Betas$cluster
  }
  
  if(all(is.null(pi))) pi = .params_Pi(x)
  
  ############################ Sample mixing prop.
  pi.samples = sample(names(pi), size = n, prob = pi, replace = TRUE)
  
  ############################ Tail samples
  n.tail = table(pi.samples)['Tail']
  sample.tails = NULL
  
  if(!is.na(n.tail) & n.tail > 0) 
  {
    N = n.tail
    
    # as the model is an approximation to the VAF we reject values > 1.
    repeat{
      # Get `N` Pareto Type I samples
      new.sample.tails = sads::rpareto(
        N,
        shape = shape,
        scale = scale)
      
      # Store only values <= 1, which might be less than N
      new.sample.tails = new.sample.tails[new.sample.tails <= tail.cutoff]
      
      sample.tails = c(sample.tails, new.sample.tails)
      
      # Stop if we accumulated `n.tail` samples
      if(length(sample.tails) == n.tail) break;
      
      # Otherwise repeat drawing N-length(sample.tails)
      N = N - length(new.sample.tails)
    }
  }

  ############################ Beta samples
  n.Betas = table(pi.samples)[names(table(pi.samples)) != 'Tail']
  
  sample.Betas = sapply(
    names(n.Betas),
    function(w)
    {
      rbeta(n.Betas[w], shape1 = a[w], shape2 = b[w])
    }
  )
  
  return(c(unlist(sample.Betas), unlist(sample.tails)))
}  
  
 # USED FOR NUMERICAL SOLUTION OF THE MLE ESTIMATORES FOR EACH BETA COMPONENT
# a, b parameters of the mixture, k the component ID, and x data
# Creates a functional that depends only on a and b
# .NLLBetaMix = function(k, x, z_nk, pi)
# {
#   f = function(a,b) # NLL
#   {
#     # z[x,k] * { log pi[k] + log beta(x | a,b) }
#     -1 * z_nk[, k] %*% (log(pi[k]) + dbeta(x, shape1 = a, shape2 = b, log = TRUE))
#   }
# 
#   return(f)
# }

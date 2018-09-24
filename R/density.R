

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

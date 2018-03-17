

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
ddbpmm = function(x, data = NULL, components = 1:x$K, init = FALSE, log = TRUE) {

  stopifnot(length(components) > 0)

  if(is.null(data)) data = x$X
  if(init){
    if(!is.na(x$shape) & !is.na(x$scale)) {
      x$shape = x$tail$shape.init
      x$scale = x$tail$scale.init
    }

    x$pi = x$pi.init
    x$a = x$beta['a.init', ]
    x$b = x$beta['b.init', ]
  }

  tail = 0
  if(1 %in% components) # If one wants the tail
  {
    # If the parameters are defined, one gets the logLikelihood
    if(!is.na(x$shape) & !is.na(x$scale))
    {
        tail = base::log(x$pi[1]) + sads::dpareto(data, shape = x$shape, scale = x$scale, log = TRUE)
    }
    else tail = rep(-Inf, length(data)) # Otherwise -Inf is the min achievable
  }

  betas = 0
  for (k in 2:x$K) {
    if(k %in% components) {
      betas = betas + base::log(x$pi[k]) + dbeta(data, shape1 = x$a[k-1], shape2 = x$b[k-1], log = TRUE)
    }
  }

  logLik = betas + tail

  if(!log) logLik = exp(logLik)

  logLik
}



# UNUSED FUNCTION: numerical nlogLik for a Pareto component (the actual solution is analytical)
# a, b parameters of the mixture, k the component ID, and x data
# Creates a functional that depends only on a and b
.NLLParetoMix = function(x, z_nk, pi, scale)
{
  f = function(shape) # NLL
  {
    # z[x,k] * { log pi[k] + log Pareto(x | shape, scale) }
    -1 * z_nk[, 1] %*% (log(pi[1]) + sads::dpareto(x, shape = shape, scale = scale, log = TRUE))
  }

  return(f)
}


# USED FOR NUMERICAL SOLUTION OF THE MLE ESTIMATORES FOR EACH BETA COMPONENT
# a, b parameters of the mixture, k the component ID, and x data
# Creates a functional that depends only on a and b
.NLLBetaMix = function(k, x, z_nk, pi)
{
  f = function(a,b) # NLL
  {
    # z[x,k] * { log pi[k] + log beta(x | a,b) }
    -1 * z_nk[, k] %*% (log(pi[k]) + dbeta(x, shape1 = a, shape2 = b, log = TRUE))
  }

  return(f)
}

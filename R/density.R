#' Density function for MOBSTER.
#'
#' @description The density function for a MOBSTER mixture. The function
#' allows to compute the density from the data and parameters stored
#' inside an object of class \code{dbpmm}, or from
#' a custom set of parameters and data that can be passed as input. The
#' function can also compute the density only for a subset of components
#' of the mixture.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param data Data to compute the (log)-likelihood, if this is null the data stored
#' inside \code{x$data} is used.
#' @param components A vector specifying which components should be used to compute
#' the density. Position 1 is reserved for the tail (regardless that is fit or not
#' to the data), the other positions refer to Beta components.
#' @param log Boolean value to select the log-likelihood, or the likelihood.
#' @param a Vector of parameters for the Beta components; if this is null values stored
#' inside \code{x} are used.
#' @param b Vector of parameters for the Beta components; if this is null values stored
#' inside \code{x} are used.
#' @param pi Mixing proportions for the mixture; if this is null values stored
#' inside \code{x} are used.
#' @param scale Scale of the power law; if this is null values stored
#' inside \code{x} are used.
#' @param shape Shape of the power law; if this is null values stored
#' inside \code{x} are used.
#'
#' @return The density for the data. Notice that to compute the negative log-likelihood
#' used during fit one needs to use \code{log = TRUE} and change sign.
#' 
#' @export
#' 
#' @importFrom  sads dpareto
#'
#' @examples 
#' library(ggplot2)
#' data('fit_example', package = 'mobster')
#' 
#' # Use the full mixture, and its internal data
#' ddbpmm(fit_example$best)
#' 
#' # Use only some of the mixture components, and pass some data
#' ddbpmm(fit_example$best, data = .4, components = 1)
#' 
#' # An internal function to get f(x) with x the [0,1] range.
#' ggplot(mobster:::template_density(fit_example$best, reduce = TRUE),
#'        aes(x = x, y = y, color = cluster)) + geom_line()
ddbpmm = function(x,
                  data = NULL,
                  components = 1:x$K,
                  a = NULL,
                  b = NULL,
                  pi = NULL,
                  scale = NULL,
                  shape = NULL,
                  log = TRUE) 
{
  is_mobster_fit(x)
  stopifnot(length(components) > 0)
  
  # Get parameters and data if these are not passed explicitly
  if (all(is.null(data)))
    data = x$data$VAF
  
  if (is.null(scale))
    scale = mobster:::.params_Pareto(x)$Scale
  if (is.null(shape))
    shape = .params_Pareto(x)$Shape
  
  if (all(is.null(a)))
    a = .params_Beta(x)$a
  if (all(is.null(b)))
    b = .params_Beta(x)$b
  
  if (all(is.null(pi)))
    pi = .params_Pi(x)
  
  # Mask what we do not compute
  mask = (1:x$K) %in% components
  
  log_pi = log(pi)
  
  logLik = 0
  
  # Tail
  if (mask[1])
  {
    if (x$fit.tail)
      logLik = log_pi[1] + sads::dpareto(data,
                                         shape = shape,
                                         scale = scale,
                                         log = TRUE) * mask[1]
    else
      logLik = rep(-Inf, length(data)) # Otherwise -Inf is the min achievable
  }
  
  # Beta components
  log_pi = log_pi[2:x$K]
  
  for (k in 1:length(log_pi))
    if (mask[k + 1])
      logLik = logLik +
    log_pi[k] + dbeta(data,
                      shape1 = a[k],
                      shape2 = b[k],
                      log = TRUE)
  
  if (!log)
    logLik = exp(logLik)
  
  logLik
}


template_density = function(x,
                            x.axis = seq(0, 1, 0.01),
                            init = FALSE,
                            binwidth = 0.01,
                            reduce = FALSE)
{
  labels = names(mobster:::.params_Pi(x))
  
  values = lapply(1:x$K,
                  # Component wise
                  function(w)
                    data.frame(
                      cluster = labels[w],
                      x = x.axis,
                      y = mobster::ddbpmm(
                        x,
                        data = x.axis,
                        components = w,
                        # init = init,
                        log = FALSE
                      )
                    ))
  
  names(values) = labels
  
  # remove tail if not required, which is entry #1
  if (!x$fit.tail)
    values = values[2:x$K]
  
  # Scale density wrt binwidth
  values = lapply(values, function(w) {
    w$y = w$y * binwidth
    w
  })
  
  # Reduce if required
  if (reduce)
    values = Reduce(rbind, values)
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


#' Generate a random sample from a MOBSTER model.
#' 
#' @description This function can be used to generate data from a MOBSTER
#' model, or from a custom mixture. The principle is the same of function
#' \code{\link{ddbpmm}}.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param n Number of samples required (i.e., points).
#' @param a Vector of parameters for the Beta components; if this is null values stored
#' inside \code{x} are used.
#' @param b Vector of parameters for the Beta components; if this is null values stored
#' inside \code{x} are used.
#' @param pi Mixing proportions for the mixture; if this is null values stored
#' inside \code{x} are used.
#' @param scale Scale of the power law; if this is null values stored
#' inside \code{x} are used.
#' @param shape Shape of the power law; if this is null values stored
#' inside \code{x} are used.
#' @param tail.cutoff Because the Pareto Type I power law has support
#' over all the positive real line, tail values above \code{tail.cutoff}
#' are rejected (and re-sampled). 
#'
#' @return n samples from the mixture
#' 
#' 
#' @export
#'
#' @examples
#' library(ggplot2)
#' # 1 Beta component at 0.5 mean (symmetrical Beta)
#' a = b = 50
#' names(a) = names(b) = "C1"
#'
#' # 60% tail mutations
#' pi = c('Tail' = .6, 'C1' = .4)
#'
#' # Sample
#' v = data.frame(x = rdbpmm(
#'    x = NULL,
#'    n = 1000,
#'    a = a,
#'    b = b,
#'    pi = pi,
#'    shape = 2,
#'    scale = 0.05
#'  ))
#' ggplot(v, aes(x)) + geom_histogram(binwidth = 0.01)
#' 
#' # Or use the parameters of a model available
#' data('fit_example', package = 'mobster')
#' 
#' v = data.frame(x = rdbpmm(x = fit_example$best, n = 1000))
#' ggplot(v, aes(x)) + geom_histogram(binwidth = 0.01)
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
  if (all(is.null(scale)) | all(is.null(shape)))
  {
    scale = .params_Pareto(x)$Scale
    shape = .params_Pareto(x)$Shape
  }
  
  if (all(is.null(a)) | all(is.null(b)))
  {
    Betas = .params_Beta(x)
    a = Betas$a
    b = Betas$b
    names(a) = names(b) = Betas$cluster
  }
  
  if (all(is.null(pi)))
    pi = .params_Pi(x)
  
  ############################ Sample mixing prop.
  pi.samples = sample(names(pi),
                      size = n,
                      prob = pi,
                      replace = TRUE)
  
  ############################ Tail samples
  n.tail = table(pi.samples)['Tail']
  sample.tails = NULL
  
  if (!is.na(n.tail) & n.tail > 0)
  {
    N = n.tail
    
    # as the model is an approximation to the VAF we reject values > 1.
    repeat {
      # Get `N` Pareto Type I samples
      new.sample.tails = sads::rpareto(N,
                                       shape = shape,
                                       scale = scale)
      
      # Store only values <= 1, which might be less than N
      new.sample.tails = new.sample.tails[new.sample.tails <= tail.cutoff]
      
      sample.tails = c(sample.tails, new.sample.tails)
      
      # Stop if we accumulated `n.tail` samples
      if (length(sample.tails) == n.tail)
        break
      
      
      # Otherwise repeat drawing N-length(sample.tails)
      N = N - length(new.sample.tails)
    }
  }
  
  ############################ Beta samples
  n.Betas = table(pi.samples)[names(table(pi.samples)) != 'Tail']
  
  sample.Betas = sapply(names(n.Betas),
                        function(w)
                        {
                          rbeta(n.Betas[w], shape1 = a[w], shape2 = b[w])
                        })
  
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

# Compute parametric bootstrap replicate
.parametric_bootstrap_resamples = function(x, n = 100)
{
  stopifnot(n >= 1)

  data.size = nrow(x$data)

  lapply(1:n,
         function(w) {
           resample = rdbpmm(x, n = data.size)
           data.frame(id = 1:data.size, VAF = resample)
         })
}

# Compute nonparametric bootstrap replicate
.nonparametric_bootstrap_resamples = function(x, n = 100)
{
  stopifnot(n >= 1)

  data.size = nrow(x$data)

  lapply(1:n,
         function(w) {
           ids = sample(1:data.size, data.size, replace = TRUE)
           data.frame(
             id = 1:data.size,
             VAF = x$data$VAF[ids],
             original.id = ids
           )
         })
}


#' Bootstrap a MOBSTER fit.
#' 
#' @description Parmaetric and nonparametric implementation of
#' bootstrap estimates for MOBSTER fits. This computation is parallel
#' and uses  \code{?easypar}.
#' 
#' @param x An object of class \code{"dbpmm"}.
#' @param cores.ratio  Ratio of cores to use for the parallel; see \code{?easypar}.
#' @param bootstrap Type of boostrap: \code{"parametric"} or \code{"nonparametric"}
#' @param n.resamples Number of boostrap resamples.
#' @param cache Cache for the computation; see \code{?easypar}.
#' @param ... fit parameters for \code{mobster_fit}
#'
#' @return Data from the fits, resamples and a plottable figure.
#'
#' @import easypar
#'
#' @export
#'
#' @examples
#' TODO
mobster_bootstrap = function(x,
                             n.resamples = 100,
                             bootstrap = 'parametric',
                             cores.ratio = 0.8,
                             cache = paste0(bootstrap, '_cache.rds'),
                             ...)
{
  pio::pioHdr(
    "MOBSTER bootstrap",
    toPrint = c(`Resamples` = paste(n.resamples),
                `Type of bootstrap` = bootstrap),
    prefix = '\t-'
  )

  stopifnot(inherits(x, "dbpmm"))
  stopifnot(bootstrap %in% c('parametric', 'nonparametric'))

  pio::pioTit(paste0("Bootstrapping for this MOBSTER model"))
  print(x)

  pio::pioTit(paste0("Creating ", bootstrap, " bootstrap resamples"))

  resamples = NULL

  # Get resamples -- n datasets with same size of x
  if (bootstrap == 'parametric')
    resamples = .parametric_bootstrap_resamples(x, n.resamples)

  if (bootstrap == 'nonparametric')
    resamples = .nonparametric_bootstrap_resamples(x, n.resamples)

  resamples = lapply(resamples, list)

  save(resamples, file = 'resamples.RData')

  pio::pioTit("Running fits (might take some time)")

  # easypar
  fits = easypar::run(
    FUN = function(w) {
      mobster_fit(x = w,
                  parallel = FALSE,
                  seed = NULL,
                  ...)$best
    },
    PARAMS = resamples,
    packages = c("crayon", "mobster"),
    export = ls(globalenv(), all.names = TRUE),
    cores.ratio = cores.ratio,
    parallel = TRUE,
    cache = cache
  )

  
  return(list(resamples = resamples, fits = fits))
}



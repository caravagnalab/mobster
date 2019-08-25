
#' Bootstrap a MOBSTER fit.
#' 
#' @description Parametric and non-parametric implementation of
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
#' data('fit_example', package = 'mobster')
#' 
#' # Just 5 resamples of a nonparametric bootstrap run
#' mb = mobster_bootstrap(fit_example$best, n.resamples = 5)
#' 
#' # The resample data is wrapped as lists of lists
#' lapply(mb$resamples, function(x) x[[1]] %>% as_tibble)
#' 
#' # The best fits are returned
#' mb$fits
mobster_bootstrap = function(x,
                             n.resamples = 100,
                             bootstrap = 'nonparametric',
                             cores.ratio = 0.8,
                             cache = paste0(bootstrap, '_cache.rds'),
                             save_data = NULL,
                             ...)
{
  pio::pioHdr(
    paste0("MOBSTER bootstrap ~ ", n.resamples, ' resamples from ',
           bootstrap, ' bootstrap')
  )
  
  stopifnot(inherits(x, "dbpmm"))
  stopifnot(bootstrap %in% c('parametric', 'nonparametric'))
  
  pio::pioTit(paste0("Bootstrapping for this MOBSTER model"))
  print(x)
  
  pio::pioTit(paste0("Creating ", bootstrap, " bootstrap resamples"))
  
  resamples = NULL
  
  # Get resamples -- n datasets with same size of x
  if (bootstrap == 'parametric')
    resamples = mobster:::.parametric_bootstrap_resamples(x, n.resamples)
  
  if (bootstrap == 'nonparametric')
    resamples = mobster:::.nonparametric_bootstrap_resamples(x, n.resamples)
  
  resamples = lapply(resamples, list)
  
  # Save data if required
  if(!is.null(save_data) & is.character(save_data))
  {
    pio::pioStr("Resamples saved to file ", paste0(save_data, '.RData'), suffix = '\n')
    
    save(resamples, file = paste0(save_data, '.RData'))
  }
  
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
  
  
  return(list(resamples = resamples, fits = fits, bootstrap = bootstrap))
}

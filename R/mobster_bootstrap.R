
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
#' @export
#'
#' @examples
#' # Random small dataset
#' dataset = random_dataset(N = 200, seed = 123, Beta_variance_scaling = 100)
#' x = mobster_fit(dataset$data, auto_setup = 'FAST')
#' 
#' # Just 5 resamples of a nonparametric bootstrap run, disabling the parallel engine
#' options(easypar.parallel = FALSE)
#' boot_results = mobster_bootstrap(x$best, n.resamples = 5, auto_setup = 'FAST')
#' 
#' # The resample data is available in a list
#' print(boot_results$resamples[[1]])
#' 
#' # The best fits are returned
#' print(boot_results$fits)
mobster_bootstrap = function(x,
                             n.resamples = 100,
                             bootstrap = 'nonparametric',
                             cores.ratio = 0.8,
                             cache = NULL,
                             save_data = NULL,
                             ...)
{
  pio::pioHdr(
    paste0("MOBSTER bootstrap ~ ", n.resamples, ' resamples from ',
           bootstrap, ' bootstrap')
  )
  cat('\n')
  
  is_mobster_fit(x)
  stopifnot(bootstrap %in% c('parametric', 'nonparametric'))
  
  # pio::pioTit(paste0("Bootstrapping for this MOBSTER model"))
  print(x)
  cat('\n')
  
  cli::cli_process_start(paste0("Creating ", bootstrap, " bootstrap resamples"))
  
  resamples = NULL
  
  # Get resamples -- n datasets with same size of x
  if (bootstrap == 'parametric')
    resamples = mobster:::.parametric_bootstrap_resamples(x, n.resamples)
  
  if (bootstrap == 'nonparametric')
    resamples = mobster:::.nonparametric_bootstrap_resamples(x, n.resamples)
  
  resamples = lapply(resamples, list)
  
  cli::cli_process_done()
  
  # Save data if required
  if(!is.null(save_data) & is.character(save_data))
  {
    mobster:::m_inf("Resamples will be saved to file {.field {save_data}.RData}")
    
    save(resamples, file = paste0(save_data, '.RData'))
  }
  
  cat('\n')
  cli::cli_rule("Running fits", right = "Might take some time ... ")

  # easypar
  fits = easypar::run(
    FUN = function(w) {
      mobster_fit(x = w,
                  parallel = FALSE,
                  seed = NULL,
                  ...)$best
    },
    PARAMS = resamples,
    packages = c("dplyr", "tidyr", "mobster"),
    export = ls(globalenv(), all.names = TRUE),
    cores.ratio = cores.ratio,
    parallel = TRUE,
    cache = cache
  )
  
  # Check for errors
  # Polish errors if any
  nerrs = easypar::numErrors(fits)
  
  if(nerrs == n.resamples) {
    
    lapply(fits, function(w) print(w$message))
    
    stop("All task returned errors, no fit available, raising error.")
  }
  
  errors = NULL
  if(nerrs > 0) {
    cat(crayon::red(nerrs, "tasks raised an error, filtering them out.\n"))
    cat(crayon::red(n.resamples - nerrs, "bootstrap(s) available.\n"))
    
    errs = sapply(fits, function(w) inherits(w, 'simpleError') | inherits(w, 'try-error'))
    errors = fits[errs]
  }
  
  fits = easypar::filterErrors(fits)
  
  return(list(resamples = resamples, fits = fits, bootstrap = bootstrap, errors = errors))
}

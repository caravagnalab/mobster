# lapply(.parametric_bootstrap_resamples(x, n = 2), function(w) hist(w$VAF, breaks = seq(0, 1, 0.01)))
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


#' Parallel implementation of parametric/ nonparametric 
#' bootstrap for MOBSTER.
#'
#' @param x a fit by MOBSTER
#' @param n number of parametric resamples (drawn via \code{rdbpmm})
#' @param cores.ratio  ratio of cores to use
#' @param ... fit parameters for \code{mobster_fit}
#' @param palette plot palette used for Beta clusters
#' @param alpha transparency
#' @param tail.color for the tail's boxplot
#' @param bootstrap "parametric" or "nonparametric"
#'
#' @return Data from the fits, resamples and a plottable figure.
#' 
#' @import easypar
#' 
#' @export
#'
#' @examples
mobster_bootstrap = function(x,
                             n.resamples = 100,
                             bootstrap = 'parametric',
                             palette = 'Set1',
                             alpha = 1,
                             tail.color = 'darkgray',
                             cores.ratio = 0.8,
                             cache = NULL,
                             ...)
{
  pio::pioHdr(
    "MOBSTER bootstrap",
    toPrint = c(`Resamples` = paste(n.resamples),
                `Type of bootstrap` = bootstrap),
    prefix = '\t-'
  )
  
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
  
  # # Setup clusters for parallel computing
  # cl = .setup_parallel(cores.ratio = cores.ratio)
  # 
  # 
  # 
  # # perform parallel inferences with custom parameters
  # # -- disable internal parallelism
  # # -- randomize seed
  # fits = foreach(
  #   num = 1:n.resamples,
  #   .packages = c("crayon", "mobster"),
  #   .export = ls(globalenv(), all.names = TRUE)
  # ) %dopar%
  # {
  #   # best fit from resample -- with control for errors
  #   tryCatch({
  #     
  #     fit = mobster_fit(x = resamples[[num]],
  #               parallel = FALSE,
  #               seed = NULL,
  #               ...)$best
  #     
  #     # Incremental saves of the computation so that we get some results as 
  #     # soon as they are computed
  #     if(!is.null(incremental)) 
  #     {
  #       pio::pioStr("CACHING RESULTS : ", incremental)
  #       
  #       obj = NULL
  #       if(file.exists(incremental)) obj = readRDS(incremental)
  # 
  #       obj = append(obj, list(fit))
  #       saveRDS(obj, file = incremental)
  #     }
  #     
  #     fit
  #     
  #   }, error = function(e) NULL)
  # }
  # 
  # .stop_parallel(cl)
  
  return(list(resamples = resamples, fits = fits))
}



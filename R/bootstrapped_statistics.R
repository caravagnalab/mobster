#' Compute boostrap statistics from a bootstrap run.
#' 
#' @description From a bootstrap run several statistics are computed.
#' These are the model frequency, the bootstrap statistics of each parameter
#' and, if the bootstrap is nonparametric, also the co-clustering probability
#' of the input points. All the bootstrapped values are also returned. For the
#' bootstrapped statistics the confidence intervals are returned for a 
#' required  alpha value.
#'
#' @param x A MOBSTER fit.
#' @param bootstrap_results The results of running function \code{mobster_bootstrap}
#' @param bootstrap Type of boostrap: \code{"parametric"} or \code{"nonparametric"}
#' @param alpha Alpha-level empirical CIs for the bootstrap distribution. 
#'
#' @return A list with the bootstrapped values, the model frequency, the bootstraped statistics and
#' the co-clustering probability - non-null only for nonparametric boostrap.
#' 
#' @importFrom cli cli_process_start cli_process_done make_spinner
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
#' boot_stats = bootstrapped_statistics(x$best, boot_results)
#' print(boot_stats)
bootstrapped_statistics = function(x, bootstrap_results, bootstrap = 'nonparametric', alpha = 0.025)
{
  is_mobster_fit(x)
  stopifnot(bootstrap %in% c('parametric', 'nonparametric'))
  
  fit = x
  resamples = bootstrap_results$resamples
  bootstrap.fits = bootstrap_results$fits
  
  n = length(bootstrap.fits)
  
  mobster:::m_inf(paste("Bootstrap observations n =", n))
  
  # cli::cli_process_start(paste("Computing model frequency"))
  cli::cli_h1("Computing model frequency")
  
  res = NULL
  
  # Overall statistic: "model frequency"
  models.tab = lapply(seq(bootstrap.fits),
                      function(w)
                        data.frame(
                          `resample` = w,
                          `tail` = bootstrap.fits[[w]]$fit.tail,
                          `K` =  bootstrap.fits[[w]]$Kbeta,
                          `Model` = paste0(
                            'K = ',
                            bootstrap.fits[[w]]$Kbeta,
                            ifelse(bootstrap.fits[[w]]$fit.tail,
                                   " with tail",
                                   " without tail")
                          ),
                          stringsAsFactors = FALSE
                        ))
  
  models.tab = Reduce(rbind, models.tab)
  
  model_fr = data.frame(table(models.tab$Model)/n, stringsAsFactors = F)
  model.frequency = model_fr %>% as_tibble %>% dplyr::arrange(dplyr::desc(Freq))
  colnames(model.frequency) = c("Model", "Frequency")
  
  this.model = paste0(
    'K = ',
    fit$Kbeta,
    ifelse(fit$fit.tail,
           " with tail",
           " without tail")
  )
  
  model.frequency$fit.model = FALSE
  model.frequency$fit.model[model.frequency$Model == this.model] = TRUE
  
  # Parameter statistics : alpha-level empirical CIs for the bootstrap distribution
  bootstrap.values = lapply(seq(bootstrap.fits),
                            function(w)
                              cbind(bootstrap.fits[[w]]$Clusters, `resample` = w))
  
  bootstrap.values = tibble::as_tibble(Reduce(rbind, bootstrap.values)) %>%
    dplyr::rename(statistics = type) 
  
  # Stats table
  fit$Clusters = fit$Clusters %>% rename(statistics = type)
  
  stats = bootstrap.values %>%
    dplyr::group_by(cluster, statistics) %>%
    dplyr::summarise(
      min = min(fit.value),
      lower_quantile = quantile(fit.value, alpha, na.rm = TRUE),
      higher_quantile = quantile(fit.value, 1-alpha, na.rm = TRUE),
      max = max(fit.value)
    ) %>%
    dplyr::left_join(fit$Clusters, by = c('cluster', 'statistics')) %>%
    dplyr::ungroup()
  
  pio::pioDisp(model.frequency)
  # cli::cli_process_done()
  
  cli::cli_h1("Computing confidence Intervals (CI) for empirical quantiles")
  
  
  # cli::cli_process_start("Confidence Intervals (CI) for empirical quantiles")
  # pio::pioDisp(stats)
  
  pio::pioStr("\nMixing proportions", "\n")
  print(stats %>%
          dplyr::filter(statistics == 'Mixing proportion'))
  
  pio::pioStr("\nTail shape/ scale", "\n")
  print(stats %>%
          dplyr::filter(cluster == 'Tail' & statistics %in% c('Shape', 'Scale')))
  
  pio::pioStr("\nBeta peaks", "\n")
  print(stats %>%
          dplyr::filter(statistics %in% c('Mean', 'Variance') & cluster != 'Tail'))
  
  # cli::cli_process_done()
  
  
  
  # Co-clustering only for nonparametric bootstrap
  co_clustering_matrix = co_clustering_labels = NULL
  if(bootstrap == 'nonparametric') {
    
    cli::cli_h1("Computing co-clustering probability for nonparametric bootstrap")
    
    # cli::cli_process_start("Co-clustering probability from nonparametric bootstrap")
    
    co_clustering = mobster:::compute_co_clustering(x, resamples, bootstrap.fits)
    co_clustering_matrix = co_clustering$co.clustering
    co_clustering_labels = co_clustering$ordered.labels
    
    # cli::cli_process_done()
  }
  
  ret = list(
    bootstrap_values = bootstrap.values,
    bootstrap_model = model.frequency,
    bootstrap_statistics = stats,
    bootstrap_co_clustering = co_clustering_matrix,
    bootstrap_co_clustering_ordered_labels = co_clustering_labels
  )
  
  return(ret)
}
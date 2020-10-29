#' Plot the boostrapped model frequency
#' 
#' @description From the results of a call to \code{mobster_boostrap}, and
#' the results of a call to \code{bootstrapped_statistics}, a barplot with
#' the model frequency is produced
#'
#' @param bootstrap_results Results of a call to \code{mobster_boostrap}.
#' @param bootstrap_statisticsResults of a call to \code{bootstrapped_statistics}.
#'
#' @return A barplot of the boostrapped model frequency
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
#' plot_bootstrap_model_frequency(boot_results, boot_stats)
plot_bootstrap_model_frequency = function(bootstrap_results, bootstrap_statistics)
{
  is_bootstrap_results(bootstrap_results)
  is_bootstrap_statistics(bootstrap_statistics)
  
  # plot 
  n = length(bootstrap_results$fits)
  type = bootstrap_results$bootstrap
  
  ggplot(
    bootstrap_statistics$bootstrap_model,
    aes(x = Model, y = Frequency, fill = fit.model)
  ) +
    geom_bar(stat = 'identity') +
    ylim(0, 1) +
    scale_fill_manual(values = c(`TRUE` = 'forestgreen', `FALSE` = 'darkred')) +
    mobster:::my_ggplot_theme() + 
    labs(
      title = 'Bootstrap model frequency',
      subtitle = paste0('n = ', n, ' ', type, ' bootstraps.')
    ) +
    guides(fill = guide_legend('Bootstrapped model'))
}

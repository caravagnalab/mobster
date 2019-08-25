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

#' Plot the boostrapped mixing proportions
#' 
#' @description From the results of a call to \code{mobster_boostrap}, and
#' the results of a call to \code{bootstrapped_statistics}, a boxplot with
#' the mixing proportions is produced.
#'
#' @param x A MOBSTER fit.
#' @param bootstrap_results Results of a call to \code{mobster_boostrap}.
#' @param bootstrap_statisticsResults of a call to \code{bootstrapped_statistics}.
#' @param colors If provided, these colours will be used for each cluster.
#' If a subset of colours is provided, palette Set1 from \code{RColorBrewer} is used.
#' By default the tail colour is provided as 'gainsboro'.
#'
#' @return A barplot of the bootstrapped mixing proportions.
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
#' plot_bootstrap_mixing_proportions(x$best, boot_results, boot_stats)
plot_bootstrap_mixing_proportions = function(x,
                                             bootstrap_results, 
                                             bootstrap_statistics,
                                             colors = c(`Tail` = 'gainsboro')
                                             )
{
  is_mobster_fit(x)
  is_bootstrap_results(bootstrap_results)
  is_bootstrap_statistics(bootstrap_statistics)
  
  # plot 
  n = length(bootstrap_results$fits)
  type = bootstrap_results$bootstrap
  
  mixing = bootstrap_statistics$bootstrap_values %>%
    filter(statistics == 'Mixing proportion')
  
  fit = data.frame(fit = mobster:::.params_Pi(x)) %>% as_tibble()
  fit$cluster = names(mobster:::.params_Pi(x))
  
  bplot = ggplot(
    mixing,
    aes(x = cluster, y = fit.value, fill = cluster)
  ) +
    geom_violin(color = NA, alpha = .8) +
    geom_boxplot(alpha = 1, width = 0.05) +
    geom_jitter(alpha = 1, size = .5, height = 0, show.legend = F) +
    # geom_count() +
    # ggbeeswarm::geom_quasirandom(size = .5, show.legend = F) +
    geom_point(data = fit, inherit.aes = F, aes(x = cluster, y = fit), shape = 2, size =3) +
    ylim(0, 1) +
    # scale_fill_manual(values = c(`TRUE` = 'forestgreen', `FALSE` = 'darkred')) +
    mobster:::my_ggplot_theme() + 
    labs(
      title = 'Bootstrap mixing proportions',
      x = 'Cluster',
      y = 'Fit value',
      subtitle = paste0('n = ', n, ' ', type, ' bootstraps.')
    ) +
    guides(fill = guide_legend('Cluster')) +
    guides(color = guide_legend('Fit'))
  
  mobster:::add_fill_color_pl(x, bplot, colors)
}

#' Plot the boostrapped tail parameters.
#'
#' @description From the results of a call to \code{mobster_boostrap}, and
#' the results of a call to \code{bootstrapped_statistics}, a boxplot with
#' the tail parameters is produced.
#'
#' @param x A MOBSTER fit.
#' @param bootstrap_results Results of a call to \code{mobster_boostrap}.
#' @param bootstrap_statisticsResults of a call to \code{bootstrapped_statistics}.
#' @param colors If provided, these colours will be used for each cluster.
#' If a subset of colours is provided, palette Set1 from \code{RColorBrewer} is used.
#' By default the tail colour is provided as 'gainsboro'.
#'
#' @return A barplot of the bootstrapped tail parameters.
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
#' plot_bootstrap_tail(x$best, boot_results, boot_stats)
plot_bootstrap_tail = function(x,
                               bootstrap_results,
                               bootstrap_statistics,
                               colors = c(`Tail` = 'gainsboro'))
{
  is_mobster_fit(x)
  is_bootstrap_results(bootstrap_results)
  is_bootstrap_statistics(bootstrap_statistics)

  # plot
  n = length(bootstrap_results$fits)
  type = bootstrap_results$bootstrap

  tail = bootstrap_statistics$bootstrap_values %>%
    dplyr::filter(
      statistics == 'Shape' | statistics == 'Scale',
      cluster == 'Tail'
      )

  if(nrow(tail) == 0)
  {
    cli::cli_alert_warning("No tail have been fit in these bootstraps, returning an empty plot")
    return(CNAqc:::eplot())
  }

  scale = ggplot(data = tail,
                 aes(x = statistics,
                     y = fit.value,
                     fill = cluster)) +
    geom_violin(color = NA,
                alpha = .8,
                trim = F) +
    geom_boxplot(alpha = 1, width = 0.2) +
    # geom_count() +
    geom_jitter(
      alpha = 1,
      size = .5,
      height = 0,
      width = 0.05,
      show.legend = F
    ) +
    mobster:::my_ggplot_theme() +
    facet_wrap(~ statistics, scale = 'free') +
    labs(
      title = bquote("Bootstrap tail parameters"),
      y = "Fit value",
      x = "Parameter",
      subtitle = paste0('n = ', n, ' ', type, ' bootstraps.')
    ) +
    guides(fill = FALSE)

  scale = mobster:::add_fill_color_pl(x, scale, colors)
  scale
}

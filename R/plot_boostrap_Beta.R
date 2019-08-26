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
plot_bootstrap_Beta = function(x,
                               bootstrap_results,
                               bootstrap_statistics,
                               colors = c(`Tail` = 'gainsboro'))
{
  is_bootstrap_results(bootstrap_results)
  is_bootstrap_statistics(bootstrap_statistics)
  
  # plot
  n = length(bootstrap_results$fits)
  type = bootstrap_results$bootstrap
  
  plt = function(s)
  {
    betas = bootstrap_statistics$bootstrap_values %>%
      filter(statistics == s, 
             cluster != 'Tail')
      
    fit = x$Clusters %>% filter(type == s, 
                                cluster != 'Tail')
  
    rg = c(min(betas$fit.value), max(betas$fit.value)) * c(0.25, 1.75)
    
    bt = ggplot(data = betas,
                aes(fit.value,
                    fill = cluster)) +
      geom_histogram(bins = 100) +
      geom_vline(data = fit, aes(xintercept = fit.value), size = 0.3) +
      mobster:::my_ggplot_theme() +
      xlim(rg[1], rg[2]) +
      # facet( cluster ~ statistics, scale = 'free') +
      labs(
        title = bquote("Bootstrap Beta parameters"),
        y = "Fit value",
        x = paste("Beta", s),
        subtitle = paste0('n = ', n, ' ', type, ' bootstraps.')
      ) +
      guides(fill = FALSE)
    
    bt = mobster:::add_fill_color_pl(x, bt, colors)
    bt
  }
  
  mp = suppressWarnings(plt('Mean') + xlim(0, 1))
  vp = plt('Variance') + labs(title = bquote(" "), subtitle = ' ')
  ggpubr::ggarrange(mp, vp, ncol = 2)
}

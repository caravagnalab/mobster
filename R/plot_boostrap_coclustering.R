#' Plot the boostrapped co-clustering probability.
#'
#' @description From the results of a call to \code{mobster_boostrap}, and
#' the results of a call to \code{bootstrapped_statistics}, a figure (heatmap-style)
#' for the boostrapped co-clustering probability is produced assembling plots 
#' with \code{cowplot::plot_grid}. 
#'
#' @param x A MOBSTER fit.
#' @param bootstrap_results Results of a call to \code{mobster_boostrap}.
#' @param bootstrap_statisticsResults of a call to \code{bootstrapped_statistics}.
#' @param colors If provided, these colours will be used for each cluster.
#' If a subset of colours is provided, palette Set1 from \code{RColorBrewer} is used.
#' By default the tail colour is provided as 'gainsboro'.
#'
#' @return A \code{cowplot} figure of the bootstrapped co-clustering probability.
#'
#' @export
#'
#' @importFrom reshape2 melt
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
#' plot_bootstrap_coclustering(x$best, boot_results, boot_stats)
plot_bootstrap_coclustering = function(x,
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
  
  # Normalised matrix
  cocl = bootstrap_statistics$bootstrap_co_clustering/n
  ordered_labels = bootstrap_statistics$bootstrap_co_clustering_ordered_labels
  
  # cocl = bootstrap_statistics$bootstrap_co_clustering[ord_cl, ord_cl] /
  #   n
  # colnames(cocl) = rownames(cocl) = 1:ncol(cocl)
  # cocl[upper.tri(cocl)] = 0
  
  bt = ggplot(data = cocl %>% reshape2::melt(),
              aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    mobster:::my_ggplot_theme() +
    labs(
      title = bquote("Bootstrap co-clustering probability"),
      y = "Mutation",
      x = "Mutation",
      subtitle = paste0('n = ', n, ' ', type, ' bootstraps.')
    ) +
    scale_fill_gradientn(colours = c(alpha('steelblue', .3), 'forestgreen', 'darkorange', 'darkred'), 
                         values = c(0, 0.4, 0.8, 1),
                         aes(alpha = .7), limits = c(0, 1)) +
    # scale_fill_distiller(palette = 'Purples', direction = 1, aes(alpha = .7), limits = c(0, 1)) +
    # scale_fill_viridis_c(direction = -1, limits = c(0, 1), option = "plasma") +
    guides(fill = guide_colorbar(title = 'Probability', barwidth = unit(3, 'cm'))) 
  
  splits = cumsum(table(x$data$cluster))
  splits = data.frame(cluster = names(splits), 
                      n = as.vector(table(x$data$cluster)),
                      sum = splits, stringsAsFactors = FALSE)
  
  # add squares
  bt = bt + 
    geom_segment(data = splits %>% mutate(x = sum, y = 0, xend = sum, yend = sum),
                 inherit.aes = FALSE,
                 aes(x = x, y = y, xend = xend, yend = yend, color = cluster),
                 size = 1) +
    geom_segment(data = splits %>% mutate(x = sum - n, y = sum - n, xend = sum, yend = sum),
                 inherit.aes = FALSE,
                 aes(x = x, y = y, xend = xend, yend = yend, color = cluster),
                 size = 1) +
    geom_segment(data = splits %>% mutate(x = sum(n), y = sum, xend = sum, yend = sum),
                 inherit.aes = FALSE,
                 aes(x = x, y = y, xend = xend, yend = yend, color = cluster),
                 size = .5, linetype = 'dashed')+
    guides(color = FALSE) 
  
  bt = mobster:::add_color_pl(x, bt, colors)
  
  # Cluster assignments bar
  label = ggplot(data = x$data %>% mutate(id = row_number()),
                 aes(x = 1, y = id, fill = cluster)) +
    geom_tile() +
    mobster:::my_ggplot_theme() +
    labs(
      title = ' ',
      y = "MOBSTER Cluster",
      x = " ",
      subtitle = ' '
    ) +
    theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
      ) +
    guides(fill = guide_legend(''))+
    scale_y_continuous(position = "right")
  
  label = mobster:::add_fill_color_pl(x, label, colors)
  
  cowplot::plot_grid(bt, label, ncol = 2, rel_widths = c(0.9, 0.2), align = 'h')
}

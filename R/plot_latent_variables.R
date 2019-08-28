#' Plot the latent variables of the mixture.
#' 
#' @description It renders a heatmap where the latent variables
#' (reponsibilities) are shown and colured according to their value.
#' This function also calls function \code{Clusters}, using a parameter
#' that determines if a point is not to be assigned its best cluster
#' based on a cutoff.
#' 
#'
#' @param x A MOBSTER fit.
#' @param cutoff_assignment The parameter used to call function
#' \code{Clusters}, which does not assign a point to its best cluster
#' if the value of the corresponding latent variable is not above the cutoff.
#'
#' @return A plot of the latent variables of the mixture.
#' 
#' @export
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' plot_latent_variables(fit_example$best)
#' plot_latent_variables(fit_example$best, cutoff_assignment = .9)
plot_latent_variables = function(x, cutoff_assignment = 0)
{
  stopifnot(inherits(x, "dbpmm"))
  
  # assignments
  assignments = Clusters(x, cutoff_assignment) %>%
    dplyr::select(VAF, cluster) %>%
    data.frame(stringsAsFactors = FALSE)
  
  not_assign = is.na(assignments$cluster)
  n = sum(not_assign)
  p = (n/nrow(assignments)) * 100 
  ordering = order(assignments$cluster, na.last = TRUE)
  
  # Reshape and cut
  lv = reshape2::melt(x$z_nk[ordering, , drop = FALSE])
  
  # lv$value = cut(lv$value,
  #                breaks = c(-Inf, seq(0, 1, 0.05), Inf))
  # 
  colnames(lv) = c('Point', "Cluster", "Value")
  
  # 
  # colors = colorRampPalette(RColorBrewer::brewer.pal(8, 'YlGnBu'))(20)
  cuts_below = cuts_below = c()
  
  if(cutoff_assignment - 0.05 > 0) cuts_below  = seq(0, cutoff_assignment - 0.05, 0.05)
  if(cutoff_assignment - 0.05 < 1) cuts_above  = seq(cutoff_assignment, 1, 0.05)
  
  lv$Value = cut(lv$Value, breaks = c(-Inf, cuts_below, cuts_above, Inf))
  
  lblues = RColorBrewer::brewer.pal(5, 'Blues')
  lreds = RColorBrewer::brewer.pal(5, 'Reds')
  lreds = c('darkorange2', 'darkred')
  
  colors_below = colorRampPalette(lblues[1:3])(length(cuts_below))
  colors_above = colorRampPalette(lreds)(length(cuts_above))
  
  ggplot(lv, aes(x = Cluster, y = Point, fill = Value)) +
    geom_raster() +
    scale_fill_manual(values = c(colors_below, colors_above)) +
    my_ggplot_theme() +
    theme(
      legend.text = element_text(size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    guides(fill = guide_legend('')) +
    labs(
      title = bquote("Latent variables"),
      subtitle = bquote(
        .(n) ~" non assignable ("* .(p) *'%) with cutoff ' * z['nk']  ~' > ' * .(cutoff_assignment) ),
      y = paste0("Points (n =", x$N, ')')
    )
}
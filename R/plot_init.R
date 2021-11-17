#' Plot the initial density of a fit.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param colors If provided, these colours will be used for each cluster.
#' If a subset of colours is provided, palette Set1 from \code{RColorBrewer} is used.
#' By default the tail colour is provided as 'gainsboro'.
#'
#' @return A ggplot object for the plot.
#' @export
#'
#' @examples
#' data(fit_example)
#' plot_init(fit_example$best)
plot_init = function(x,
                     colors = c(`Tail` = 'gainsboro')
)
{
  is_mobster_fit(x)
  
  # Simple plot.
  #
  # We just get the desnity with the usual functions
  # after assinging the initial parameters appropriately.
  binwidth = 0.01
  domain = seq(0, 1, binwidth)
  
  n = x
  n$Clusters$fit.value = n$Clusters$init.value
  
  initial.densities = template_density(n,
                                       x.axis = domain[2:(length(domain) - 1)],
                                       # Restricted for numerical errors
                                       binwidth = 0.01,
                                       reduce = TRUE)
  
  den_init_pl = ggplot() +
    labs(title = bquote("Initialization"),
         x = "Observed Frequency",
         y = "Density") +
    guides(fill = FALSE) +
    ylim(0, max(initial.densities$y)) +
    my_ggplot_theme() +
    geom_line(data = initial.densities, aes(y = y, x = x, color = cluster)) +
    guides(color = FALSE)
  
  
  return(add_color_pl(x, den_init_pl, colors))
}
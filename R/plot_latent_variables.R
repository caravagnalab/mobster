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
#' @importFrom reshape2 melt
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' plot_latent_variables(fit_example$best)
#' plot_latent_variables(fit_example$best, cutoff_assignment = .9)
plot_latent_variables = function(x, cutoff_assignment = 0)
{
  mobster:::is_mobster_fit(x)

  # assignments with LV
  assignments = mobster::Clusters(x, cutoff_assignment) %>%
    dplyr::arrange(cluster)
  
  clusters_names = names(x$pi)
  
  # Statistics about non-assigned mutations
  not_assign = is.na(assignments$cluster)
  n = sum(not_assign)
  p = round((n/nrow(assignments)) * 100) 

  # Reshape and cut, preserving ordering on the y-axis
  lv = reshape2::melt(
    assignments %>% dplyr::select(clusters_names) %>% dplyr::mutate(pos = row_number()),
    id = 'pos'
  )
  
  colnames(lv) = c('Point', "Cluster", "Value")
  
  # TODO If there are drivers, we add the label
  # if(mobster:::has_drivers_annotated(x))
  # {
  #   drv_points = which(x$data$is_driver, arr.ind = TRUE)
  #   
  #   if(length(drv_points) > 0)
  #   {
  #     lv = lv %>%
  #       dplyr::mutate(
  #         label = ifelse(Point %in% drv_points, x$data$driver_label[Point], NA)
  #       )
  #   }
  # }
  

  ggplot(lv, aes(x = Cluster, y = Point, fill = Value)) +
    geom_raster() +
    scale_fill_viridis_c(direction = -1) +
    mobster:::my_ggplot_theme() +
    guides(fill = guide_colorbar(bquote(z['nk'] ~ ' '), barwidth = unit(3, 'cm'))) +
    labs(
      title = bquote("Latent variables"),
      subtitle = bquote(
        z['nk']  ~' > ' * .(cutoff_assignment) ~ ': n =' ~.(n) ~" NAs ("* .(p) *'%)'),
      y = paste0("Points (n =", x$N, ')')
    )
}

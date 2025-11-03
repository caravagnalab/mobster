#' Scatterplots of a MOBSTERm fit.
#'
#' @param x A MOBSTERm fit.
#'
#' @return A list of ggplot objects for the scatter plots.
#' @export
#'
#' @examples TBD
plot.dbpmm_m = function(x){

  cm = combn(x$sample_names, 2) # Generate all combinations of the elements of sample_names taken 2 at a time
  
  plots <- apply(
    cm,
    2,
    function(w) mobster:::plotm_2D(x, d1 = w[1], d2 = w[2]) # w is the sample name
  )
  
  plots

}

my_ggplot_theme = function(cex = 1)
{
  theme_light(base_size = 10 * cex) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      panel.background = element_rect(fill = 'white')
    )
}

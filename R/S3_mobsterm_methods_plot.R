#' Scatterplots of a MOBSTERm fit.
#'
#' @param x A MOBSTERm fit.
#'
#' @return A list of ggplot objects for the scatter plots.
#' @export
#'
#' @examples TBD
plot.dbpmm_m = function(x, color_palette=NA)
{
  if (is.null(color_palette) || all(is.na(unlist(color_palette)))) {
    color_palette = c(
      "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#a65628",
      "#FFD700",  "#999999", "#000000", "#f781bf", # First 10 colors (Set1)  
      "#46f0f0", "#f032e6", "#bcf60c", "#fabed4", "#008080", "#e6beff",  
      "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",  
      "#000075", "#808080", "#d3a6f3", "#ff9cdd", "#73d7b0"  
    ) 
  }
  cm = combn(x$sample_names, 2) # Generate all combinations of the elements of sample_names taken 2 at a time
  
  plots <- apply(
    cm,
    2,
    function(w) mobster:::plotm_2D(x, d1 = w[1], d2 = w[2], color_palette) # w is the sample name
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

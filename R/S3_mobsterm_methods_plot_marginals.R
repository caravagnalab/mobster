#' Marginal plots of a MOBSTERm fit.
#'
#' @param x A MOBSTERm fit.
#'
#' @return A list of ggplot objects for the marginal plots.
#' @export
#'
#' @examples TBD
plot_marginals.dbpmm_m = function(x, color_palette=NA)
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
  plots = mobster:::plotm_1D(x, color_palette)
  
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

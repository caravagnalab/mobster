#' Marginal plots of a MOBSTERm fit.
#'
#' @param x A MOBSTERm fit.
#'
#' @return A list of ggplot objects for the marginal plots.
#' @export
#'
#' @examples TBD
plot_marginals.dbpmm_m = function(x){
  
  plots <- lapply(x$sample_names,
    function(s) mobster:::plotm_1D(x, s) # s is the sample name
    # function(s) plotm_1D(x, s) # s is the sample name
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

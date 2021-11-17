#' Plot summary for model selection.
#'
#' @description This plot is usefull to understand the different
#' fits and their rank with respect to some scoring used to
#' select the best model. The plot shows alternative solutions
#' and their rank as well.
#'
#' @param x A list of fits computed via \code{mobster_fit}.
#' @param TOP The first \code{TOP} fits are used (by ranking).
#'
#' @return A complex figure with all plots arranged using both
#' \code{ggpubr} and \code{cowplot}.
#' 
#' @importFrom cowplot plot_grid
#'
#' @export
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' plot_model_selection(fit_example)
plot_model_selection = function(x,
                                TOP = 6,
                                nx = 3,
                                ny = 2,
                                ...)
{
  # =-=-=-=-=-=-=-=-
  # Top and best
  # =-=-=-=-=-=-=-=-
  models = plot_grid(x, TOP, nx, ny)
  
  # =-=-=-=-=-=-=-=-
  # Scores rank best
  # =-=-=-=-=-=-=-=-
  rank = plot_fit_scores(x)
  
  # =-=-=-=-=-=-=-=-
  # Goodness of fit
  # =-=-=-=-=-=-=-=-
  GOFP = suppressWarnings(plot_gofit(x, TOP = TOP, nx, ny))
  
  # =-=-=-=-=-=-=-=-
  # Assembly
  # =-=-=-=-=-=-=-=-
  bottom = cowplot::plot_grid(
    rank,
    GOFP,
    nrow = 1,
    ncol = 2,
    align = 'h',
    axis = "bt"
  )
  
  full_figure = cowplot::plot_grid(models,
                                   bottom,
                                   nrow = 2,
                                   ncol = 1,
                                   align = 'v')
  full_figure
}

plot_grid = function(x,
                     TOP = 6,
                     nx = 3,
                     ny = 2)
{
  is_list_mobster_fits(x)
  stopifnot((nx * ny) >= TOP)
  
  N = length(x$runs)
  TOP = min(TOP, N)
  
  plots_g = lapply(1:TOP, function(y) {
    plot(x$runs[[y]]) +
      labs(title = NULL,
           subtitle = NULL,
           caption = paste0(' ', y, ' / ', N)) +
      guides(fill = FALSE) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
      )
  })
  
  figure = ggpubr::ggarrange(plotlist = plots_g,
                             nrow = ny,
                             ncol = nx)
  
  ggpubr::ggarrange(plot(x$best),
                    figure,
                    # labels = c("1) ", ""),
                    nrow = 1,
                    ncol = 2)
  
}
# 
# plot_model_selection = function(x,
#                                 TOP = 5,
#                                 ...)
# {
#   is_list_mobster_fits(x)
#   
#   # =-=-=-=-=-=-=-=-
#   # Best model fit
#   # =-=-=-=-=-=-=-=-
#   best = plot.dbpmm(x$best)
#   
#   # =-=-=-=-=-=-=-=-
#   # Scores rank best
#   # =-=-=-=-=-=-=-=-
#   rank = plot_fit_scores(x)
#   
#   # =-=-=-=-=-=-=-=-
#   # Goodness of fit
#   # =-=-=-=-=-=-=-=-
#   
#   # SSE plot
#   GOFP = suppressWarnings(plot_gofit(x, TOP = TOP))
#   
#   # =-=-=-=-=-=-=-=-=-=-
#   # Low-rank solutions
#   # =-=-=-=-=-=-=-=-=-=-
#   
#   # Other solutions ranked below top best -- maximum TOP 4 fixed
#   TOP.plot = min(4, nrow(x$fits.table))
#   
#   other.best = NULL
#   
#   if (TOP > 1)
#     other.best = lapply(2:TOP,
#                         function(w)
#                           plot.dbpmm(x$runs[[w]]) +
#                           labs(title = paste("Solution rank:", w)) +
#                           my_ggplot_theme(cex = .6))
#   
#   
#   # =-=-=-=-=-=-=-=-
#   # Final layout
#   # =-=-=-=-=-=-=-=-
#   bottom = ggarrange(rank, GOFP, ncol = 2, nrow = 1)
#   
#   left_panel = ggarrange(best,
#                          bottom,
#                          nrow = 2,
#                          ncol = 1)
#   
#   
#   k = ceiling(sqrt(length(other.best)))
#   
#   right_panel = ggarrange(plotlist = other.best,
#                           ncol = k,
#                           nrow = k)
#   
#   figure = ggarrange(
#     left_panel,
#     right_panel,
#     widths = c(1.5, 1),
#     nrow = 1,
#     ncol = 2
#   )
#   
#   figure = annotate_figure(figure,
#                            top = text_grob(
#                              "MOBSTER model selection",
#                              color = "black",
#                              face = "bold",
#                              size = 18
#                            ))
#   
#   return(figure)
# }
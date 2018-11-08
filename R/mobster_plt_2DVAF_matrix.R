
#' Title
#'
#' @param x 
#' @param title 
#' @param ... 
#' @param samples
#' @param lower.cluster 
#' @param lower.cluster.label 
#' @param upper.cluster 
#' @param upper.cluster.label 
#' @param palette 
#' @param MOBSTER 
#'
#' @return
#' @import pio
#'
#' @export
#' 
#' @examples
mobster_plt_2DVAF_matrix = function(
  lower.x,
  lower.cluster = NULL,
  lower.cluster.label = NULL,  
  lower.exclude = NULL,
  upper.x = NULL,
  upper.cluster = NULL,
  upper.cluster.label = NULL,  
  upper.exclude = NULL,
  MOBSTER = NULL,
  samples = lower.x$samples,
  lower.description = lower.x$description,
  upper.description = ifelse(is.null(upper.x), "" , upper.x$description),
  palette = list(MOBSTER = 'Set1', lower = "Set2", upper = "Spectral"), 
  bg = list(MOBSTER = 'gainsboro', lower = "darkslategray4", upper = "ivory3"),
  ...)
{
  # The numbe of plots to compute, and a template empty plot
  k = length(samples)
  layout = matrix(1:(k*k), ncol = k, byrow = T)
  
  eplot = ggplot() + geom_blank() + theme_void()
  plots = NULL
  
  # check input for MOBSTER diagonoal plots
  if(!is.null(MOBSTER)) {
    stopifnot(is.list(MOBSTER))
    stopifnot(length(MOBSTER) == k)
  }
  
  # Matrix layout made explicit
  for (s in 1:k) {
    for (w in 1:k) {
      fig = eplot
      
      # Diagonal is MOBSTER if it exists
      if(s == w & !is.null(MOBSTER)) 
        fig = plot(
          MOBSTER[[s]]$best,
          palette = palette$MOBSTER,
          histogram.main = paste("MOBSTER ", samples[w]),
          silent = TRUE, ...)$mainHist  +
            theme(plot.background = element_rect(fill = bg$MOBSTER))
      
      # Lower triangular is a 2D VAF plot with lower.cluster
      if(s > w & !is.null(lower.cluster)) 
        fig = mobster_plt_2DVAF(
          obj = lower.x, 
          x = samples[s], 
          y = samples[w],
          cluster = lower.cluster, 
          cluster.label = lower.cluster.label,
          palette = palette$lower,
          exclude = lower.exclude,
          ...) +
            theme(plot.background = element_rect(fill = bg$lower))
      
      # Upper triangular is a 2D VAF plot with upper.cluster
      if(s < w & !is.null(upper.cluster)) 
        fig = mobster_plt_2DVAF(
          obj = upper.x, 
          x = samples[s], 
          y = samples[w],
          cluster = upper.cluster, 
          cluster.label = upper.cluster.label,
          palette = palette$upper,
          exclude = upper.exclude,
          ...) +
            theme(plot.background = element_rect(fill = bg$upper))
      
      # queue the figure
      plots = append(plots, list(fig))
    }
  }
  
  # Now arrange the plots in a k x k matrux
  figure = ggarrange(
    plotlist = plots,
    ncol = k,
    nrow = k
  )
  
  figure = annotate_figure(figure, 
                           top =  text_grob("2D VAF plot", color = "black", hjust = 1, x = 1, face = "bold", size = 25), 
                           left = text_grob(
                             paste0("Lower triangular: ", lower.description, '\n'), 
                             color = bg$lower, face = "bold", size = 20, rot = 90),
                           right = text_grob(
                             paste0("Upper triangular: ", upper.description, '\n'), 
                             color = bg$upper, face = "bold", size = 20, rot = 90),
                           bottom = text_grob(
                             paste0("\nDiagonal: ", ifelse(!is.null(MOBSTER), "MOBSTER fits", "")), 
                             color = bg$MOBSTER, face = "bold", size = 20))
  
  figure
  
}



# 
# MB.figure = best.MOBSTER.plots = NULL
# if(MOBSTER)
# {
#   best.MOBSTER = lapply(x$fit.MOBSTER, function(w) w$best)
#   best.MOBSTER.plots = mobster:::plot_diagonal_MOBSTER(best.MOBSTER, samples, ...)
#   
#   MB.figure = ggpubr::ggarrange(
#     plotlist = best.MOBSTER.plots,
#     nrow = 1,
#     ncol = length(best.MOBSTER),
#     labels = LETTERS[seq(best.MOBSTER)]
#   )
#   
#   panel.labels = panel.labels[-seq(best.MOBSTER)]
# }
# 
# plots = NULL
# id = 1
# 
# for (s in seq(samples)) {
#   for (w in s:length(samples)) {
#     if (s != w) {
#       
#       pl.1 = mobster_plt_2DVAF(
#         obj = x, 
#         x = samples[s], 
#         y = samples[w],
#         cluster = cluster, 
#         cluster.label = cluster.label,
#         ...)
# 
#       fig = ggpubr::ggarrange(pl.1, nrow = 1, ncol = 1, labels = panel.labels[id])
#       
#       plots = append(plots, list(fig))
#       id = id + 1
#     }
#   }
# }
# 
# k = ceiling(sqrt(length(plots)))
# k = length(plots)/3
# 
# if(k < 1) k = 1
# 
# twoBtwo = ggpubr::ggarrange(
#   plotlist = plots,
#   ncol = 3,
#   nrow = k
# )
# 
# figure = ggpubr::ggarrange(
#   MB.figure, 
#   twoBtwo,
#   ncol = 1,
#   nrow = 2,
#   heights = c(.25, 1)
# )


# S = length(samples)
# 
# layout = matrix(0, ncol = S, nrow = S)
# if(MOBSTER) diag(layout) = 1:S
# 
# combs = S * (S-1) / 2
# layout[lower.tri(layout)] = (1:combs) + max(layout)
# 
# mobster:::.multiplot(plotlist = append(best.MOBSTER.plots, plots), layout = layout, title = title)
# # plotlist = append(list(MB.figure), plots)

# figure = ggpubr::ggarrange(
#   plotlist = plotlist,
#   nrow = length(plotlist),
#   ncol = 1,
#   common.legend = T,
#   heights = c(.75, rep(1, length(plotlist)))
# )
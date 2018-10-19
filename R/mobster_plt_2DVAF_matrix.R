
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
  upper.x = NULL,
  samples = lower.x$samples,
  MOBSTER = !is.null(lower.x$fit.MOBSTER),
  lower.cluster = NULL,
  lower.cluster.label = 'cluster',  
  upper.cluster = NULL,
  upper.cluster.label = 'cluster',  
  title = lower.x$description,
  palette = list(MOBSTER = 'Set1', lower = "Set2", upper = "Spectral"), 
  ...)
{
  # pio::pioHdr("MOBSTER - Multivariate fits plot",
  #             c(`Cluster 1` = cluster1,
  #               `Cluster 2` = cluster2,
  #               `Samples` = paste(samples, collapse = ',')
  #             )
  # )
  # 
  # 
  require(ggplot2)
  
  k = length(samples)
  layout = matrix(1:(k*k), ncol = k, byrow = T)
  
  eplot = ggplot() + geom_blank() + theme_void()
  plots = NULL
  
  for (s in 1:k) {
    for (w in 1:k) {
      fig = eplot
      
      # Diagonal is MOBSTER if it exists
      if(s == w & MOBSTER) fig = plot(
        lower.x$fit.MOBSTER[[s]]$best,
        palette = palette$MOBSTER,
        histogram.main = paste("MOBSTER ", samples[w]),
        silent = TRUE, ...)$mainHist
      
      # Lower triangular is a 2D VAF plot with lower.cluster
      if(s > w & !is.null(lower.cluster)) fig = mobster_plt_2DVAF(
        obj = lower.x, 
        x = samples[s], 
        y = samples[w],
        cluster = lower.cluster, 
        cluster.label = lower.cluster.label,
        palette = palette$lower,
        ...)
      
      # Upper triangular is a 2D VAF plot with upper.cluster
      if(s < w & !is.null(upper.cluster)) fig = mobster_plt_2DVAF(
        obj = upper.x, 
        x = samples[s], 
        y = samples[w],
        cluster = upper.cluster, 
        cluster.label = upper.cluster.label,
        palette = palette$upper,
        ...)
      
      # queue the figure
      plots = append(plots, list(fig))
    }
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
  
  figure = ggpubr::ggarrange(
    plotlist = plots,
    ncol = k,
    nrow = k
  )
  
  figure = ggpubr::annotate_figure(figure, 
                                   top = title, 
                                   left = ifelse(is.null(lower.cluster), "" , lower.cluster.label),
                                   right = ifelse(is.null(upper.cluster), "" , upper.cluster.label))
  
  figure
  
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
}


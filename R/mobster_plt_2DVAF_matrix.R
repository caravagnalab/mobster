
#' Title
#'
#' @param x 
#' @param cluster 
#' @param cluster.label 
#' @param title 
#' @param ... 
#' @param samples
#' @param MOBSTER 
#' @param cex 
#' @param panel.labels 
#'
#' @return
#' @import pio
#'
#' @export
#' 
#' @examples
mobster_plt_2DVAF_matrix = function(
  x,
  samples = x$samples,
  MOBSTER = !is.null(x$fit.MOBSTER),
  cluster = NULL,
  cluster.label = 'cluster',  
  title = x$description,
  cex = 1,
  panel.labels = LETTERS,
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
  
  
  MB.figure = best.MOBSTER.plots = NULL
  if(MOBSTER)
  {
    best.MOBSTER = lapply(x$fit.MOBSTER, function(w) w$best)
    best.MOBSTER.plots = mobster:::plot_diagonal_MOBSTER(best.MOBSTER, samples, cex = cex)
    
    MB.figure = ggpubr::ggarrange(
      plotlist = best.MOBSTER.plots,
      nrow = 1,
      ncol = length(best.MOBSTER),
      labels = LETTERS[seq(best.MOBSTER)]
    )
    
    panel.labels = panel.labels[-seq(best.MOBSTER)]
  }

  plots = NULL
  id = 1
  
  for (s in seq(samples)) {
    for (w in s:length(samples)) {
      if (s != w) {
        
        pl.1 = mobster_plt_2DVAF(
          obj = x, 
          x = samples[s], 
          y = samples[w],
          cluster = cluster, 
          cluster.label = cluster.label)

        fig = ggpubr::ggarrange(pl.1, nrow = 1, ncol = 1, labels = panel.labels[id])
        
        plots = append(plots, list(fig))
        id = id + 1
      }
    }
  }
  
  k = ceiling(sqrt(length(plots)))

  twoBtwo = ggpubr::ggarrange(
    plotlist = plots,
    ncol = k,
    nrow = k
  )
  
  figure = ggpubr::ggarrange(
    MB.figure, 
    twoBtwo,
    ncol = 1,
    nrow = 2,
    heights = c(.15, 1)
  )
  
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


#' #' Title
#' #'
#' #' @param data
#' #' @param samples
#' #' @param VAF.range
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' mobster_plt_2DVAF_matrix = function(
#'   x,
#'   samples = x$samples,
#'   cluster1 = 'sciClone.cluster',
#'   cluster2 = 'MOBSTER.sciClone.cluster',
#'   palette.cluster1 = 'Set1',
#'   palette.cluster2 = 'Spectral',
#'   cex = 1,
#'   VAF.range = c(0.05, 0.6))
#' {
#'   
#'   pio::pioHdr("MOBSTER - Multivariate fits plot",
#'               c(`Cluster 1` = cluster1,
#'                 `Cluster 2` = cluster2,
#'                 `Samples` = paste(samples, collapse = ',')
#'               )
#'   )
#'   
#'   
#'   require(ggplot2)
#'   
#'   best.MOBSTER = lapply(x$fit.MOBSTER, function(w) w$best)
#'   
#'   #  Diagonal
#'   MB.figure = ggpubr::ggarrange(
#'     plotlist = mobster:::plot_diagonal_MOBSTER(best.MOBSTER, samples, cex = cex),
#'     nrow = 1,
#'     ncol = length(best.MOBSTER),
#'     labels = LETTERS[seq(best.MOBSTER)]
#'   )
#'   
#'   labels = LETTERS[-seq(best.MOBSTER)]
#'   
#'   plots = NULL
#'   id = 1
#'   
#'   for (s in seq(samples)) {
#'     for (w in s:length(samples)) {
#'       if (s != w) {
#'         
#'         pl.1 = mobster_plt_2DVAF(
#'           fit, 
#'           x = samples[s], 
#'           y = samples[w],
#'           cluster = SClusters(fit), 
#'           cluster.label = )
#'         
#'         pl.1 = plot_2DVAF(
#'           data$data,
#'           x = paste0('VAF.', samples[s]),
#'           y = paste0('VAF.', samples[w]),
#'           cluster = cluster1,
#'           palette = palette.cluster1,
#'           VAF.range = VAF.range,
#'           cex = cex
#'         )
#'         
#'         pl.2 = plot_2DVAF(
#'           data$data,
#'           x = paste0('VAF.projected.', samples[s]),
#'           y = paste0('VAF.projected.', samples[w]),
#'           cluster = cluster2,
#'           palette = palette.cluster2,
#'           VAF.range = VAF.range,
#'           cex = cex
#'         )
#'         
#'         fig = ggpubr::ggarrange(pl.1,
#'                                 pl.2,
#'                                 nrow = 1,
#'                                 ncol = 2,
#'                                 labels = labels[id])
#'         
#'         plots = append(plots, list(fig))
#'         id = id + 1
#'       }
#'     }
#'   }
#'   
#'   twoBtwo = ggpubr::ggarrange(
#'     plotlist = plots,
#'     ncol = 2,
#'     nrow = length(plots)/2
#'   )
#'   
#'   figure = ggpubr::ggarrange(
#'     MB.figure,
#'     twoBtwo,
#'     ncol = 1,
#'     nrow = 2,
#'     heights = c(.15, 1)
#'   )
#'   
#'   figure
#'   
#'   # plotlist = append(list(MB.figure), plots)
#'   # figure = ggpubr::ggarrange(
#'   #   plotlist = plotlist,
#'   #   nrow = length(plotlist),
#'   #   ncol = 1,
#'   #   common.legend = T,
#'   #   heights = c(.75, rep(1, length(plotlist)))
#'   # )
#' }

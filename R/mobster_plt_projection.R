#' Title
#'
#' @param x 
#' @param type 
#'
#' @return
#' @export
#'
#' @examples
mobster_plt_projection = function(x, type = 'global')
{
  if (type == "global")
  {
    gclusters = Clusters(x, 'M')
    
    gclusters$projected = NA
    for(i in 1:nrow(gclusters)) 
    {
      # samples where it is tail
      where = colnames(gclusters)[which(gclusters[i, ] == 'Tail')]
      # if(length(where) > 0) {
      #   labels = sort(gsub('cluster.', '', where))
      #   gclusters$projected[i] = paste0(labels, collapse = ', ')
      # }
      gclusters$projected[i] = paste0(length(where), '/', length(x$samples))
    }
    
    figure = mobster_plt_2DVAF_matrix(
      lower.x = x,
      lower.cluster = gclusters,
      lower.cluster.label = 'projected',
      palette = list(lower = 'Spectral', MOBSTER = 'Set1')
    )
    figure = ggpubr::annotate_figure(figure, bottom = paste0('Projection plot: ', type))
  }
  else {
    stop('TODO')
    
    # Local ones require a different approach
    # we have to explicit the structure of the clusters
    require(ggplot2)
    
    k = length(x$samples)
    layout = matrix(1:(k*k), ncol = k, byrow = T)
    
    eplot = ggplot() + geom_blank() + theme_void()
    plots = NULL
    
    for (s in 1:k) {
      for (w in 1:k) {
        fig = eplot
        
        # Diagonal is MOBSTER if it exists
        if(s == w & !is.null(x$fit.MOBSTER)) 
          fig = plot(
            x$fit.MOBSTER[[s]]$best,
            palette = palette$MOBSTER,
            histogram.main = paste("MOBSTER ", samples[w]),
            silent = TRUE)$mainHist
        
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
  }
}


# ggsave(
#   plot = mobster_plt_2DVAF_matrix(lower.x = x,
#                                   # upper.x = fitlp,
#                                   lower.cluster = clusters, lower.cluster.label = 'projected',
#                                   palette = list(lower = 'Spectral', MOBSTER = 'Set1')),
#   filename = 'a.pdf',
#   width = 15,
#   height = 15
# )
# 


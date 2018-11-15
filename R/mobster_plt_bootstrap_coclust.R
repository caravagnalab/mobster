
#' Co-clustering probability from non-parametric bootstrap 
#' 
#' @description Compute co-clustering probability from non-parametric bootstrap resamples
#' 
#' This function can be applied only to nonparametric bootstrap resamples.
#'
#' @param resamples Data resampled
#' @param fits Fits (one per resample)
#'
#' @return Nothing
#' @export
#'
#' @examples
mobster_plt_bootstrap_cocluster = function(x, resamples, fits, background.color = 'white') 
{
  # sum up occurrences 
  .coocc = function(l, M) 
  {
    cluster.labels = unique(l)
    
    for (cl in cluster.labels) 
    {
      # A unique is for nonparametric bootstrap
      # where we resample the same samples twice
      cl.assignments = unique(names(l[l == cl]))
      if(is.null(cl.assignments) | length(cl.assignments) == 1) next;
      
      pairs = combn(cl.assignments, 2, simplify = F)
      
      for (p in 1:length(pairs)) 
      {
        M[pairs[[p]][1], pairs[[p]][2]] = M[pairs[[p]][1], 
                                            pairs[[p]][2]] + 1
        M[pairs[[p]][2], pairs[[p]][1]] = M[pairs[[p]][2], 
                                            pairs[[p]][1]] + 1
      }
    }
    M
  }
  
  ########################## Analyze outputs
  # -- Co-clustering probability
  #
  # Note: only with nonparametric bootstrap
  
  # number of resamples
  n = length(resamples)
  
  # number of mutations
  N = nrow(fits[[1]]$data)
  
  if(!('original.id' %in% colnames(fits[[1]]$data)))
    stop('Missing the original.id column in the bootstrap data! \n\nAre you sure this is
         a result from non parametric bootstrap?')
  
  # Still assuming that id is the row index of the VAF
  co.clustering = matrix(0, nrow = N, ncol = N)
  
  rownames(co.clustering) =
    colnames(co.clustering) = 1:N
  
  # Extract co-clustering labels
  pb = txtProgressBar(0, length(fits), style = 3)
  
  for(w in seq(fits))
  {    
    setTxtProgressBar(pb, w)
    
    cluster.results = fits[[w]]$data
    
    cluster.labels = cluster.results$cluster
    names(cluster.labels) = cluster.results$original.id
    
    co.clustering = .coocc(cluster.labels, co.clustering)
  }
  
  ordered.data = x$data[, c('VAF', 'cluster')]
  rownames(ordered.data) = 1:N
  
  # sort heatmap by cluster
  ordering = order(ordered.data$cluster)
  
  ordered.data = ordered.data[ordering, ]
  co.clustering = co.clustering[ordering, ordering]
  
  # normalize values
  co.clustering = co.clustering/n
  
  colors = mobster:::scols(seq(0.1, 1, 0.05), palette = 'YlGnBu')
  colors = c(`0` = background.color, colors)
  
  ann.colors = mobster:::getColors_model(x)
  ann.VAF =  mobster:::scols(seq(0, 1, 0.05), palette = 'Oranges')
  
  # pheatmap::pheatmap(co.clustering[1:500, 1:500],
  figure = pheatmap::pheatmap(co.clustering,
                     main = 'Co-clustering (non parametric bootstrap)',
                     cluster_rows = FALSE,
                     cluster_cols = FALSE, 
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     border = NA,
                     color = colors,
                     annotation_row = ordered.data,
                     annotation_col = ordered.data, 
                     silent = TRUE,
                     annotation_colors = list(cluster = ann.colors, VAF = ann.VAF)
  )
  
  return(list(figure = figure, co.clustering = co.clustering))
}


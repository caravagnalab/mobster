# Compute parametric bootstrap replicate
.parametric_bootstrap_resamples = function(x, n = 100)
{
  stopifnot(n >= 1)
  
  data.size = nrow(x$data)
  
  lapply(1:n,
         function(w) {
           resample = rdbpmm(x, n = data.size)
           data.frame(id = 1:data.size, VAF = resample)
         })
}

# Compute nonparametric bootstrap replicate
.nonparametric_bootstrap_resamples = function(x, n = 100)
{
  stopifnot(n >= 1)
  
  data.size = nrow(x$data)
  
  lapply(1:n,
         function(w) {
           ids = sample(1:data.size, data.size, replace = TRUE)
           data.frame(
             id = 1:data.size,
             VAF = x$data$VAF[ids],
             original.id = ids
           )
         })
}

# Compute bootstrap co-clustering probability from nonparametric bootstraps
compute_co_clustering = function(x, resamples, fits)
{
  # sum up occurrences 
  .coocc = function(l, M) 
  {
    cluster.labels = unique(l)
    
    for (cl in cluster.labels) 
    {
      # A unique is for nonparametric bootstrap
      # where we resample the same point twice
      cl.assignments = unique(names(l[l == cl]))
      if(is.null(cl.assignments) | length(cl.assignments) == 1) next;
      
      # Sort them so the matrix will be lower diagonal
      cl.assignments = sort(cl.assignments %>% as.numeric)
      
      # Use an index to optimise (July 2020)
      for (p in 1:(length(cl.assignments) - 1))  
      {
        # The p-th is what we look at, we take that element 
        # and all those coming after it (p+1), (p+2), ...
        pointer = cl.assignments[p]
        targets = cl.assignments[(p+1):length(cl.assignments)]
        
        # We increase them by 1
        M[pointer, targets] = M[pointer, targets] + 1
      }
      
      # pairs = combn(cl.assignments, 2, simplify = F)
      # 
      # for (p in 1:length(pairs))  
      # {
      #   M[pairs[[p]][1], pairs[[p]][2]] = M[pairs[[p]][1], 
      #                                       pairs[[p]][2]] + 1
      #   M[pairs[[p]][2], pairs[[p]][1]] = M[pairs[[p]][2], 
      #                                       pairs[[p]][1]] + 1
      # }
    }
    M
  }
  
  ########################## Analyze outputs
  # -- Co-clustering probability
  #
  # Note: only with nonparametric bootstrap
  
  # number of resamples
  n = length(fits)
  
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
  # sp1 <- make_spinner()
  pr_bar = dplyr::progress_estimated(length(fits))
  
  for(w in seq(fits))
  {    
    # sp1$spin();
    pr_bar$tick()$print()
    
    cluster.results = mobster::Clusters(fits[[w]], cutoff_assignment = 0)
    
    cluster.labels = cluster.results$cluster
    names(cluster.labels) = cluster.results$original.id
    
    co.clustering = .coocc(l = cluster.labels, M = co.clustering)
  }
  
  # sp1$finish()
  
  ordered.data = x$data[, c('VAF', 'cluster')]
  ordered.data$input_point = 1:N
  
  # sort heatmap by cluster
  ordering = order(ordered.data$cluster)
  
  ordered.data = ordered.data[ordering, ]
  co.clustering = co.clustering[ordering, ordering]
  
  co.clustering[co.clustering > n] = n
  return(list(co.clustering = co.clustering, ordered.labels = ordered.data))
}

# compute_co_clustering = function(x, resamples, fits)
# {
#   if(!('original.id' %in% colnames(fits[[1]]$data)))
#     stop('Missing the original.id column in the bootstrap data! \n\nAre you sure this is
#          a result from non parametric bootstrap?')
#   
#   N = nrow(fits[[1]]$data)
#   
#   # Still assuming that id is the row index of the VAF
#   co.clustering = matrix(0, nrow = N, ncol = N)
#   
#   rownames(co.clustering) =
#     colnames(co.clustering) = 1:N
#   
#   f = fits[[1]]
#   
#   clusters  = Clusters(f, cutoff_assignment = 0)
#   clusters_label = unique(clusters$cluster)
#   
#   for(cl in clusters_label){
#     ids = clusters %>% filter(cluster == cl) %>% pull(original.id) 
#     ids = unique(ids)
#     
#     co_assign = 43
#     
#   }
#   
#   # sum up occurrences 
#   .coocc = function(l, M) 
#   {
#     cluster.labels = unique(l)
#     
#     for (cl in cluster.labels) 
#     {
#       # A unique is for nonparametric bootstrap
#       # where we resample the same samples twice
#       cl.assignments = unique(names(l[l == cl]))
#       if(is.null(cl.assignments) | length(cl.assignments) == 1) next;
#       
#       pairs = combn(cl.assignments, 2, simplify = F)
#       
#       for (p in 1:length(pairs)) 
#       {
#         M[pairs[[p]][1], pairs[[p]][2]] = M[pairs[[p]][1], 
#                                             pairs[[p]][2]] + 1
#         M[pairs[[p]][2], pairs[[p]][1]] = M[pairs[[p]][2], 
#                                             pairs[[p]][1]] + 1
#       }
#     }
#     M
#   }
#   
#   ########################## Analyze outputs
#   # -- Co-clustering probability
#   #
#   # Note: only with nonparametric bootstrap
#   
#   # number of resamples
#   n = length(resamples)
#   
#   # number of mutations
#   N = nrow(fits[[1]]$data)
#   
#   if(!('original.id' %in% colnames(fits[[1]]$data)))
#     stop('Missing the original.id column in the bootstrap data! \n\nAre you sure this is
#          a result from non parametric bootstrap?')
#   
#   # Still assuming that id is the row index of the VAF
#   co.clustering = matrix(0, nrow = N, ncol = N)
#   
#   rownames(co.clustering) =
#     colnames(co.clustering) = 1:N
#   
#   # Extract co-clustering labels
#   pb = txtProgressBar(0, length(fits), style = 3)
#   
#   for(w in seq(fits))
#   {    
#     setTxtProgressBar(pb, w)
#     
#     cluster.results = fits[[w]]$data
#     
#     cluster.labels = cluster.results$cluster
#     names(cluster.labels) = cluster.results$original.id
#     
#     co.clustering = .coocc(cluster.labels, co.clustering)
#   }
#   
#   ordered.data = x$data[, c('VAF', 'cluster')]
#   ordered.data$id = 1:N
#   
#   # sort heatmap by cluster
#   ordering = order(ordered.data$cluster)
#   
#   ordered.data = ordered.data[ordering, ]
#   co.clustering = co.clustering[ordering, ordering]
#   
#   co.clustering
# }


is_bootstrap_results = function(x)
{
  stopifnot('bootstrap' %in% names(x))
  stopifnot('resamples' %in% names(x))
  stopifnot('fits' %in% names(x))
}

is_bootstrap_statistics = function(x)
{
  stopifnot('bootstrap_model' %in% names(x))
  stopifnot('bootstrap_statistics' %in% names(x))
  stopifnot('bootstrap_co_clustering' %in% names(x))
}

get_mixing_proportions = function(data, cluster)
{
  pi = as.vector(table(data[, cluster])) /
    sum(!is.na(data[, cluster]))

  pi.df = data.frame(sort(pi, decreasing = TRUE))
  colnames(pi.df) = cluster
  pi.df$Cluster = 1:nrow(pi.df)

  pi.df
}



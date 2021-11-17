is_reasonable_clonal_cluster = function(x, cluster, scale = 3) 
{
  if(x$Kbeta == 1) return(TRUE)
  if(cluster == 'Tail') return(FALSE)
  
  Clonal_peak = x$Clusters %>% filter(type == "Mean", cluster == !!cluster) %>% pull(fit.value)
  
  # Beta peak above 0.5 is generally aneuploidy with VAF
  if(Clonal_peak > 0.5) return(FALSE)
  
  Clonal_size = x$Clusters %>% filter(type == "Mixing proportion", cluster == !!cluster) %>% pull(fit.value)
  
  subclonal = names(x$pi)
  subclonal = subclonal[subclonal != 'Tail']
  subclonal = subclonal[subclonal != cluster]
  
  subclonal_size = x$Clusters %>% filter(type == "Mixing proportion", cluster %in% subclonal) %>% pull(fit.value)
  
  if(any(Clonal_size * scale< subclonal_size)) return(FALSE)
  
  return(TRUE)
}

# med_variance = x$Clusters %>% filter(type == 'Variance') %>% pull(fit.value) %>% median(na.rm = TRUE)



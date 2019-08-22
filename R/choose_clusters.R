

#' Filter MOBSTER output clusters.
#' 
#' @description This function can filter out the clusters computed by MOBSTER
#' based on two criteria: the mixing proportion value, and the number of mutations
#' assigned. For both criteria a scalar should be given as input. The return object
#' will contain only the clusters that pass both filters. If any cluster is dropped
#' the latent variables are re-computed, as well as the clustering assignments and the
#' mixing proportions (all mutations will be still assigned after clusters' removal).
#'
#' @param x A MOBSTER fit object.
#' @param pi_cutoff The cutoff on the mixing proportions, default is 2%.
#' @param N_cutoff  The cutoff on the number of mutations assigned to a cluster, default is 10.
#'
#' @return A MOBSTER fit object where clusters are larger than \code{pi_cutoff} and contain
#' at least \code{N_cutoff}. If no such cluster exists an error is generated.
#' 
#' @export
#'
#' @examples
#' TODO
choose_clusters = function(x, 
                           pi_cutoff = 0.02,
                           N_cutoff = 10)
{
  stopifnot(inherits(x, "dbpmm"))
  
  pio::pioTit(paste0("Selecting MOBSTER clusters (F1,2-heuristic)."))
  pio::pioStr("\nF1.       Cluster size (proportion) > ", pi_cutoff)
  pio::pioStr("\nF2.        Cluster size (mutations) > ", N_cutoff, '\n\n')

  # Cluster size - remove clusters smaller than pi_cutoff, or less than N_cutoff mutations
  pass_clusters_picutoff = x$pi > pi_cutoff
  pass_clusters_Ncutoff = x$N.k > N_cutoff
  
  tab_clusters = tibble(cluster = names(x$pi), pi = x$pi, N = x$N.k)
  tab_clusters$F1 = pass_clusters_picutoff
  tab_clusters$F2 = pass_clusters_Ncutoff
  
  # CLusters that will be cancelled
  clusters_to_cancel = !pass_clusters_picutoff | !pass_clusters_Ncutoff 
  
  # Report to screen
  tab_clusters$remove = clusters_to_cancel
  
  print(tab_clusters %>% arrange(desc(remove)))
  
  # Report easy cases
  if(sum(tab_clusters$remove) == nrow(tab_clusters)) stop("All clusters filtered out, this seems an error.")
  if(sum(tab_clusters$remove) == 0) return(x)
  
  # Clusters that remain, which we remove from a copy of x
  remaining_clusters = tab_clusters %>% filter(!remove) %>% pull(cluster)
  remaining_beta_clusters = remaining_clusters[grepl('C', remaining_clusters)]
  
  # We modify a copy of the input
  y = x
  
  # Remove from the clusters table any information regarding old clusters,
  # preserving the tail's proportion as a valid entry
  y$Clusters = y$Clusters %>% filter(cluster %in% remaining_clusters)
  y$fit.tail = 'Tail' %in% remaining_clusters
  
  if(!('Tail' %in% remaining_clusters)) 
    y$Clusters = y$Clusters %>%
        bind_rows(
          tibble::tribble(
            ~cluster, ~type, ~fit.value, ~init.value,
            'Tail', "Mixing proportion", 0, NA
          )
        )
    
  # Update counts for Beta components
  y$Kbeta = length(remaining_beta_clusters)
  y$K = y$Kbeta + 1
  
  # Renormalize the mixing proportions -- adding 0% Tail if it's no longer fit
  y$pi = x$pi[remaining_clusters]
  y$pi = y$pi/sum(y$pi)
  if(!('Tail' %in% remaining_clusters)) y$pi['Tail'] = 0
    
  y = mobster:::.set_params_Pi(y, y$pi)
  
  # If we need to cancel a tail
  if(!('Tail' %in% remaining_clusters))
  {
    y$scale = y$shape = NA
    y = mobster:::.set_params_Pareto(y, y$scale, y$shape) # This should be useless
  }
  
  # For the Beta that remains, we need to subset their parameters
  if(length(remaining_beta_clusters) > 0)
  {
    y$a = y$a[remaining_beta_clusters]
    y$b = y$b[remaining_beta_clusters]
    y = mobster:::.set_params_Beta(y, y$a, y$b)
  }
  
  new_lv = mobster:::latent_vars(y)
  
  y$NLL = new_lv$NLL
  y$z_nk = new_lv$z_nk
  y$pdf.w = new_lv$pdf.w
  
  # Update clustering assignments and summary numbers
  y$data$cluster = mobster:::latent_vars_hard_assignments(new_lv)
  
  y$N.k = rep(0, y$K)
  names(y$N.k) = names(y$pi)
  obFreq = table(y$data$cluster)
  
  y$N.k[names(obFreq)] = obFreq
  
  # Update scores for model selection   
  y$scores = mobster:::latent_vars_scores(
    mobster:::latent_vars(y), # Extract latent variables
    y$K,
    y$fit.tail,
    y$data$cluster)
  
  pio::pioTit("Selected clusters")
  print(y)
  
  y
}

rename_Beta_cluster = function(x)
{
  params = x$Clusters %>%
    filter(type == 'Mean', cluster != 'Tail') %>%
    arrange(desc(fit.value)) %>%
    mutate(new.name = paste0('C', row_number())) 
  
  mapping = params$new.name
  names(mapping) = params$cluster
  
  mapping['Tail'] = 'Tail'
  mapping = mapping[c('Tail', params$cluster)]
  
  # Copy of x
  y = x

  # clusters
  y$data$cluster = mapping[y$data$cluster]
  
  # numbers from clustering
  names(y$N.k) = mapping[names(y$N.k)]
  
  # LV - forced assumed the order has been mantained
  colnames(y$z_nk) =  colnames(y$pdf.w) = mapping 
  
  # mixing
  names(y$pi) = mapping[names(y$pi)]
  y$pi = mobster:::.params_Pi(y)

  # Clusters table
  y$Clusters$cluster = mapping[y$Clusters$cluster]
  
  # Beta parmeters
  names(y$a) = mapping[names(y$a)]
  names(y$b) = mapping[names(y$b)]
  
  y
}

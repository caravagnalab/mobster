#' Filter MOBSTER output clusters.
#' 
#' @description This function can filter out the clusters computed by MOBSTER
#' based on two criteria: the mixing proportion value,  the number of mutations
#' assigned and the variance of the Beta clusters. 
#' 
#' For all criteria a scalar should be given as input. The return object
#' will contain only the clusters that pass all filters. If any cluster is dropped
#' the latent variables are re-computed, as well as the clustering assignments and the
#' mixing proportions (all mutations will be still assigned after clusters' removal).
#' 
#' @param x A MOBSTER fit object.
#' @param pi_cutoff The cutoff on the mixing proportions, default is 0.02.
#' @param N_cutoff The cutoff on the number of mutations assigned to a cluster, default is 10.  
#' @param Beta_variance_cutoff Minimum variance for a Beta peak.
#' @param verbose If outputs should be reported to screen or not, default is no.
#'
#' @return A MOBSTER fit object where clusters are larger than \code{pi_cutoff} and contain
#' at least \code{N_cutoff}. If no such cluster exists an error is generated.
#' 
#' @export
#' 
#' @examples
#' data('fit_example', package = 'mobster')
#'
#' # Does not change anything (no filter triggered)
#' choose_clusters(fit_example$best)
#' 
#' # Remove one Beta component because it has less than 100 points (renders the fit very poor)
#' choose_clusters(fit_example$best, N_cutoff = 100)
choose_clusters = function(x, 
                           pi_cutoff = 0.02,
                           N_cutoff = 10,
                           Beta_variance_cutoff = 0.0001,
                           verbose = FALSE)
{
  mobster:::is_mobster_fit(x)
  
  if(verbose)
  {
    cli::cli_rule("Filtering MOBSTER clusters")
    cli::cli_text(mobster:::m_inf(paste0("pi > {.field {", pi_cutoff, '}} (mixing proportions)')))
    cli::cli_text(mobster:::m_inf(paste0(" n > {.field {", N_cutoff, '}} (mutations per cluster)')))
    cli::cli_text(mobster:::m_inf(paste0(" v > {.field {", Beta_variance_cutoff, '}} (minimum Betta variance)')))
  }
  
  # Cluster size - remove clusters smaller than pi_cutoff, or less than N_cutoff mutations
  pi = mobster:::.params_Pi(x)
  varb = x$Clusters %>% dplyr::filter(type == 'Variance', cluster != 'Tail')
  
  pass_clusters_picutoff = pi > pi_cutoff
  pass_clusters_Ncutoff = x$N.k > N_cutoff
  pass_clusters_Bvcutoff = pio:::nmfy(varb$cluster, varb$fit.value > Beta_variance_cutoff)
  
  tab_clusters = tibble::tibble(cluster = names(pi), pi = pi, N = x$N.k[names(pi)])
  tab_clusters$F1 = pass_clusters_picutoff
  tab_clusters$F2 = pass_clusters_Ncutoff
  tab_clusters = tab_clusters %>% 
    dplyr::left_join(
      data.frame(stringsAsFactors = FALSE, 
                 cluster = names(pass_clusters_Bvcutoff), 
                 Beta_var = varb$fit.value,
                 F3 = pass_clusters_Bvcutoff),
      by = 'cluster'
      ) %>%
    dplyr::select(cluster, pi, N, Beta_var, starts_with('F'))
  
  tab_clusters$retain = apply(tab_clusters, 1, function(x) all(as.logical(x[5:length(x)]), na.rm = TRUE))
  
  # CLusters that will be cancelled
  clusters_to_cancel = tab_clusters$cluster[which(!tab_clusters$retain)]
    
  # Report to screen
  if(verbose) print(tab_clusters %>% arrange(desc(retain)))
  
  # Report easy cases
  if(sum(!tab_clusters$retain) == nrow(tab_clusters)) stop("All clusters filtered out, this seems an error.")
  if(sum(tab_clusters$retain) == 0) return(x)
  
  # Clusters that remain, which we remove from a copy of x
  remaining_clusters = tab_clusters$cluster[which(tab_clusters$retain)]
  remaining_beta_clusters = remaining_clusters[grepl('C', remaining_clusters)]
  
  # We modify a copy of the input
  y = x
  
  # Remove from the clusters table any information regarding old clusters,
  # preserving the tail's proportion as a valid entry
  y$Clusters = y$Clusters %>% dplyr::filter(cluster %in% remaining_clusters)
  y$fit.tail = 'Tail' %in% remaining_clusters
  
  if(!(y$fit.tail)) 
    y$Clusters = y$Clusters %>%
    dplyr::bind_rows(
      tibble::tribble(
        ~cluster, ~type, ~fit.value, ~init.value,
        'Tail', "Mixing proportion", 0, NA
      )
    )
  
  # Update counts for Beta components
  y$Kbeta = length(remaining_beta_clusters)
  y$K = y$Kbeta + 1
  
  # Renormalize the mixing proportions -- adding 0% Tail if it's no longer fit
  # and ensure that the ordering is 
  y$pi = x$pi[remaining_clusters]
  y$pi = y$pi/sum(y$pi)
  if(!(y$fit.tail)) 
  {
    y$pi['Tail'] = 0
    y$pi = y$pi[c('Tail', remaining_clusters)]
  }
  
  y = mobster:::.set_params_Pi(y, y$pi)
  
  # If we need to cancel a tail
  if(!(y$fit.tail)) 
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
  
  if(verbose)
  {
    pio::pioTit("Selected clusters")
    print(y)
  }
  
  y
}

# Rename Beta clusters so that C1 is the one with highest mean etc.
rename_Beta_clusters = function(x)
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
  y$data$cluster = mapping[y$data$cluster] %>% as.vector
  
  # numbers from clustering
  names(y$N.k) = mapping[names(y$N.k)] %>% as.vector
  
  # LV - forced assumed the order has been mantained
  colnames(y$z_nk) =  mapping[colnames(y$z_nk)]
  colnames(y$pdf.w) = mapping[colnames(y$pdf.w)]
  names(colnames(y$z_nk)) = names(colnames(y$pdf.w)) = NULL 
  
  # Clusters table
  y$Clusters$cluster = mapping[y$Clusters$cluster]  %>% as.vector
  
  # mixing
  # names(y$pi) = mapping[names(y$pi)]  %>% as.vector
  y$pi = .params_Pi(y)
  
  # Beta parmeters
  names(y$a) = mapping[names(y$a)]  %>% as.vector
  names(y$b) = mapping[names(y$b)]  %>% as.vector
  
  y
}
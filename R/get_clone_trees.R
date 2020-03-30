get_clone_trees = function(x, ...)
{
  if(!(c('driver', 'gene') %in% colnames(x$data)))
    stop("Your data should have a logical 'driver' and 'gene' column to annotate driver events, cannot build a ctree otherwise.")
  
  mobster:::is_mobster_fit(x)
  
  # Patient ID
  patientID = ifelse(is.null(x$description), "MOBSTER dataset", x$description)
  patientID = gsub(pattern = ' ', replacement = '_', patientID)
  
  # Get clusters table: cluster and fit
  cluster_table = x$Clusters %>%
    dplyr::filter(cluster != 'Tail', type == 'Mean') %>%
    dplyr::select(cluster, fit.value) %>%
    rename(R1 = fit.value)
  
  # Cluster size
  cluster_table$nMuts = x$N.k[cluster_table$cluster]
  
  # Clonality status - maximum fit is the clonal
  cluster_table$is.clonal = FALSE
  cluster_table$is.clonal[which.max(cluster_table$R1)] = TRUE
  
  # Detect presence/ absence of drivers
  drivers_collapse = x$data %>%
    dplyr::filter(driver) %>%
    pull(cluster) %>%
    unique
  
  cluster_table$is.driver = FALSE
  cluster_table$is.driver[which(cluster_table$cluster %in% drivers_collapse)] = TRUE
  
  # Create drivers table
  drivers_table = x$data %>%
    as_tibble() %>%
    dplyr::filter(driver) %>%
    dplyr::rename(variantID = gene, is.driver = driver) %>%
    dplyr::mutate(patientID = patientID, R1 = VAF) 
  
  drivers_table$is.clonal = FALSE
  drivers_table$is.clonal[which(
    drivers_table$cluster == cluster_table %>% dplyr::filter(is.clonal) %>% dplyr::pull(cluster)
  )] = TRUE
  
  drivers_table %>%  dplyr::select(patientID, variantID, is.driver, is.clonal, cluster, R1, dplyr::everything())
    
  trees = ctree::ctrees(
    CCF_clusters =  cluster_table,
    drivers = drivers_table,
    samples = 'R1',
    patient = patientID,
    ...
  )

  return(trees)  
}
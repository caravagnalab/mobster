#' Return clone trees from the fit.
#' 
#' @description This function uses the output fit of MOBSTER
#' to create a call to \code{ctree} (\url{https://caravagn.github.io/ctree/}),
#' a package to create clone trees for cancer evolution models.
#' 
#' Creation of a clone tree requires annotations that are not usually 
#' necessary for just a plain MOBSTER analyses. These annotations report the status of \code{driver} 
#' and \code{gene} for each one of the input datapoints, and should
#' be part of data given in input for MOBSTER (so they should be in \code{x$data}). 
#' 
#' MOBSTER clusters are only used if the come from a Beta distribtutions; that is
#' the tail is removed. The clonal cluster is estimated from the cluster with the highest parameter
#' value for the Beta peak.
#' 
#' The output is the result of calling the constructor \code{ctree::cetrees}
#' on the input clustering results \code{x}.
#'
#' @param x A MOBSTER fit.
#' @param ... Extra parameters passed to the constructor \code{ctree::cetrees}, which
#' affect the sampling of the trees.
#'
#' @return The output of the constructor \code{ctree::cetrees}.
#' 
#' @export
#'
#' @examples
#' # We take one of the released datasets
#' x = mobster::PD4120a_breast_sample$best
#' 
#' # Genes are already annotated
#' head(x$data$gene)
#' 
#' # Drivers we add, just taken at random from the exonic mutations.
#' # Just print: x$data %>% filter(region == 'exonic')
#' # SETD2 could really be a driver.
#' require(dplyr)
#' genes_list = c('ARHGAP31', 'ABCC5', 'SETD2')
#' x$data = x$data %>% 
#'   dplyr::mutate(
#'     driver = ifelse(region == 'exonic' & gene %in% genes_list, TRUE, FALSE)
#'     )
#' 
#' # Get the trees
#' trees = get_clone_trees(x)
#' 
#' # Print and plot the first tree (top rank)
#' ctree:::print.ctree(trees[[1]])
#' ctree::plot.ctree(trees[[1]])
#' ctree::plot_CCF_clusters(trees[[1]])
#' ctree::plot_icon(trees[[1]])
#' ctree::plot_clone_size(trees[[1]])
get_clone_trees = function(x, ...)
{
  if(!all(c('driver', 'gene') %in% colnames(x$data)))
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
  
  drivers_table = drivers_table %>%  dplyr::select(patientID, variantID, is.driver, is.clonal, cluster, R1, dplyr::everything())
    
  trees = ctree::ctrees(
    CCF_clusters =  cluster_table,
    drivers = drivers_table,
    samples = 'R1',
    patient = patientID,
    ...
  )

  return(trees)  
}

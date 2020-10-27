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
#' @importFrom ctree ctrees
#'
#' @examples
#' # We take one of the released datasets
#' x = mobster::LU4_lung_sample$best
#' 
#' # Annotate some random mutation as driver, we need that to build the trees with ctree
#' x$data$is_driver = FALSE
#' x$data$is_driver[1:3] = TRUE
#' 
#' x$data$driver_label = ""
#' x$data$driver_label[1] = "Fake_driver_1"
#' x$data$driver_label[2] = "Fake_driver_2"
#' x$data$driver_label[3] = "Fake_driver_3"
#' 
#' # Get the trees
#' trees = get_clone_trees(x)
#' 
#' # Print and plot the first tree (top rank)
#' library(ctree)
#' ctree::print.ctree(trees[[1]])
#' ctree::plot.ctree(trees[[1]])
#' ctree::plot_CCF_clusters(trees[[1]])
#' ctree::plot_icon(trees[[1]])
#' ctree::plot_clone_size(trees[[1]])
get_clone_trees = function(x, ...)
{
  mobster:::is_mobster_fit(x)
  
  if(!mobster:::has_drivers_annotated(x))
    stop("Your data should have driver events annotated, cannot use 'ctree' otherwise.")
  
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
    dplyr::filter(is_driver) %>%
    pull(cluster) %>%
    unique
  
  if(length(drivers_collapse) == 0) 
    stop("Your data should have driver events annotated, cannot use 'ctree' otherwise.")
  
  cluster_table$is.driver = FALSE
  cluster_table$is.driver[which(cluster_table$cluster %in% drivers_collapse)] = TRUE
  
  # Create drivers table
  drivers_table = x$data %>%
    as_tibble() %>%
    dplyr::filter(is_driver) %>%
    dplyr::rename(variantID = driver_label, is.driver = is_driver) %>%
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

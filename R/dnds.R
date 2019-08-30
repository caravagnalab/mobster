#' Run a dN/dS analysis on Mobster clusters
#'
#' @description This function takes a MOBSTER fit and runs dndscv (https://github.com/im3sanger/dndscv) to calculate dN/dS
#' values per cluster. It computes global dN/dS and per gene dN/dS values and makes a plot. dN/dS values are computed with
#' the best fitting MOBSTER model.
#'
#' @param x A MOBSTER fit object
#' @param gene_list An optional vector of gene names to infer dN/dS values,
#' default is to use the whole exome
#' @param colors If provided, these colours will be used for each cluster.
#' If a subset of colours is provided, palette Set1 from \code{RColorBrewer} is used.
#' By default the tail colour is provided as 'gainsboro'.
#' @param refdb The genome referene to use, default is to use hg19. Other references are available from
#' https://github.com/im3sanger/dndscv_data
#' @param dndscv_plot What of the dndscv scores should be visualized in a plot, by default (`wall`) is
#' only the global dnds value.
#' @param mapping The groups used to compute this statistics are defined by this variable. If
#' `mapping  = c(`A` = 'G1', `B` = 'G1', `C` = 'G2')`, then mutations from clusters `A` and 
#' `B` will be pooled into one group (`G1`), while mutations from cluster `C` will constitute
#' a group themselves. By default, each cluster in the fit is a own group. 
#'
#' @return The fit object is a list with the summary table and the observation counts reported
#' by package \code{dndscv}, together with a \code{ggplot} plot.
#'
#' @export
#'
dnds <- function(x,
                 gene_list = NULL,
                 colors = c(`Tail` = 'gray'),
                 refdb = "hg19",
                 dndscv_plot = 'wall',
                 mapping = pio:::nmfy(mobster:::.get_clusters_labels(x),
                                      mobster:::.get_clusters_labels(x))
                 )
{
  # Check(s): dndscv installationand mobster fit
  check_dnds_package()
  mobster:::is_mobster_fit(x)
  
  # Getter -- checks for the mapping correctness and apply it
  dnds_input = get_dnds_input(x, mapping, refdb)
  cl = dnds_input$clusters
  clusters = dnds_input$clusters_labels
  
  pio::pioTit("Running dndscv")
  
  labels_outputs = unique(mapping)
  
  result_fit = wrapper_dndsfit(clusters = cl,
                               groups = labels_outputs,
                               gene_list,
                               mode = 'Mapping')
  
  pio::pioStr("Generating ouptut plot", '\n')
  
  plot_results = wrapper_plot(
    result_fit,
    mode = result_fit$dndstable$run[1],
    gene_list,
    dndscv_plot,
    mask_colors = TRUE
  )
  
  results <- list(
    dnds_summary = result_fit$dndstable %>% as_tibble(),
    dndscv_table = result_fit$dndscvtable  %>% as_tibble(),
    plot = plot_results
  )
 
  return(results)
}



dnds_multifits <- function(x,
                           gene_list = NULL,
                           colors = c(`Tail` = 'gray'),
                           refdb = "hg19",
                           dndscv_plot = 'wall',
                           mapping = pio:::nmfy(
                             mobster:::.get_clusters_labels(x),
                             mobster:::.get_clusters_labels(x)
                             )
                           )
{
  
  stop("TODO  - implement this")
  
  return()
}
#' Run a dN/dS analysis on MOBSTER clusters.
#'
#' @description This function takes a MOBSTER fit and runs `dndscv` (https://github.com/im3sanger/dndscv) to calculate dN/dS
#' values per cluster. It computes global dN/dS and per gene dN/dS values and makes a plot. dN/dS values are computed with
#' the best fitting MOBSTER model.
#'
#' @param x A MOBSTER fit object.
#' @param mapping The groups used to compute this statistics are defined by this variable. If
#' `mapping  = c(`A` = 'G1', `B` = 'G1', `C` = 'G2')`, then mutations from clusters `A` and 
#' `B` will be pooled into one group (`G1`), while mutations from cluster `C` will constitute
#' a group themselves. By default, with `mapping = NULL`, each cluster is a group. 
#' @param gene_list An optional vector of gene names to infer dN/dS values,
#' default (`NULL`) is to use \code{dndscv} default (whole-exome. This package provides lists
#' genes that can be used for this value (essential genes, cancer genes, etc.); see package data. 
#' @param colors If provided, these colours will be used for each cluster.
#' If a subset of colours is provided, palette Set1 from \code{RColorBrewer} is used.
#' By default the tail colour is provided as 'gainsboro'.
#' @param refdb The genome referene to use, default is to use hg19. Other references are available from
#' https://github.com/im3sanger/dndscv_data
#' @param dndscv_plot What of the dndscv scores should be visualized in a plot, by default all the statistcs
#' are reported. One can use `dndscv_plot = wall` to get only the global dnds value.
#' @param ... Extra parameters forwarded to \code{dndscv}.
#'
#' @return The fit object is a list with the summary table and the observation counts reported
#' by package \code{dndscv}, together with a \code{ggplot} plot for the results.
#'
#' @importFrom dndscv dndscv
#' @export
#' 
#' @examples
#' 
#' # Example run with real data
#' data('LUFF76_lung_sample', package = 'mobster')
#' 
#' clusters = Clusters(LUFF76_lung_sample$best)
#' 
#' dnds_stats = dnds(clusters, gene_list = NULL)
dnds <- function(x,
                 mapping = NULL,
                 gene_list = NULL,
                 colors = c(`Tail` = 'gray'),
                 refdb = "hg19",
                 dndscv_plot = c('wall', 'wmis', 'wnon', 'wspl', 'wtru'),
                 ...
)
{
  # Check(s): dndscv installationand mobster fit
  mobster:::crash_ifnotinstalled(c('dndscv'))
  
  if(!is.data.frame(x))
  {
    mobster:::is_mobster_fit(x)
    x = mobster::Clusters(x)
  }
  
  # Getter -- checks for the mapping correctness and apply it
  dnds_input = mobster:::get_dnds_input(x, mapping, refdb, gene_list)
  clusters = unique(dnds_input$dnds_group)
  
  cli::cli_h1("Running dndscv")

  result_fit = mobster:::wrapper_dndsfit(clusters = dnds_input,
                               groups = clusters,
                               gene_list,
                               mode = 'Mapping',
                               refdb = refdb,
                               ...)
  
  cli::cli_rule(left = "dndscv results", right = paste0(dndscv_plot, collapse = ', '))
  pio::pioDisp(result_fit$dndstable %>% filter(name %in% dndscv_plot))
  
  plot_results = mobster:::wrapper_plot(
    result_fit,
    mode = result_fit$dndstable$run[1],
    gene_list,
    dndscv_plot,
    colors,
    mask_colors = TRUE
  )
  
  results <- list(
    dnds_summary = result_fit$dndstable %>% as_tibble(),
    dndscv_table = result_fit$dndscvtable  %>% as_tibble(),
    dndscv_output = result_fit$dndscvout,
    plot = plot_results
  )
  
  return(results)
}





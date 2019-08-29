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
#' @param plot Which dndscv result should be plot, by default is the global dnds `wall`
#' 
#' @return the fit object is returned with additional dnds list to the best fit model.
#'
#'
#'
#' @export
dnds <- function(x, 
                 gene_list = NULL, 
                 colors = c(`Tail` = 'gray'),
                 refdb = "hg19",
                 dndscv_plot = 'wall',
                 clonal_cluster = 'C1')
  {
  
  # check dndscv
  package = installed.packages() %>% 
    as_tibble() %>%
    filter(Package == 'dndscv')
  
  if(nrow(package) == 0)
  {
    stop(
      "dndscv is not installed, you should install it.\nSee https://github.com/im3sanger/dndscv"
    )
  }
  else{
    pio::pioTit("Running dndscv wrapper from MOBSTER")
    pio::pioStr("dndscv version", package$Version)
  }
  
  # Check input: a mobster fit
  mobster:::is_mobster_fit(x)
  
  fit = x
  
  # Check input: required columns
  required_columns = c('chr', 'from', 'ref', 'alt')
  if(!all(required_columns %in% names(fit$data)))
    stop(
      "Required columns are missing: ", paste(required_columns, collapse = ', '), '\nCannot compute dnds...'
    )
  
  # Clustering assignments are used to find coding mutations
  cl <- Clusters(fit) %>%
    dplyr::mutate(dummysample = "sample") %>%
    dplyr::select(dummysample, chr, from, ref, alt, everything())
  
  # Checkings for the reference
  if (refdb == "hg19" & stringr::str_detect(cl$chr[1], "chr")){
    message("Removing chr from chromosome names for hg19 reference compatability")
    cl$chr <- stringr::str_sub(cl$chr, 4)
  }
  
  clusters = unique(cl$cluster)
  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Run by Cluster
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pio::pioTit("Running dndscv by cluster")
  
  by_cluster = wrapper_dndsfit(clusters = cl, groups = clusters, gene_list, mode = 'By cluster')
    
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Run Clonal vs Subclonal
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # check input
  if(!(clonal_cluster %in% clusters)) {
    stop("Putative clonal cluster", clonal_cluster, 'does not exist?') 
  }
  
  pio::pioTit("Running dndscv clonal versus subclonal")
  
  subclonal_clusters = setdiff(clusters, clonal_cluster)
  label_clonal = paste0('Clonal (', clonal_cluster,')')
  label_subclonal = paste0('Subclonal (', paste0(subclonal_clusters, collapse = ', '),')')
  
  by_cln_subcl = wrapper_dndsfit(
    clusters = cl %>%
      mutate(
        cluster = ifelse(cluster == clonal_cluster, label_clonal, label_subclonal)
      ), 
    groups = c(label_clonal, label_subclonal), 
    gene_list, 
    mode = 'By clonality')
  
  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Run Tail vs Non-Tail
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  by_tail = NULL

  # check input
  if('Tail' %in% clusters) 
  {
    pio::pioTit("Running dndscv by tail status")
  
    non_tail = setdiff(clusters, 'Tail')
    label_nontail = paste0('Non-tail (', paste0(non_tail, collapse = ', '), ')')

    by_tail = wrapper_dndsfit(
      clusters = cl %>%
        mutate(
          cluster = ifelse(cluster == 'Tail', 'Tail', label_nontail)
        ), 
      groups = c(label_nontail, "Tail"), 
      gene_list, 
      mode = 'By tail')
  
  }
  else
  {
    message("\n> Will not run by tail status because there is no tail in this fit.")
  }
    
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plots all of them
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pby_cluster = wrapper_plot(
    by_cluster,
    mode = by_cluster$dndstable$run[1],
    gene_list,
    dndscv_plot,
    mask_colors = TRUE
  )
    
  pby_clon = wrapper_plot(
    by_cln_subcl,
    mode = by_cln_subcl$dndstable$run[1],
    gene_list,
    dndscv_plot,
    mask_colors = FALSE
  )
  
  pby_tail = ggplot() + geom_blank()
  
  if('Tail' %in% clusters) 
    pby_tail = wrapper_plot(
      by_tail,
      mode = by_tail$dndstable$run[1],
      gene_list,
      dndscv_plot,
      mask_colors = FALSE
    )
  
  figure = cowplot::plot_grid(
    pby_cluster,
    pby_clon,
    pby_tail,
    ncol = 3, 
    nrow = 1,
    align = 'h'
  )
  
  # Assembled results
  dndstable = bind_rows(by_cluster$dndstable, by_cln_subcl$dndstable) %>% as_tibble()
  dndscvtable = bind_rows(by_cluster$dndscvtable, by_cln_subcl$dndscvtable) %>% as_tibble()
  
  results <- list(dnds_summary = dndstable, dndscv_table = dndscvtable, 
                  plot_by_cluster = pby_cluster, plot_by_clonality = pby_clon)
  return(fit)
}

wrapper_dndsfit = function(clusters, groups, gene_list, mode)
{
  globaldndstable = dndscvtable = NULL
  for (i in groups)
  {
    pio::pioStr("\nCalculating dN/ds values:", i, "@", mode, '\n')
    
    dndsout <- clusters %>%
      dplyr::filter(cluster == i) %>%
      dndscv::dndscv(., gene_list = gene_list)
    
    globaldndstable <- dplyr::bind_rows(globaldndstable, dndsout$globaldnds %>%
                                          mutate(group = i, run = mode))
    dndscvtable <- dplyr::bind_rows(dndscvtable, dndsout$sel_cv %>%
                                      mutate(group = i, run = mode))
  }
  
  return(list(dndstable = globaldndstable, dndscvtable = dndscvtable))
}

wrapper_plot = function(results, mode, gene_list, dndscv_plot, mask_colors = FALSE)
{
  # Plotting
  gene_list_label = ifelse(
    is.null(gene_list),
    paste0('No input genes (default dndscv)'),
    paste0(length(gene_list), ' input genes')
  )

  dndsplot <- 
    results$dndstable %>% 
    dplyr::filter(name %in% dndscv_plot) %>% 
    ggplot2::ggplot(ggplot2::aes(x = group, y = mle, ymin = cilow, ymax = cihigh)) +
    mobster:::my_ggplot_theme() +
    facet_wrap(~ name, ncol = 1, scales = 'free_y') +
    ggplot2::xlab("") +
    ggplot2::ylab("dN/dS") +
    labs(
      title = paste0("dN/dS values ", mode),
      subtitle = gene_list_label
    ) +
    ggplot2::geom_hline(yintercept = 1.0, lty = 2, size = .3) +
    guides(color = F, fill = F)
  
  
  # Add or not the colours...
  if (mask_colors) 
    {
      dndsplot = dndsplot + geom_pointrange(aes(color = group))
      dndsplot = suppressMessages(mobster:::add_color_pl(x, dndsplot, colors))
    }
    else
    {
      dndsplot = dndsplot + geom_pointrange(color = 'black')
    }
    
    
  dndsplot
}

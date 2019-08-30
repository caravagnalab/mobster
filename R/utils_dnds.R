# Check package
check_dnds_package = function()
{
  package = installed.packages() %>%
    as_tibble() %>%
    filter(Package == 'dndscv')
  
  if (nrow(package) == 0)
  {
    stop(
      "dndscv is not installed, you should install it.\nSee https://github.com/im3sanger/dndscv"
    )
  }
  else{
    pio::pioTit("Running dndscv wrapper from MOBSTER")
    pio::pioStr("dndscv version", package$Version, '\n')
  }
}

# Return the input for dndscv, checking that parameters make sense
get_dnds_input = function(x, mapping, refdb)
{
  # Check input: required columns
  required_columns = c('chr', 'from', 'ref', 'alt')
  if (!all(required_columns %in% names(x$data)))
    stop(
      "Required columns are missing: ",
      paste(required_columns, collapse = ', '),
      "\nYour columns are: ", paste(colnames(x$data), collapse = ', '),
      '\nCannot compute dnds if your data does not have the required columns...'
    )
  
  # Clustering assignments are used to find coding mutations
  cl <- mobster::Clusters(x) %>%
    dplyr::mutate(dummysample = "sample") %>%
    dplyr::select(dummysample, chr, from, ref, alt, everything())
  
  pio::pioStr("Using these data", '\n')
  pio::pioDisp(cl)
  
  # Checkings for the reference
  if (refdb == "hg19" & stringr::str_detect(cl$chr[1], "chr")) {
    message("Removing chr from chromosome names for hg19 reference compatability")
    cl$chr <- stringr::str_sub(cl$chr, 4)
  }
  
  clusters = unique(cl$cluster)
  
  # Check on the mapping
  # m: Clusters -> Groups
  if (!all(clusters %in% names(mapping))) {
    missing = clusters[!(clusters %in% names(mapping))]
    stop(
      'Mapping is not exauhstive, cannot use it.\n',
      "There should be one entry for each one of ", paste(clusters, collapse = ', '), ' ~ ',
      'Cluster(s) ', paste(missing,  collapse = ', '), " are missing!"
    )
  }

  # Apply mapping
  pio::pioStr("Using the following mapping", '\n')
  print(mapping)
  
  cl$dnds_group = mapping[cl$cluster]
  
  return(cl)
}

# Fits via dndscv
wrapper_dndsfit = function(clusters, groups, gene_list, mode)
{
  globaldndstable = dndscvtable = NULL
  for (i in groups)
  {
    pio::pioStr("\ndndscv @ dnds_group ", i, '\n')
    
    dndsout = NULL
    tryCatch({
      dndsout <- clusters %>%
        dplyr::filter(dnds_group == i) %>%
        dndscv::dndscv(., gene_list = gene_list)
    },
    error = function(e)
    {
      # Intercepted error
      message('Intercepted error from dndscv')
      message(e)
      
      dndsout = NULL
    })
  
    
    if(!is.null(dndsout))
    {
      globaldndstable <- dplyr::bind_rows(globaldndstable, dndsout$globaldnds %>%
                                          mutate(dnds_group = i))
    dndscvtable <- dplyr::bind_rows(dndscvtable, dndsout$sel_cv %>%
                                      mutate(dnds_group = i))
    }
  }
  
  return(list(dndstable = globaldndstable, dndscvtable = dndscvtable))
}

# Plotting function for dndscv results
wrapper_plot = function(results, mode, gene_list, dndscv_plot, colors, mask_colors = FALSE)
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
    ggplot2::ggplot(ggplot2::aes(x = dnds_group, y = mle, ymin = cilow, ymax = cihigh)) +
    mobster:::my_ggplot_theme() +
    facet_wrap(~ name, ncol = 1, scales = 'free_y') +
    ggplot2::xlab("") +
    ggplot2::ylab("dN/dS") +
    labs(
      title = paste0("dN/dS values via dndscv"),
      subtitle = gene_list_label
    ) +
    ggplot2::geom_hline(yintercept = 1.0, lty = 2, size = .3) +
    guides(color = F, fill = F)
  
  
  # Add or not the colours...
  if (mask_colors) 
  {
    dndsplot = dndsplot + geom_pointrange(aes(color = dnds_group))
    dndsplot = suppressMessages(mobster:::add_color_pl(x, dndsplot, colors))
  }
  else
  {
    dndsplot = dndsplot + geom_pointrange(color = 'black')
  }
  
  
  dndsplot
}

tail_non_tail_mapping = function(n=10)
{
  c(
    `Tail` = "Tail",
    pio:::nmfy(paste0("C", 1:n), rep("Non-tail", n))
  )
}

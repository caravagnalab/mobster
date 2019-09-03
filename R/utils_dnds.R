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
    pio::pioTit("MOBSTER wrapper for dndscv version", package$Version)
  }
}

# Return the input for dndscv, checking that parameters make sense
get_dnds_input = function(x, mapping, refdb, gene_list)
{
  # Check input: required columns
  required_columns = c('chr', 'from', 'ref', 'alt')
  if (!all(required_columns %in% names(x)))
    stop(
      "Required columns are missing: ",
      paste(required_columns, collapse = ', '),
      "\nYour columns are: ", paste(colnames(x), collapse = ', '),
      '\nCannot compute dnds if your data does not have the required columns...'
    )
  
  # Clustering assignments are used to find coding mutations
  if(!('sample' %in% colnames(x)))
 {   
    message("Missing sample column.\n" ,
            "> Assuming these are mutations from a single patient, adding dummy sample id.\n",
            "> If this is not the case, label each mutation with a sample id.")
    
    x <- x %>%
      dplyr::mutate(dummysample = "sample") %>%
      dplyr::select(dummysample, chr, from, ref, alt, everything())
  }
  else
    x <- x %>%
      dplyr::select(sample, chr, from, ref, alt, everything())
  
  
  pio::pioStr("\n  Mutations ", nrow(x))
  pio::pioStr("\n      Genes ", length(gene_list))
  pio::pioStr("\n   Clusters ", length(unique(x$clusters)))
  pio::pioStr("\nDnds groups ", length(unique(mapping)), '\n')
  
  # pio::pioDisp(x)
  
  # Checkings for the reference
  if (refdb == "hg19") {
    message("[refdb = hg19] \n" ,
            "> Removing chr from chromosome names for hg19 reference compatability\n")
    
    x$chr <- gsub(pattern = 'chr', replacement = '', x$chr)
  }
  
  clusters = unique(x$cluster)
  
  pio::pioStr("Mapping clusters to dnds_groups\n")
  
  # Special mapping: identity
  if(all(is.null(mapping))) 
  {
    message("[mapping = NULL]\n", "> Creating mapping by cluster.")
    
    mapping = pio:::nmfy(clusters, clusters)
  }
  
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

  x$dnds_group = mapping[x$cluster]
  
  print(table(x$dnds_group))
  
  return(x)
}

# Fits via dndscv
wrapper_dndsfit = function(clusters, groups, gene_list, mode, ...)
{
  globaldndstable = dndscvtable = NULL
  for (i in groups)
  {
    pio::pioStr("\ndndscv @ dnds_group ", i, '\n')
    
    dndsout = NULL
    tryCatch({
      dndsout <- clusters %>%
        dplyr::filter(dnds_group == i) %>%
        dndscv::dndscv(., gene_list = gene_list,  ...)
    },
    error = function(e)
    {
      # Intercepted error
      cat(crayon::red('BEGIN Intercepted error from dndscv\n'))
      message(e)
      cat(crayon::red('\n'))
      
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
  counts = results$dndscvtable %>% 
    group_by(dnds_group) %>%
    summarise(wsyn = sum(n_syn), wmis = sum(n_mis), wnon = sum(n_non), wspl = sum(n_spl), wtru = sum(n_spl) + sum(n_non), wall = wtru + sum(n_syn) + sum(n_mis)) %>%
    reshape2::melt(id = 'dnds_group') %>%
    rename(name = variable, n = value) 
  
  syn_label = counts %>% filter(name == 'wsyn') %>%
    mutate(label = paste0(dnds_group, ' (n = ', n, ')')) %>%
    pull(label) %>%
    paste(collapse =', ')
  
  counts =  counts %>% filter(name != 'wsyn') %>%
    full_join(
      results$dndstable %>% 
        mutate(name = paste(name)) %>%
        select(name, mle, dnds_group)
      )
  
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
    facet_wrap(~ name, nrow = 1, scales = 'free_y') +
    ggplot2::xlab("") +
    ggplot2::ylab("dN/dS") +
    labs(
      title = paste0("dN/dS values via dndscv"),
      subtitle = gene_list_label,
      caption = paste('Synonimous mutations:', syn_label)
    ) +
    ggplot2::geom_hline(yintercept = 1.0, lty = 2, size = .3) +
    guides(color = F, fill = F)
  
  
  # Add or not the colours...
  if (mask_colors) 
  {
    dndsplot = dndsplot + geom_pointrange(aes(color = dnds_group))
    dndsplot = suppressMessages(mobster:::add_color_pl(unique(results$dndstable$dnds_group), dndsplot, colors))
  }
  else
  {
    dndsplot = dndsplot + geom_pointrange(color = 'black')
  }
  
  dndsplot +
    ggrepel::geom_label_repel(data = counts %>% mutate(n = paste0('n = ', n)),
                              inherit.aes = FALSE, aes(x = dnds_group, y = mle, label = n, color = dnds_group), size = 3)
}

tail_non_tail_mapping = function(n=10)
{
  c(
    `Tail` = "Tail",
    pio:::nmfy(paste0("C", 1:n), rep("Non-tail", n))
  )
}

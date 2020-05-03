# Return the input for dndscv, checking that parameters make sense
get_dnds_input = function(x, mapping, refdb, gene_list)
{
  # Check input: required columns
  required_columns = c('chr', 'from', 'ref', 'alt')
  if (!all(required_columns %in% names(x)))
    stop(
      "Required columns are missing: ",
      paste(required_columns, collapse = ', '),
      "\nYour columns are: ",
      paste(colnames(x), collapse = ', '),
      '\nCannot compute dnds if your data does not have the required columns...'
    )
  
  # Clustering assignments are used to find coding mutations
  if (!('sample' %in% colnames(x)))
  {
    message(
      "Missing 'sample' column, assuming mutations from a single patient (adding a sample label otherwise)."
    )
    
    x <- x %>%
      dplyr::mutate(dummysample = "sample") %>%
      dplyr::select(dummysample, chr, from, ref, alt, everything())
  }
  else
    x <- x %>%
      dplyr::select(sample, chr, from, ref, alt, everything())
  
  ngrp = ifelse(
    all(is.null(mapping)),
    "'by cluster'",
    length(unique(mapping))
  )
  
  gl = ifelse(
    all(is.null(gene_list)),
    "no genes (default dndscv)",
    paste0(length(gene_list), " genes (custom list)")
  )
  
  nm = nrow(x)
  ns = length(unique(x$sample))
  
  cli::cli_alert_info(
      "{.value {nm}} mutations; {.field {ngrp}} groups in {.value {ns}} samples, with {.field {gl}}."
    )
  
  # Checkings for the reference
  if (refdb == "hg19") {
    message(
      "[refdb = hg19] Removing chr from chromosome names for hg19 reference compatability"
    )
    
    x$chr <- gsub(pattern = 'chr', replacement = '', x$chr)
  }
  
  clusters = unique(x$cluster)
  
  # pio::pioStr("Mapping clusters to dnds_groups\n")
  
  # Special mapping: identity
  if (all(is.null(mapping)))
  {
    # message("[mapping = NULL]\n", "> Creating mapping by cluster.")
    
    mapping = pio:::nmfy(clusters, clusters)
  }
  
  # Check on the mapping
  # m: Clusters -> Groups
  if (!all(clusters %in% names(mapping))) {
    missing = clusters[!(clusters %in% names(mapping))]
    stop(
      'Mapping is not exauhstive, cannot use it.\n',
      "There should be one entry for each one of ",
      paste(clusters, collapse = ', '),
      ' ~ ',
      'Cluster(s) ',
      paste(missing,  collapse = ', '),
      " are missing!"
    )
  }
  
  x$dnds_group = mapping[x$cluster]
  
  print(table(x$dnds_group))
  
  return(x)
}

# Fits via dndscv
wrapper_dndsfit = function(clusters, groups, gene_list, mode, refdb, ...)
{
  globaldndstable = dndscvtable = NULL
  dndsoutputs = NULL
  for (i in groups)
  {
    cli::cli_h3(paste0("Group ", i))

    dndsout = NULL
    tryCatch({
      dndsout <- clusters %>%
        dplyr::filter(dnds_group == i) %>%
        dndscv::dndscv(., gene_list = gene_list, refdb = refdb, ...)
    },
    error = function(e)
    {
      # Intercepted error
      cat(crayon::red('dndscv error\n'))
      cat(crayon::blue(e))
      cat("\n")
      
      dndsout = NULL
    })
    
    
    if (!is.null(dndsout))
    {
      globaldndstable <-
        dplyr::bind_rows(globaldndstable,
                         dndsout$globaldnds %>%
                           mutate(
                             name = paste(name),
                             dnds_group = i))
      dndscvtable <-
        dplyr::bind_rows(dndscvtable, dndsout$sel_cv %>%
                           mutate(dnds_group = i))
    }
    
    dndsoutputs = append(dndsoutputs, list(dndsout))
  }
  
  names(dndsoutputs) = groups
  
  return(list(dndstable = globaldndstable, dndscvtable = dndscvtable, dndscvout = dndsoutputs))
}

# Plotting function for dndscv results
wrapper_plot = function(results,
                        mode,
                        gene_list,
                        dndscv_plot,
                        colors,
                        mask_colors = FALSE)
{
  counts = results$dndscvtable %>%
    group_by(dnds_group) %>%
    summarise(
      wsyn = sum(n_syn),
      wmis = sum(n_mis),
      wnon = sum(n_non),
      wspl = sum(n_spl),
      wtru = sum(n_spl) + sum(n_non),
      wall = wtru + sum(n_syn) + sum(n_mis)
    ) %>%
    reshape2::melt(id = 'dnds_group') %>%
    dplyr::rename(name = variable, n = value) %>%
    mutate(name = paste(name))
  
  syn_label = counts %>% filter(name == 'wsyn') %>%
    mutate(label = paste0(dnds_group, ' (n = ', n, ')')) %>%
    pull(label) %>%
    paste(collapse = ', ')
  
  counts = counts %>% 
    filter(name != 'wsyn') %>%
    filter(name %in% dndscv_plot) %>%
    full_join(
      results$dndstable %>%
        dplyr::filter(name %in% dndscv_plot) %>%
        mutate(name = paste(name)) %>%
        select(name, mle, dnds_group),
      by = c("dnds_group", "name")
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
    ggplot2::ggplot(ggplot2::aes(
      x = dnds_group,
      y = mle,
      ymin = cilow,
      ymax = cihigh
    )) +
    mobster:::my_ggplot_theme() +
    facet_wrap( ~ name, nrow = 1, scales = 'free_y') +
    ggplot2::xlab("") +
    ggplot2::ylab("dN/dS") +
    labs(
      title = paste0("dN/dS values via dndscv"),
      subtitle = gene_list_label,
      caption = paste('Synonymous mutations:', syn_label)
    ) +
    ggplot2::geom_hline(yintercept = 1.0,
                        lty = 2,
                        size = .3) +
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
    ggrepel::geom_label_repel(
      data = counts %>% mutate(n = paste0('n = ', n)),
      inherit.aes = FALSE,
      aes(
        x = dnds_group,
        y = mle,
        label = n,
        color = dnds_group
      ),
      size = 3
    )
}

tail_non_tail_mapping = function(n = 10)
{
  c(`Tail` = "Tail",
    pio:::nmfy(paste0("C", 1:n), rep("Non-tail", n)))
}

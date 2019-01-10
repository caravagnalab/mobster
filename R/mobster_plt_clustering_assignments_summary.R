#' Title
#'
#' @param fit_with_MOBSTER 
#' @param fit_without_MOBSTER 
#' @param CCF.cutoff 
#' @param title 
#' @param cex_legend 
#'
#' @return
#' @export
#'
#' @examples
mobster_plt_clustering_assignment_summary = function(
  fit_with_MOBSTER, 
  fit_without_MOBSTER, 
  CCF.cutoff = .70,
  title = fit_with_MOBSTER$description,
  cex_legend = 1)
{
  x = fit_with_MOBSTER
  y = fit_without_MOBSTER

  pio::pioTit("Plotting ...")
  
  pio::pioStr("Cutoff to binarize adjusted VAF", paste0(CCF.cutoff/2, " (CCF ", CCF.cutoff, ')' ))
  
  CCF.cutoff = CCF.cutoff/2

  sample = x$samples[1]
  k = length(x$samples)

# mutations = full_join(VAF_table(x), MOBSTER_clusters, by = 'id')

# M_clusters = MClusters(x)
# colnames(M_clusters)[1+1:k] = paste0("MOBSTER.", colnames(M_clusters)[1+1:k])
# 
# mutations = full_join(VAF_table(x), M_clusters, by = 'id')
# 
# mutations = VAF(x) %>% select(-id)
# colnames(mutations) = sub('.VAF',  '', colnames(mutations))
# 
# pheatmap::pheatmap(mutations)
# 
# bin.mutations = as.data.frame( %>% select(-id))
# rownames(bin.mutations) = VAF_table(x)$id
# bin.mutations[bin.mutations > 0] = 1



  suffix = paste0(" >", CCF.cutoff)
  
  # VAF data, to which we add the median VAF (for sorting)
  VAF_data = VAF_table(y) 
  VAF_data$VAF_median = apply(VAF_data[1+1:k], 1, mean) 
  
  # Binarize data from VAF
  is = function(x){ is.numeric(x) & x > CCF.cutoff} # We binarize to 1 everything above CCF.cutoff
  
  bin.mutations = VAF_data %>%
    select(id, ends_with('.VAF')) %>% 
    mutate_if(~ any(is(.x)),~ if_else(is(.x),1,0))
  
  # We count also how many =1 binarized mutations per sample are there
  bin.mutations$Above_cutoff = rowSums(bin.mutations[1+1:k])
  
  # We assemble mutations, with binary data as well
  mutations = full_join(bin.mutations, VAF_data, by = 'id', suffix = c(suffix, ""))
  
  # Prepare clustering assignments from MOBSTER, and Binomial data 
  M_clusters = MClusters(x)
  colnames(M_clusters)[1+1:k] = gsub(pattern = 'cluster.', replacement = 'MOBSTER_', colnames(M_clusters)[1+1:k])
  
  # paste0(colnames(M_clusters)[1+1:k], '.MOBSTER')
  B_clusters = BClusters(x) %>% rename(Binomial = cluster.Binomial)
  
  # All clustering outputs
  all_clusters = full_join(B_clusters, M_clusters, by = 'id') %>% 
    select(
      starts_with('Binomial'), 
      starts_with('MOBSTER'),
      'id')
  
  is.tail = apply(all_clusters, 1, function(w) any(w == 'Tail', na.rm = TRUE))
  
  all_clusters$tail_mutation = paste(is.tail)

  mutations = full_join(all_clusters, mutations, by = 'id')
  # mutations = mutations %>% arrange(desc(k))
  
  cl_cols = colnames(mutations)[startsWith(colnames(mutations), 'MOBSTER')]
  vf_cols = colnames(mutations)[endsWith(colnames(mutations), '.VAF')]
  
  # View(mutations %>% mutate(medianVAF = median(!!vf_cols)))
  
  # mutations = mutations %>% arrange(
  #   desc(k),
  #   `Binomial`,
  #   `MOBSTER Set7_55`,
  #   `MOBSTER Set7_57`,
  #   `MOBSTER Set7_59`,
  #   `MOBSTER Set7_62`
  # )
  
  mutations = mutations %>% arrange_(.dots = c('desc(Above_cutoff)', 'Binomial', cl_cols))
  
  mutations = mutations %>%
    group_by_(.dots =  c('Above_cutoff', 'Binomial', cl_cols)) %>%
    arrange(VAF_median, .by_group = TRUE) %>%
    ungroup()
  
  
  # mutations = mutations[
  #   with(mutations, order(cl_cols)),
  #   ]
  
  require(pheatmap)

  core_pheat = as.data.frame(mutations %>% select(ends_with('.VAF')))
  rownames(core_pheat) = mutations$id
  
  
  
  core_pheat_colors = rev(mobster:::scols(1:50, 'Spectral'))
  
  core_pheat_colors = viridisLite::viridis(50, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
  core_pheat_colors[1] = 'gray95'
  
  # viridisLite::viridis(50, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
  M = max(mutations %>% select(ends_with('.VAF')))
  bin_color = M/50
  core_pheat_colors[floor(CCF.cutoff/bin_color)] = 'red'
  
  M = max(mutations %>% select(ends_with('.VAF')), na.rm = TRUE)
  bin_color = M/50
  above = floor(abs((M-CCF.cutoff)/0.02))
  below = floor(abs(CCF.cutoff/0.02))

  core_pheat_colors = 
    c(
      'gray95',
      (mobster:::scols(1:below, 'YlGnBu')),
      'red',
      mobster:::scols(1:above, 'Oranges')
    )
  
  col_below = mobster:::scols(1:below, 'YlGnBu')
  core_pheat_colors = 
    c('gray95', col_below, rep(col_below[length(col_below)], above))
  core_pheat_colors[floor(CCF.cutoff/(M/35))] = 'red'
  
  
  
  
  # Annotate binarized VAF
  annotations_binarized_vaf = as.data.frame(mutations %>% select(ends_with(suffix), Above_cutoff))
  rownames(annotations_binarized_vaf) = mutations$id
  
  # annotations_binarized_vaf[which(annotations_binarized_vaf == 1, arr.ind = T)] = "TRUE"
  # annotations_binarized_vaf[which(annotations_binarized_vaf != "TRUE", arr.ind = T)] = "FALSE"
  
  # Annotate clusters
  annotations_clusters = as.data.frame(mutations %>% select(Binomial, starts_with('MOBSTER'), tail_mutation))
  rownames(annotations_clusters) = mutations$id
  
  # Color Binomial clusters
  colors_B_clusters = mobster:::scols(unique(annotations_clusters$Binomial), 'Set1')
  
  colors_B_clusters[is.na(names(colors_B_clusters))] = 'gray'
  
  annotation_colors = list(Binomial = colors_B_clusters)
  
  
  
  # Color MOBSTER clusters
  for(sample in x$samples) {
    # assgn = paste0('cluster.', sample, '.MOBSTER')
    assgn = paste0('MOBSTER_', sample)
    
    lbl = unique(annotations_clusters[, assgn])
    lbl = sort(lbl, na.last = T)
    
    pal = rev(wesanderson::wes_palette("FantasticFox1"))
    colors_M_clusters = mobster:::scols(lbl, 'Set2')
    
    colors_M_clusters = pal[1:length(colors_M_clusters)]
    names(colors_M_clusters) = lbl
    
    colors_M_clusters[is.na(names(colors_M_clusters))] = 'gray'
    colors_M_clusters[names(colors_M_clusters) == 'Tail'] = 'darkgray'
    
    l = list(colors_M_clusters)
    names(l) = assgn
    
    annotation_colors = append(annotation_colors, l)
  }
  
  tail_yesno_color = list(c(`TRUE`='darkgray', `FALSE` = 'white'))
  names(tail_yesno_color) = 'tail_mutation'
  
  annotation_colors = append(annotation_colors, tail_yesno_color)
  
  annotations_binarized_vaf = annotations_binarized_vaf[mutations$id, , drop = F]
  annotations_clusters = annotations_clusters[mutations$id, , drop = F]
  
  
  gaps_Above_cutoff = unique(annotations_binarized_vaf$Above_cutoff)
  gaps_Above_cutoff = lapply(gaps_Above_cutoff, function(s) which(annotations_binarized_vaf$Above_cutoff == s))
  gaps_Above_cutoff = sapply(gaps_Above_cutoff, function(w)w[1])
  
  gaps_Above_cutoff = gaps_Above_cutoff[!is.na(gaps_Above_cutoff)]
  
  p = pheatmap(
    main = title,
    core_pheat, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    show_rownames = FALSE,
    color = core_pheat_colors, 
    gaps_row = gaps_Above_cutoff,
    # annotation_legend = FALSE,
    annotation_row = cbind(annotations_binarized_vaf, annotations_clusters),
    annotation_colors = annotation_colors,
    cellwidth = 20,
    border_color = NA,
    silent = TRUE
  )
  
  # library(grid)
  
  # grid.gedit("GRID.rect.531", gp = gpar(cex=.1))
  # if(cex_legend != 1)
  # grid::grid.gedit("annotation_legend", gp = gpar(cex=cex_legend))
  p$gtable$grobs[[6]] = grid::editGrob(p$gtable$grobs[[6]], gp = gpar(cex=cex_legend))
  
  p
}



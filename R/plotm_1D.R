
plotm_1D = function(x, 
                    sample_name,
                    cex = 1){
  
  NV_df <- as.data.frame(x$NV)
  DP_df <- as.data.frame(x$DP)
  
  vaf = x$NV/x$DP
  vaf_df <- as.data.frame(vaf)
  
  # rename columns using sample_names
  colnames(vaf_df) <- paste0("vaf_", x$sample_names)
  
  # Create the dataframe
  df <- cbind(mutation_id = x$mutation_id,
              vaf_df,
              cluster_id = x$cluster_id)
  
  df = df %>% mutate(cluster_id=paste0("C",cluster_id+1))
  
  data = mobster:::getm_1D_points(df, sample_name)
  # data = getm_1D_points(df, sample_name)
  
  cluster = df$cluster_id
  
  color_palette = c(
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#a65628",
    "#FFD700",  "#999999", "#000000", "#f781bf", # First 10 colors (Set1)  
    "#46f0f0", "#f032e6", "#bcf60c", "#fabed4", "#008080", "#e6beff",  
    "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",  
    "#000075", "#808080", "#d3a6f3", "#ff9cdd", "#73d7b0"  
  ) %>% setNames(str_sort(unique(cluster), numeric=T))
  
  data$cluster = cluster
  
  # Now the actual plot
  cluster_order = data %>%
    count(cluster) %>%       # Count the number of data points per cluster
    arrange(desc(n)) %>%           # Sort clusters by size (descending order)
    pull(cluster) 
  
  # data <- data %>%
  #   mutate(cluster = factor(cluster, levels = cluster_order))
  data = data %>% filter(.data[[sample_name]] > 0)

  data = data %>%
    filter(get(sample_name) > 0)
  
  plot = ggplot() +
    geom_histogram(data=data, aes(x=eval(parse(text = sample_name)), fill=cluster), 
                   position="identity", 
                   alpha=1, 
                   bins=100) +
    labs(
      title = bquote(bold(.(sample_name))),
      x = 'VAF',
      y = 'Count'
    ) +
    scale_color_manual(values=color_palette, name="Cluster",
                       breaks=str_sort(names(color_palette), numeric=TRUE)) +
    scale_fill_manual(values=color_palette, name="Cluster",
                      breaks=str_sort(names(color_palette), numeric=TRUE)) +
    guides(color = guide_legend(title = 'Cluster', override.aes = list(alpha = 1))) +
    my_ggplot_theme()
  
  plot
}

getm_1D_points = function(df, s)
{
  cols <- paste0("vaf_", s)
  
  data = df[ ,cols, drop=FALSE] %>% as_tibble()
  
  colnames(data) = s
  
  data
}

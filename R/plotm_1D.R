plotm_1D = function(x,
                    color_palette,
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
  
  data <- df %>%
    pivot_longer(
      cols = starts_with("vaf_"),        # columns to pivot
      names_to = "sample",               # new column for sample names
      names_prefix = "vaf_",             # remove this prefix
      values_to = "VAF"                  # new column for values
    )
  data = data %>% mutate(cluster=cluster_id)
  
  cluster = data$cluster_id
  
  color_palette = color_palette %>% setNames(str_sort(unique(cluster), numeric=T))
  
  # data$cluster = cluster
  
  # Now the actual plot
  cluster_order = data %>%
    count(cluster) %>%       # Count the number of data points per cluster
    arrange(desc(n)) %>%           # Sort clusters by size (descending order)
    pull(cluster) 
  
  # data <- data %>%
  #   mutate(cluster = factor(cluster, levels = cluster_order))
  # data = data %>% filter(.data[[x$sample_names]] > 0)
  # 
  # data = data %>%
  #   filter(get(x$sample_names) > 0)
  data <- data %>%
    filter(VAF>0)
  
  plot = data %>% ggplot() +
    geom_histogram(aes(x=VAF, fill=cluster), position="identity", alpha=1, bins=100) +
    facet_grid(~sample) +
    labs(
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

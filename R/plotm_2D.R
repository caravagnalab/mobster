
plotm_2D = function(x, 
                    d1, 
                    d2,
                    color_palette,
                    cex = 1){
  
  NV_df <- as.data.frame(x$NV)
  DP_df <- as.data.frame(x$DP)
  
  vaf = x$NV/x$DP
  vaf_df <- as.data.frame(vaf)
  
  # rename columns using sample_names
  colnames(NV_df) <- paste0("NV_", x$sample_names)
  colnames(DP_df) <- paste0("DP_", x$sample_names)
  colnames(vaf_df) <- paste0("vaf_", x$sample_names)
  
  # Create the dataframe
  df <- cbind(mutation_id = x$mutation_id,
              NV_df, 
              DP_df, 
              vaf_df,
              cluster_id = x$cluster_id)
  
  df = df %>% mutate(cluster_id=paste0("C",cluster_id+1))
  
  data = mobster:::getm_2D_points(df, d1, d2)
  
  cluster = df$cluster_id
  
  color_palette = color_palette %>% setNames(str_sort(unique(cluster), numeric=T))
  
  
  # Now the actual plot
  plot = ggplot() + 
    geom_point(data=data,
             aes(
               x = eval(parse(text = d1)),
               y = eval(parse(text = d2)),
               colour = cluster,
               fill=cluster
             ), size=1 * cex, alpha=1) +
    labs(
      title = bquote(bold(.(d1)) ~ "vs" ~ bold(.(d2))),
      x = d1,
      y = d2
    ) +
    scale_color_manual(values=color_palette, name="Cluster",
                       breaks=str_sort(names(color_palette), numeric=TRUE)) +
    scale_fill_manual(values=color_palette, name="Cluster",
                      breaks=str_sort(names(color_palette), numeric=TRUE)) +
    guides(color = guide_legend(title = 'Cluster', override.aes = list(alpha = 1))) +
    my_ggplot_theme() +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      legend.text = element_text(size = 8 * cex)
    ) +
    geom_vline(xintercept = 0,
               colour = "darkgray") +
    geom_hline(yintercept = 0,
               colour = "darkgray") +
    guides(fill = 'none')
    
    plot
}

getm_2D_points = function(df, d1, d2)
{
  cols <- paste0("vaf_", c(d1, d2))
  
  data = df[ ,cols, drop=FALSE] %>% as_tibble()
  
  colnames(data) = c(d1, d2)
  
  data
}

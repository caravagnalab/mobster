# # Driver annotations 
# require(CNAqc)
# obj = init(CNAqc::example_dataset_CNAqc$snvs, CNAqc::example_dataset_CNAqc$cna, CNAqc::example_dataset_CNAqc$purity)
# obj = subset_by_segment_totalcn(obj, 2)
# inp = obj$snvs
# 
# x = mobster::mobster_fit(inp, auto_setup = 'FAST')$best
# plot(x)

add_extra_plot_annotations = function(x, annotation_extras, base_plot, binwidth = 0.01)
{
  all_list = NULL
  
  # Check for drivers already available inside x
  if(all(c('is_driver', 'driver_label') %in% colnames(x$data)))
    all_list = all_list %>%
      dplyr::bind_rows(x$data %>% dplyr::filter(is_driver))

  # Extra(s)
  if(!is.null(annotation_extras))
    all_list = all_list %>% dplyr::bind_rows(annotation_extras)
  
  # If there is anything to add
  if(!is.null(all_list) & nrow(all_list) > 0)
  {
    all_list = all_list %>% dplyr::select(VAF, driver_label)
    
    # Get the density of the driver's VAF
    points_df = mobster::Clusters_denovo(x, all_list)
    points_df$density = sapply(
      points_df$VAF, 
      function(v){
        
        mobster:::template_density(x,
                                   x.axis = v,
                                   binwidth = binwidth,
                                   reduce = TRUE) %>%
          dplyr::summarise(d = sum(y)) %>%
          pull(d)
      })
    
    nudge = max(points_df$density, na.rm = TRUE)/5
    
    base_plot = base_plot + 
      geom_point(
        data = points_df, 
        aes(x = VAF, y = density), 
        show.legend = F
      ) +
      ggrepel::geom_label_repel(
        data = points_df, 
        aes(x = VAF, y = density, label = driver_label),
        nudge_y = nudge,
        nudge_x = nudge,
        direction = 'x', 
        angle = 45,
        vjust = 0,
        segment.size = 0.2,   
        show.legend = F,
        size = 2
      ) 
  }
  
  return(base_plot)
}
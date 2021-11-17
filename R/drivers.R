has_drivers_annotated = function(x)
{
  return(all(c('is_driver', 'driver_label') %in% colnames(x$data)))
}

add_extra_plot_annotations = function(x,
                                      annotation_extras,
                                      base_plot,
                                      binwidth = 0.01)
{
  has_drivers = mobster:::has_drivers_annotated(x)
  has_annotations = !all(is.null(annotation_extras))
  
  all_list = NULL
  
  if (has_drivers & has_annotations)
    all_list = dplyr::bind_rows(x$data %>% dplyr::filter(is_driver), annotation_extras)
  
  if (has_drivers & !has_annotations)
    all_list = x$data %>% dplyr::filter(is_driver)
  
  if (!has_drivers & has_annotations)
    all_list = annotation_extras
  
  # If there is anything to add
  if (!is.null(all_list))
  {
    if (nrow(all_list) == 0)
      return(base_plot)
    
    all_list = all_list %>% dplyr::select(VAF, driver_label)
    
    # Get the density of the driver's VAF
    points_df = mobster::Clusters_denovo(x, all_list)
    points_df$density = sapply(points_df$VAF,
                               function(v) {
                                 mobster:::template_density(x,
                                                            x.axis = v,
                                                            binwidth = binwidth,
                                                            reduce = TRUE) %>%
                                   dplyr::summarise(d = sum(y)) %>%
                                   pull(d)
                               })
    
    nudge = max(points_df$density, na.rm = TRUE) / 5
    
    base_plot = base_plot +
      geom_point(data = points_df,
                 aes(x = VAF, y = density),
                 show.legend = F) +
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
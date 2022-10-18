plot_random_samples_overlay <- function(obj, predictive_samples, type = c("density", "ecdf")){
  
  used_mutations <- obj$used_mutations
  
  data_real <- obj$data %>% mutate(mutation_id = paste(chr,from,to, sep = ":")) %>% filter(mutation_id %in% used_mutations) %>% select(NV,DP, karyotype)
  data_real <- split(data_real, data_real$karyotype, drop = T)
  data_real <- lapply(data_real, function(x) x[,1])
  plot_list = vector(length = length(names(obj$model_parameters)) , mode = "list")
  names(plot_list) <- names(obj$model_parameters)
  
  for(k in names(obj$model_parameters)) {
    
    if(type == "density"){
      
      density_samples <- lapply(predictive_samples, function(x) density(x[[k]]))
      
      x_samples <- lapply(density_samples, function(x) x$x) %>% do.call(c,.)
      y_samples <- lapply(density_samples, function(x) x$y) %>% do.call(c,.)
      
      
      density_real <- density(data_real[[k]])
      
      x_real <- density_real$x
      y_real <- density_real$y
      
    
    } else {
      
      x <- seq.default(0,max(data_real[[k]]) + 10,by = 1)
      
      
      density_samples <- lapply(predictive_samples, function(x) ecdf(x[[k]]))
      
      x_samples <- lapply(density_samples, function(y) return(x)) %>% do.call(c,.)
      y_samples <- lapply(density_samples, function(y) y(x)) %>% do.call(c,.)
      
      density_real <- ecdf(data_real[[k]])
      x_real <- x
      y_real <- density_real(x)
      
    }
    df_samples <- data.frame(x = x_samples, y = y_samples, alpha = 0.3,  color = "aquamarine1")
    df_real <- data.frame(x = x_real, y = y_real, alpha = 1, color = "black")
    df_plot <- rbind(df_samples, df_real)
    
    plot_list[[k]] <- ggplot(df_plot, aes(x = x, y = y, alpha = paste(alpha), color = color, size = color)) + geom_line() + theme_bw() + 
      ggtitle(paste0(type, " plot for karyotype ", k)) + scale_alpha_manual("", labels = c("simulated data", "real data"), values = c(0.25, 1), guide = "none") + 
      scale_color_manual("", labels = c("simulated data", "real data"), values = c("darkred", "black")) + theme() + xlab("") + ylab("")+ 
      scale_size_manual("", labels = c("simulated data", "real data"), values = c(0.5, 1.5), guide = "none")
    
  }
  
  cowplot::plot_grid(
    plotlist = plot_list,
    nrow = 1,
    align = 'h',
    axis = 'tb'
  )
  
}


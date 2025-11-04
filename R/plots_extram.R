#' Plot the delta parameter of a MOBSTERm fit.
#' 
#' @description Returns a plot of the delta parameter.
#'
#' @param x An object of class \code{"dbpmm_m"}.
#'
#' @return A ggplot object for the plot.
#' 
#' @export
#'
#' @examples TBD
#' 
plot_deltas = function(x)
{
  deltas = x$model_parameters$delta_param # dim(fit$best_fit$model_parameters$delta_param): K x D x 3
  sample_names = x$sample_names
  distribution_names = c("ParetoBinomial", "BetaBinomial", "Dirac")
  # cluster = x$cluster_id
  
  deltas_df <- as.data.frame.table(deltas, responseName = "value") %>%
    rename(cluster = Var1, sample = Var2, distribution = Var3) %>%
    mutate(
      cluster = paste0("C",as.integer(cluster)),
      sample = sample_names[as.integer(sample)],
      distribution = distribution_names[as.integer(distribution)]
    )
  
  pl_deltas = deltas_df %>% 
    group_by(cluster) %>% 
    mutate(distribution=ifelse(distribution=="BetaBinomial", "Beta-\nBinomial",
                               ifelse(distribution=="ParetoBinomial", "Pareto-\nBinomial", distribution))) %>% 
    ggplot() +
    geom_raster(aes(y=factor(distribution, levels=c("Dirac","Beta-\nBinomial","Pareto-\nBinomial")), 
                    x=sample, fill=value)) +
    facet_grid(~ factor(cluster, levels=str_sort(unique(cluster), numeric=TRUE))) +
    scale_fill_distiller(palette="Oranges", limits=c(0,1), breaks=c(0,0.5,1), direction=1,
                         name=expression(delta*" value")) +
    xlab("Sample") +
    my_ggplot_theme() + theme(axis.title.y=element_blank(), 
                     axis.text.x=element_text(size=8, angle=90))
  
  pl_deltas
}


#' Plot the mixing proportions a MOBSTERm fit.
#' 
#' @description Returns a plot of the mixing proportions.
#'
#' @param x An object of class \code{"dbpmm_m"}.
#'
#' @return A ggplot object for the plot.
#' 
#' @export
#'
#' @examples TBD
#' 
plot_mixing_proportions = function(x, color_palette=NA)
{
  if (is.null(color_palette) || all(is.na(unlist(color_palette)))) {
    color_palette = c(
      "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00","#a65628",
      "#FFD700",  "#999999", "#000000", "#f781bf", # First 10 colors (Set1)  
      "#46f0f0", "#f032e6", "#bcf60c", "#fabed4", "#008080", "#e6beff",  
      "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1",  
      "#000075", "#808080", "#d3a6f3", "#ff9cdd", "#73d7b0"  
    ) 
  }
  weights = x$model_parameters$weights_param # dim(fit$best_fit$model_parameters$delta_param): K x D x 3
  sample_names = x$sample_names
  
  c = unique(x$cluster_id) %>% order()
    
  weights_df <- data.frame(
    cluster = paste0("C", c),
    value = weights
  )
  
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
  
  color_palette = color_palette %>% setNames(str_sort(unique(data$cluster_id), numeric=T))
  
  cluster_order = data %>%
    count(cluster) %>%
    arrange(desc(n)) %>%
    pull(cluster) 
  
  weights_df <- weights_df %>%
    mutate(cluster = factor(cluster, levels = str_sort(unique(cluster), numeric = TRUE)))
  
  cluster_counts <- data %>%
    count(cluster) %>%
    mutate(
      cluster = factor(cluster, levels = str_sort(unique(cluster), numeric = TRUE)),
      proportion = n / sum(n)
    )
  
  total_n <- sum(cluster_counts$n)
  
  pl_counts_prop = cluster_counts %>%
    ggplot(aes(x = cluster, y = n, fill = cluster)) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = color_palette,
      labels = paste0(levels(cluster_counts$cluster), " (", cluster_counts$n, ", ", sprintf("%.3f", cluster_counts$proportion), ")"),
      name = "Cluster"
    ) +
    scale_y_continuous(
      name = "Number of mutations",
      sec.axis = sec_axis(~ . / total_n, name = "Mixture weights")
    ) +
    labs(x = "Cluster") +
    my_ggplot_theme() +
    guides(fill = guide_legend(nrow = 2))
  
  pl_counts_prop
  
}
  
#' Plot the responsibilities a MOBSTERm fit.
#' 
#' @description Returns a plot of the responsibilities.
#'
#' @param x An object of class \code{"dbpmm_m"}.
#'
#' @return A ggplot object for the plot.
#' 
#' @export
#'
#' @examples TBD
#' 
plot_responsibilities = function(x)
{
  responsib = x$model_parameters$responsib # dim(x$model_parameters$responsib): K x N
  sample_names = x$sample_names
  
  respons_df = as.data.frame.table(responsib, responseName = "value")%>%
    rename(cluster = Var1, mutation = Var2, responsibility = value) %>% 
    mutate(
      cluster = paste0("C",as.integer(cluster)),
      mutation = as.integer(mutation)
    )
  
  pl_responsib = respons_df %>%
    ggplot() +
    geom_raster(aes(x = mutation, 
                    y = cluster, fill = responsibility)) +
    # facet_grid(~ cluster) +
    labs(x = "Mutations", y = "Cluster") +
    scale_fill_distiller(palette = "Oranges", limits = c(0, 1), breaks = c(0, 0.5, 1), direction = 1, name = "Responsibility\nvalue") +
    scale_x_continuous(breaks = round(seq(min(responsibilities$mutation), max(responsibilities$mutation), length.out = 8))) +
    my_ggplot_theme()
  
  pl_responsib
}

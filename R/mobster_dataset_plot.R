# S3 Class objects for 'mbst_data'

#' Plot a dataset.
#' 
#' @details It shows the VAF distributions of the data, with a special panel for
#' low-frequency VAF values, and panels for coverage (DP). Number of reads with the
#' mutant allele (NV) are not shown as they are as the VAF. There is an option to
#' plot the VAF value by kariotype.
#' 
#' @param x 
#' @param title Title of the plot, default is taken from \code{x}. 
#' @param scales Scales for ggplot facets, default is \code{'fixed'} but \code{'free'} can also help.
#'
#' @return A ggpubr figure object.
#' @export
#'
#' @examples
plot.mbst_data = function(x, 
                         # by_genotype = FALSE, 
                         title = x$description,
                         scales = 'fixed')
{
  by_genotype = FALSE
  
  if(!by_genotype)
  { 
    # Full VAF plot
    vaf = VAF(x) %>% filter(value > 0)
    
    max.VAF = max(vaf$value)
    if(max.VAF < 1) max.VAF = 1
    
    pl_vaf = ggplot(vaf, aes(value, fill = sample)) +
      geom_histogram(binwidth = 0.01) +
      # geom_density() +
      facet_wrap(~sample, nrow = 1, scales = scales) +
      theme_light(base_size = 8) +
      guides(fill = FALSE) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")
      ) +
      geom_vline(aes(xintercept = 0.5), colour = 'red', linetype = "longdash", size = .3) +
      geom_vline(aes(xintercept = 0.25), colour = 'red', linetype = "longdash", size = .3) +
      geom_vline(aes(xintercept = 0.33), colour = 'blue', linetype = "longdash", size = .3) +
      geom_vline(aes(xintercept = 0.66), colour = 'blue', linetype = "longdash", size = .3) +
      geom_vline(aes(xintercept = 1), colour = 'black', linetype = "longdash", size = .2) +
      # geom_density() +
      xlim(0, max.VAF) + 
      labs(
        title = "Variant Allele Frequency (VAF)",
        x = 'VAF')
    
    # Low-frequency VAF
    vaf = vaf %>% filter(value < 0.2)
    
    pl_vaf_lfreq = ggplot(vaf, aes(value, fill = sample)) +
      geom_histogram(binwidth = 0.01) +
      # geom_density() +
      facet_wrap(~sample, nrow = 1, scales = scales) +
      theme_light(base_size = 8) +
      guides(fill = FALSE) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")
      ) +
      labs(
        title = "Low-freq. VAF < 0.20",
        x = 'VAF')
    
    DP_val = DP(x) %>% filter(value > 0)
    
    stats_DP_val = quantile(DP_val %>% pull(value))
    stats_DP_val = data.frame(value = stats_DP_val, quantile = names(stats_DP_val))
    
    pl_DP = ggplot(DP_val, aes(value, fill = sample)) +
      geom_histogram(binwidth = 1) +
      geom_vline(xintercept = stats_DP_val$value, size = .3, linetype = 'dashed', show.legend = TRUE) +
      geom_vline(xintercept = median(DP_val$value), size = .5, colour = 'steelblue') +
      facet_wrap(~sample, nrow = 1,  scales = scales) +
      guides(fill = FALSE, color = guide_legend(title="Quantile")) +
      labs(
        title = "DP histogram (raw data)",
        subtitle = paste0('Quantiles (dashed):', 
                          paste0(stats_DP_val$quantile, collapse = ', '),
                          '. Median (blue)  = ', round(median(DP_val$value), 0), 
                          'x (var. ', round(var(DP_val$value), 0), ')'),
        x = 'Coverage (DP)') +
      theme_light(base_size = 8) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")
      ) 
    
    # 
    # v = DP_val %>%
    #   group_by(sample) %>%
    #   summarize(
    #     # n = n(),
    #     min = min(value),
    #     max = max(value),
    #     mean = mean(value),
    #     sd = sd(value),
    #     median = median(value)
    #   ) 
    # 
    # # v.size = v %>% select(-sample)
    # # C.v.size = apply(v.size, 2, max)
    # # v.size = data.frame(t(t(v.size)/C.v.size))
    # # v.size$sample = v$sample
    # # 
    # v = reshape2::melt(v)
    # v$value = round(v$value)
    # 
    # ggballoonplot(v, 
    #               x = "sample", y = "variable",
    #               size = "value", fill = "value", size.range = c(1, 5), ) 
    # scale_fill_gradientn(colors = my_cols) +
    # guides(size = FALSE) +
    # facet_wrap(~variable, scales = 'free')
    
    
    figure = ggpubr::ggarrange(
      ggpubr::ggarrange(pl_vaf, pl_vaf_lfreq, ncol = 1, heights = c(2, 1.5), nrow = 2),
      pl_DP,
      ncol = 1,
      nrow = 2)
    
    figure = ggpubr::annotate_figure(figure, top = title)
    
    return(figure)
  }
  
  
  values = x$VAF_cn_adjustment %>% 
    mutate(Genotype = paste0(Major, ':', minor))
  
  pio::pioStr('Available genotypes',
              paste0(unique(values$Genotype), collapse = ', '))
  
  ggplot(values, aes(value, fill = sample)) +
    geom_histogram(binwidth = 0.01) +
    facet_grid(Genotype ~ sample, scales = 'free') +
    theme_light(base_size = 8) +
    guides(fill = FALSE) +
    theme(legend.position="bottom",
          legend.text = element_text(size = 8),
          legend.key.size = unit(.3, "cm")
    ) +
    geom_vline(aes(xintercept = 0.5), colour = 'red', linetype = "longdash", size = .3) +
    geom_vline(aes(xintercept = 0.25), colour = 'red', linetype = "longdash", size = .3) +
    geom_vline(aes(xintercept = 0.33), colour = 'blue', linetype = "longdash", size = .3) +
    geom_vline(aes(xintercept = 0.66), colour = 'blue', linetype = "longdash", size = .3) +
    geom_vline(aes(xintercept = 1), colour = 'black', linetype = "longdash", size = .2) +
    # xlim(0, max.VAF) + 
    labs(
      title = "VAF histogram (raw data)",
      x = 'VAF')
  
  stop("Not yet implemented -- plot histograms by CN genotype")
  # 
  # segments = x$segments %>% group_by(sample) %>% spread(variable, value) %>% mutate(CN = paste0(Major, ':', minor))
  # 
  # for(s in x$samples)
  # {
  #   sample.segments = segments %>% filter(sample == s)
  #   
  #   
  #   sample.segments %>% group_by(id) %>% mutate(CN = minor + Major)
  #   
  #   
  # }
  
  
  
}
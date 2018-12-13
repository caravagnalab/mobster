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
    # Full VAF plot, and zoom
    vaf = VAF(x) %>% filter(value > 0)
    lf_vaf = vaf %>% filter(value < 0.2)
    
    # Plus coverage
    DP_val = DP(x) %>% filter(value > 0)
    
    stats_DP_val = quantile(DP_val %>% pull(value))
    stats_DP_val = data.frame(value = stats_DP_val, quantile = names(stats_DP_val))
    # 
    
    max.VAF = max(vaf$value)
    if(max.VAF < 1) max.VAF = 1
    
    pl_vaf = ggplot(vaf, aes(value, fill = sample)) +
      facet_wrap(~sample, nrow = 1, scales = 'free') +
      theme_light(base_size = 8) +
      guides(fill = FALSE) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")
      ) +
      scale_fill_brewer(palette = 'Set1') +
      geom_vline(aes(xintercept = 0.5), colour = 'black', size = .3) +
      geom_vline(aes(xintercept = 0.25), colour = 'darkblue', linetype = "longdash", size = .5) +
      geom_vline(aes(xintercept = 0.33), colour = 'steelblue', linetype = "longdash", size = .5) +
      geom_vline(aes(xintercept = 0.66), colour = 'blue', linetype = "longdash", size = .3) +
      geom_vline(aes(xintercept = 1), colour = 'black', linetype = "longdash", size = .2) +
      geom_histogram(binwidth = 0.01, alpha = .8) +
      # geom_density() +
      xlim(0, max.VAF) +
      labs(
        title = "Adjusted VAF",
        x = 'Data', 
        y = 'Counts'
      )
    
    # Low-frequency VAF
    vaf = vaf %>% filter(value < 0.2)
    
    pl_vaf_lfreq = ggplot(vaf, aes(value, fill = sample)) +
      # geom_density() +
      facet_wrap(~sample, nrow = 1, scales = 'free') +
      theme_light(base_size = 8) +
      scale_fill_brewer(palette = 'Set1') +
      guides(fill = FALSE) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")
      ) +
      labs(
        title = "Adjusted low-frequency VAF (<20%)",
        x = 'VAF', 
        y = 'Counts'
      ) +
      geom_vline(aes(xintercept = 0.01), colour = 'gray', size = .2) +
      geom_vline(aes(xintercept = 0.02), colour = 'gray', size = .2) +
      geom_vline(aes(xintercept = 0.03), colour = 'gray', size = .2) +
      geom_vline(aes(xintercept = 0.04), colour = 'gray', size = .2) +
      geom_vline(aes(xintercept = 0.05), colour = 'gray', size = .2) +
      geom_histogram(binwidth = 0.001, alpha = .8)
    
    
    
    pl_DP = ggplot(DP_val, aes(value, fill = sample)) +
      geom_histogram(bins = 100, alpha = .8) +
      geom_vline(xintercept = stats_DP_val$value, size = .3, linetype = 'dashed', show.legend = TRUE) +
      geom_vline(xintercept = median(DP_val$value), size = .5, colour = 'steelblue') +
      facet_wrap(~sample, nrow = 1,  scales = 'free') +
      guides(fill = FALSE, color = FALSE) +
      scale_fill_brewer(palette = 'Set1') +
      labs(
        title = "Coverage distribution",
        subtitle = paste0('Quantiles (dashed):', 
                          paste0(stats_DP_val$quantile, collapse = ', '),
                          '. Median (blue)  = ', round(median(DP_val$value), 0), 
                          'x (var. ', round(var(DP_val$value), 0), ')'),
        y = 'Counts',
        x = 'Coverage (DP)') +
      theme_light(base_size = 8) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")
      ) +
      scale_x_log10()
    
    
    pl_cov_vaf = ggplot(
      bind_rows(DP_val, vaf) %>% spread(variable, value), 
      aes(x = DP, y = VAF, color = sample)) +
      geom_point(alpha = .1) +
      geom_vline(xintercept = stats_DP_val$value, size = .3, linetype = 'dashed', show.legend = TRUE) +
      geom_vline(xintercept = median(DP_val$value), size = .5, colour = 'steelblue') +
      facet_wrap(~sample, nrow = 1,  scales = 'free') +
      guides(fill = FALSE, color = FALSE) +
      scale_color_brewer(palette = 'Set1') +
      labs(title = "Coverage versus VAF",
           subtitle = paste0('Quantiles (dashed):', 
                             paste0(stats_DP_val$quantile, collapse = ', '),
                             '. Median (blue)  = ', round(median(DP_val$value), 0), 
                             'x (var. ', round(var(DP_val$value), 0), ')'),
           x = 'Coverage (DP)') +
      theme_light(base_size = 8) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")) +
      scale_x_log10()
    
    
    figure = ggpubr::ggarrange(
      pl_vaf, 
      pl_vaf_lfreq, 
      pl_cov_vaf,
      pl_DP,
      ncol = 1,
      nrow = 4)
    
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
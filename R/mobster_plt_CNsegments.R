# load("/Users/gcaravagna/Downloads/rCGH/data/hg19.rda")
# hg19 = as_tibble(hg19)
# hg19 = hg19 %>% mutate(from = cumlen)
# hg19 = hg19 %>% mutate(to = from + length)
# hg19 = hg19 %>% mutate(centromerStart = from + centromerStart)
# hg19 = hg19 %>% mutate(centromerEnd = from + centromerEnd)
# hg19 = hg19 %>% mutate(chr = paste0('chr', chrom))
# hg19 = hg19 %>% select(chr, length, from, to, centromerStart, centromerEnd)
# hg19$chr[hg19$chr == 'chr23'] = 'chrX'
# hg19$chr[hg19$chr == 'chr24'] = 'chrY'
# chr_coordinate_hg19 = hg19
# save(chr_coordinate_hg19, file = 'data/chr_coordinate_hg19.RData')


#' Title
#'
#' @param x 
#' @param samples 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
mobster_plt_CNsegments = function(x, 
                                  samples = x$samples,
                                  ...)
{
  
  pl = lapply(samples, 
              function(w) plt_CNsegments(x, sample = w, ...))
  ggarrange(plotlist = pl, ncol = 1, nrow = length(samples))
}
  

plt_CNsegments = function(x, sample, chromosomes = paste0('chr', c(1:22, 'X', 'Y'))) 
{
  data('chr_coordinates_hg19', package = 'mobster')
  
  chr_coordinate_hg19 = chr_coordinate_hg19 %>% filter(chr %in% chromosomes)

  low = min(chr_coordinate_hg19$from)
  upp = max(chr_coordinate_hg19$to) 
  
  baseplot = ggplot(chr_coordinate_hg19) +
    theme_classic() +
    geom_rect(aes(
      xmin = centromerStart, xmax = centromerEnd, ymin = -Inf, ymax = Inf),
      alpha = .3,
      colour = 'gainsboro') +
    geom_vline(xintercept = chr_coordinate_hg19$from, size = 0.3, colour = 'black') +
    geom_label(data = chr_coordinate_hg19,
               aes(x = chr_coordinate_hg19$from,
                   y = -0.5,
                   label = gsub('chr', '', chr_coordinate_hg19$chr)
                   ),
               hjust = 0, colour = 'white', fill = 'black', size = 3) +
    geom_hline(yintercept = 0, size = 1, colour = 'gainsboro') + 
    geom_hline(yintercept = 1, size = .3, colour = 'black', linetype = 'dashed') + 
    labs(
         x = "Coordinate",
         y = "Major/ minor allele", 
         caption = x$description) +
    ggpubr::rotate_y_text() +
    xlim(low, upp) 
    
  # Segments regions
  segments = x$segments %>% 
    filter(sample == !!sample, chr %in% chromosomes) %>% 
    spread(variable, value) %>% arrange(chr)
  
  # if there are 1-1 segments, shadow them
  one_one = segments %>% filter(Major == 1, minor == 1)
  if(nrow(one_one) > 0)
    baseplot = baseplot +     
      geom_rect(
        data = one_one, 
        aes(xmin = from, xmax = to, ymin = -Inf, ymax = Inf),
        alpha = .2,
        fill = 'forestgreen') 
        
  # Segments
  baseplot = baseplot + 
    geom_segment(data = segments,aes(x = from, xend = to, y = Major, yend = Major), size = 1.5, colour = 'red')+
    geom_segment(data = segments %>% mutate(minor = minor - 0.1),
                 aes(x = from, xend = to, y = minor, yend = minor), size = 1, colour = 'steelblue')
  
  # Scale size
  if(max(segments$Major) < 5)
    baseplot = baseplot + ylim(-0.5, 5)
  
  # Mutation VAF
  mutations = x$map_mut_seg %>% 
    filter(chr %in% chromosomes)
  vaf = VAF(x, ids = mutations %>% pull(id), samples = sample) %>% filter(value > 0) %>% pull(id)
  mutations = mutations %>%
    filter(id %in% vaf)
  
  vaf_mutations = VAF(x, ids = mutations %>% pull(id), samples = sample)
  vaf_mutations = vaf_mutations %>% left_join(mutations, by = 'id')
  
  quant = quantile(vaf_mutations$value, probs = c(.1, .99))
  vaf_mutations = vaf_mutations %>% filter(value > quant[1], value < quant[2])
  
  dpplot = ggplot(vaf_mutations, aes(x = from, y = value, color = value)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, alpha = 1, n = 100) +
    scale_fill_distiller(palette= "Spectral") +
    theme_classic() +
    xlim(low, upp) + 
    # ylim()
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.x=element_blank()
    ) +
    ggpubr::rotate_y_text() +
    labs(y = 'Den.') +
    # scale_y_continuous(breaks = scales::pretty_breaks(n = 1), limits = c(0, 1)) +
    # ylim(0, 1) + 
    guides(fill = FALSE)
  
  pplot = ggplot(vaf_mutations, aes(x = from, y = value, color = value)) +
    geom_point(size = 1e-2, alpha = .3) +
    # stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, alpha = 1) +
    # scale_fill_distiller(palette= "Spectral") +
    geom_hline(yintercept = 0.5, col = 'forestgreen', size = .3, linetype = 'dashed') +
    geom_hline(yintercept = 0.25, col = 'steelblue', size = .3, linetype = 'dashed') +
    geom_hline(yintercept = 0.5+0.25, col = 'darkred', size = .3, linetype = 'dashed') +
    geom_hline(yintercept = 0, col = 'black', size = .3, linetype = 'dashed') +
    geom_hline(yintercept = 1, col = 'black', size = .3, linetype = 'dashed') +
    theme_classic() +
    xlim(low, upp) + 
    # scale_x_continuous(breaks = c(low, upp)) +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.x=element_blank()
    ) +
    ggpubr::rotate_y_text() +
    labs(y = 'VAF') +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 1), limits = c(0, 1)) +
    # ylim(0, 1) + 
    scale_color_gradientn(
      colours = c('steelblue', 'forestgreen', 'darkred'), 
      limits = c(0, 1), 
      breaks = seq(0, 1, by = 0.2)) +
    guides(colour = FALSE)
  
  # Mutation coverage
  mutations = x$map_mut_seg %>% 
    filter(chr %in% chromosomes)
  depth = DP(x, ids = mutations %>% pull(id), samples = sample) %>% filter(value > 0) %>% pull(id)
  mutations = mutations %>%
    filter(id %in% depth)
  
  depth_mutations = DP(x, ids = mutations %>% pull(id), samples = sample)
  depth_mutations = depth_mutations %>% left_join(mutations, by = 'id')
  
  q_cov = quantile(depth_mutations$value, probs = c(.01, .99))
  names(q_cov) = q_cov
  
  pplot_depth = ggplot(depth_mutations, aes(x = from, y = value, color = value)) +
    geom_point(size = 1e-2, alpha = .3) +
    geom_hline(yintercept = median(depth_mutations$value), 
               col = 'black', 
               size = .3, 
               linetype = 'dashed') +
    geom_hline(yintercept = q_cov[1], col = 'black', size = .3, linetype = 'dashed') +
    geom_hline(yintercept = q_cov[2], col = 'black', size = .3, linetype = 'dashed') +
    # geom_hline(yintercept = 0.25, col = 'steelblue', size = .3, linetype = 'dashed') +
    # geom_hline(yintercept = 0.5+0.25, col = 'darkred', size = .3, linetype = 'dashed') +
    theme_classic() +
    xlim(low, upp) + 
    # scale_x_continuous(breaks = c(low, upp)) +
    scale_y_continuous(breaks = q_cov, limits = q_cov) +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.x=element_blank()
    ) +
    ggpubr::rotate_y_text() +
    labs(y = 'DP') +
    # scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_color_gradientn(    
      colours = c('steelblue', 'orange', 'darkred'), 
      breaks = c(q_cov[1], median(depth_mutations$value), q_cov[2]),
      limits = q_cov) +
    guides(colour = FALSE)
  
  # Histogram of mutation counts with 1 megabase bins
  binsize = 1e6
  
  hplot = ggplot(mutations, aes(x = from)) +
    geom_histogram(aes(y = ..count..), binwidth = binsize, fill = 'black') 

  m = max(ggplot_build(hplot)$data[[1]]$count) + 20
  hplot = hplot +
    theme_classic() +
    xlim(low, upp) + 
    # scale_x_continuous(low, upp - upp*binsize) +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ggpubr::rotate_y_text() +
    # scale_y_continuous(
    #   breaks = c(0, m+1, m + 10), labels = c(0, m, '')) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2), limits = c(0, m)) +
    labs(y = 'n') 
    # geom_hline(yintercept = m, linetype = 'dashed', size = 0.3, color = 'gray') +
    # annotate("text", x = 0, y = m*0.9, label = ceiling(m), size = 2)

  figure = ggarrange(hplot,  pplot, pplot_depth, baseplot,
    nrow = 4, ncol = 1,
    heights = c(.2, .2, .2, 1)
  )

  annotate_figure(
    figure, 
    top = text_grob(
      bquote(bold('   ' ~ .(sample)) ~ 
               ' -  Copy Number profile (' * .(nrow(segments)) ~ 'seg., n =' ~ .(nrow(mutations)) * ')' ),
      color = "black", 
      face = "bold", 
      hjust = 0,
      x = 0)
  )
}

# ggsave(filename = 'a.pdf', mobster_plt_CNsegments(x), width = 10, height = 15)  

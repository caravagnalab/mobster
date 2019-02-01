#' Plots mutations that would be fitered with MOBSTER function.
#'
#' @description MOBSTER has functions to subset the mutations
#' that one has annotated in the data. With this function one
#' can generate plots that show which mutations would be filtered
#' for a certain set of parameters, this allows to one to tune
#' the filters available in the too.
#'
#' @param x A MOBSTER dataset.
#' @param VAF_min VAF values below this parameter are considered as 0 (undetected);
#' see function \code{\link{mobster_flt_minvaf}}. This filter applies to each
#' sample independently.
#' @param NV_min NV values below this parameter are considered as 0 (undetected);
#' see function \code{\link{mobster_flt_minnv}}. This filter applies to each
#' sample independently.
#' @param min.DP The depth value for each mutation must be higher than \code{min.DP},
#' whenever it is greater than 0 in a sample; see function \code{\link{mobster_flt_dprange}}.
#' @param max.DP The depth value for each mutation must be lower than \code{max.DP};
#' see function \code{\link{mobster_flt_dprange}}.
#' @param min.VAF The VAF value for each mutation must be higher than \code{min.VAF},
#' whenever it is greater than 0 in a sample; see function \code{\link{mobster_flt_vafrange}}.
#' @param max.VAF The depth value for each mutation must be lower than \code{min.VAF};
#' see function \code{\link{mobster_flt_vafrange}}.
#' @param x.lim X-axis limits for the plot (e.g., \code{c(0,1)}); NA is no limits.
#' @param y.lim Y-axis limits for the plot (e.g., \code{c(0,1)}); NA is no limits.
#'
#' @return A set of figures (per filter, and per sample) with all annotated mutations
#' coloured according to their status (if they would be removed by a filter, or not).
#'
#' @export
#'
#' @examples
mobster_plt_filters = function(x,
                               VAF_min,
                               NV_min,
                               min.DP,
                               max.DP,
                               min.VAF,
                               max.VAF,
                               x.lim = NA,
                               y.lim = NA)
{
  # Pairwise plotting function
  aux_fun_plt_filters = function(obj, x, y)
  {
    pp = VAF_table(obj, samples = c(x, y), suffix = '')
    
    #
    # FILTER #1 - VAF > minimum value.
    #
    pp$filter_minVAF =
      apply(pp %>% select(-id),
            1,
            function(w)
            {
              w1 = (w[1] > 0) & (w[1] < VAF_min)
              w2 = (w[2] > 0) & (w[2] < VAF_min)
              
              if (!w1 & !w2)
                return('No')
              else
                return("Yes")
            })
    
    #
    # FILTER #2 - Depth across min and max in all samples
    #
    ids = DP(obj) %>%  filter(value > 0 &
                                (value < min.DP |
                                   value > max.DP)) %>%
      select(id) %>%
      distinct() %>% pull(id)
    
    pp = pp %>%
      mutate(filter_dprange = ifelse(id %in% !!ids, 'Yes', 'No'))
    
    #
    # FILTER #3 - NV > minimum value.
    #
    ppNV = NV_table(obj, samples = c(x, y), suffix = '')
    
    ppNV$filter_minNV =
      apply(ppNV %>% select(-id),
            1,
            function(w)
            {
              w1 = (w[1] > 0) & (w[1] < NV_min)
              w2 = (w[2] > 0) & (w[2] < NV_min)
              
              if (!w1 & !w2)
                return('No')
              else
                return("Yes")
            })
    
    pp$filter_minNV = "No"
    pp[ppNV %>% filter(filter_minNV == 'Yes') %>% pull(id), 'filter_minNV'] = 'Yes'
    
    
    #
    # FILTER #4 - NV > minimum value.
    #
    ids = VAF(obj) %>%  filter(value > 0 &
                                 (value < min.VAF |
                                    value > max.VAF)) %>%
      select(id) %>%
      distinct() %>% pull(id)
    
    pp = pp %>%
      mutate(filter_VAFrange = ifelse(id %in% !!ids, 'Yes', "No"))
    
    pl = function(attr, t, cex = 1, sub) {
      g = ggplot(pp,
                 aes(
                   x = eval(parse(text = x)),
                   y = eval(parse(text = y)),
                   color = eval(parse(text = attr))
                 )) +
        geom_point(size = .5 * cex) +
        theme_light() +
        labs(
          x = x,
          y = y,
          title = t,
          subtitle = sub
        ) +
        guides(color = guide_legend("Filter")) +
        scale_color_manual(values = c(
          `Yes` = 'red',
          `No` = 'black',
          `NA` = 'gray'
        ))
      
      if (!is.na(x.lim))
        g = g + xlim(x.lim[1], x.lim[2])
      if (!is.na(y.lim))
        g = g + ylim(y.lim[1], y.lim[2])
      
      g
    }
    
    p1 = pl(attr = 'filter_minVAF',
            paste0('VAF > ', VAF_min),
            sub = "Filter type: per sample.")
    p2 = pl(attr = 'filter_dprange',
            paste0('DP range [', min.DP, '; ', max.DP , ']'),
            sub = "Filter type: all samples.")
    p3 = pl(attr = 'filter_minNV',
            paste0('NV > ', NV_min),
            sub = "Filter type: per sample.")
    p4 = pl(attr = 'filter_VAFrange',
            paste0('VAF range [', min.VAF, '; ', max.VAF , ']'),
            sub = "Filter type: all samples")
    
    pp = pp %>% mutate(
      flt = (filter_minVAF == 'Yes') |
        (filter_dprange == 'Yes') |
        (filter_minNV == 'Yes') | (filter_VAFrange == 'Yes'),
      filter = ifelse(flt, "Yes", "No")
    )
    
    pall = pl(attr = 'filter', paste0('Filter (any)'),  sub = "Mutations triggering at least one filter.") + facet_wrap(~
                                                                                                                          filter)
    # Arrange the plots
    ggarrange(
      p1,
      p3,
      p2,
      p4,
      pall,
      nrow = 1,
      ncol = 5,
      common.legend = TRUE,
      legend = 'right'
    )
  }
  
  pairs = combn(x$samples, 2)
  
  pl = apply(pairs,
             2,
             function(w) {
               aux_fun_plt_filters(obj = x, x = w[1], y = w[2])
             })
  
  ggarrange(plotlist = pl,
            ncol = 1,
            nrow = ncol(pairs))
}

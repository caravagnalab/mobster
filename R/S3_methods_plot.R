#' Plot a MOBSTER fit.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param alpha Alpha value for the colors of the histogram
#' @param colors If provided, these colours will be used for each cluster.
#' If a subset of colours is provided, palette Set1 from \code{RColorBrewer} is used.
#' By default the tail colour is provided as 'gainsboro'.
#' @param cutoff_assignment Parameters passed to run function \code{Clusters} which
#' returns the hard clustering assignments for the histogram plot if one wants to plot
#' only mutations with responsibility above this parameter.
#' @param annotation_extras A dataframe that contains a label column, and a VAF value.
#' The labels will be annotated to the corresponding clusters of the VAF values.
#' @param secondary_axis \code{NULL} to leave the second axis empty, \code{"SSE"} to 
#' report the cumulative percentage of SSE error, and \code{"N"} to report the
#' cumulative number of mutations.
#' @param ... 
#'
#' @return A ggplot object for the plot.
#' @export
#'
#' @import sads
#' @import ggplot2
#'
#' @examples 
#' data(fit_example)
#' plot(fit_example$best)
#' plot(fit_example$best, secondary_axis = 'SSE')
#' plot(fit_example$best, secondary_axis = 'N')
#' plot(fit_example$best, alpha = 1, cutoff_assignment = .7)
#' plot(fit_example$best, colors =  c(`Tail` = 'gainsboro', `C1` = 'darkorange'))
plot.dbpmm = function(x,
                      alpha = .7,
                      colors = c(`Tail` = 'gainsboro'),
                      cutoff_assignment = 0,
                      annotation_extras = NULL,
                      secondary_axis = NULL,
                      ...
                      )
{
  is_mobster_fit(x)
  
  binwidth = 0.01
  histogram.main = 'MOBSTER fit'
  
  # Prepare variables
  domain = seq(0, 1, binwidth)

  labels = names(mobster:::.params_Pi(x))
  labels.betas = mobster:::.params_Beta(x)$cluster

  pi = mobster:::.params_Pi(x)

  # Main plotting data
  plot_data = mobster:::Clusters(x, cutoff_assignment)
  clusters = sort(unique(plot_data$cluster), na.last = TRUE)
  
  # Text for the plot -- convergence
  conv.steps = length(x$all.NLL)
  conv.epsilon = 0
  if (conv.steps >= 2)
    conv.epsilon = abs(rev(x$all.NLL)[1] - rev(x$all.NLL)[2])
  conv.epsilon = formatC(conv.epsilon, format = "e", digits = 0)

  sse = max(mobster:::.compute_fit_sqerr(x, binning = binwidth)$cum.y)
  sse = formatC(sse, format = "e", digits = 0)

  label.fit = bquote(
    .(x$fit.type) *
      " (" * omega * " = " * .(conv.steps) ~ '; ' * epsilon ~ '=' ~ .(conv.epsilon) *
      "; SSE" ~ .(sse) * '), LV > ' * .(cutoff_assignment)
  )

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main ggplot object is the histogram 
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  hist_pl = ggplot(
    plot_data, 
    aes(VAF, 
        fill = factor(cluster, levels = clusters), 
        y = ..count.. /sum(..count..))) +
    geom_histogram(alpha = alpha,
                   color = NA,
                   position = 'identity',
                   binwidth = binwidth) +
    geom_vline(
      xintercept = min(x$data$VAF),
      colour = 'black',
      linetype = "longdash",
      size = .3
    ) +
    guides(fill = guide_legend(title = "Cluster")) +
    labs(
      title = bquote(.(histogram.main)),
      # subtitle = annotation,
      caption = label.fit,
      x = "Observed Frequency",
      y = "Density"
    ) +
    mobster:::my_ggplot_theme() +
    theme(
      panel.background = element_rect(fill = 'white'),
      plot.caption = element_text(color = ifelse(x$status, "darkgreen",  "red"))
    ) + coord_cartesian(clip = 'off')
  
  # Get the maximum of the histogram
  hist_yMax = max(ggplot_build(hist_pl)$data[[1]]$y)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # We add the density to the histogram
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Template density values
  densities = mobster:::template_density(x,
                               x.axis = domain,
                               binwidth = binwidth,
                               reduce = TRUE)

  # Add the trace and the mean of each component
  hist_pl = hist_pl +
    geom_line(data = densities,
              aes(
                y = y,
                x = x,
                color = factor(cluster, levels = clusters)
              ),
              size = 1) +
    # scale_color_manual(values = colors, labels = names(colors)) +
    guides(color = FALSE)
  
  # The new max, if the density is too high we set the max to the hist max
  hist_den_yMax = max(ggplot_build(hist_pl)$data[[1]]$y)
  
  if(hist_yMax < hist_den_yMax ) hist_pl = hist_pl + ylim(0, hist_yMax)

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Beta peaks (means) annotated to the plot
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  Beta_peaks = x$Clusters %>%
    dplyr::filter(type == 'Mean', cluster != 'Tail')

  hist_pl = hist_pl +
    geom_vline(data = Beta_peaks,
               aes(xintercept = fit.value, 
                   color = factor(cluster, levels = clusters)),
               linetype = "longdash")

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Overall mixture density as \sum_i pi_i f_i
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  densities = tibble::as_tibble(densities)
  densities = densities %>% group_by(x) %>% summarise(y = sum(y), cluster = 'f(x)')

  m = max(densities$y, na.rm = TRUE)

  hist_pl = hist_pl +
    geom_line(
      data = densities %>% mutate(y = y + m * 0.02),
      aes(y = y, x = x),
      color = 'black',
      alpha = .8,
      size = .5,
      linetype = 'dashed',
      inherit.aes = FALSE
    )

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Dashed SSE behing the overall density
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  error = suppressWarnings(mobster:::.compute_fit_sqerr(x, binning = binwidth))

  # scale to percentage for plotting
  if(!is.null(secondary_axis))
  {
    if(secondary_axis == "SSE")
    {
      me = max(error$cum.y)
      error = error %>% dplyr::mutate(cum.y = (cum.y / me) * m)
      
      hist_pl = hist_pl +
        geom_line(
          data = error,
          aes(y = cum.y, x = x),
          color = 'darkgray',
          alpha = 1,
          size = .2,
          linetype = 'dashed',
          inherit.aes = FALSE
        ) +
        scale_y_continuous(sec.axis = sec_axis( ~ . / m * 100, name = "SSE [cumulative %]"))
    }
    
    if(secondary_axis == "N")
    {
      xbins = hist(x$data$VAF, breaks = domain, plot = F)$counts
      xcumc = cumsum(xbins)
      
      brk = seq(0, x$N, by = x$N/4)

      # hist_pl = 
        hist_pl +
        geom_line(
          data =  data.frame(x = domain[-1], y = (xcumc/x$N) * m),
          aes(y = y, x = x),
          color = 'darkgray',
          alpha = 1,
          size = .2,
          linetype = 'dashed',
          inherit.aes = FALSE
        ) +
        scale_y_continuous(
          sec.axis = sec_axis( ~ . , 
                               name = "Counts (n)", 
                               breaks = brk / x$N * m,
                               labels = brk
                               )
          )
    }
    
  }
  

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Annotate mixting proportions
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pi = sort(pi)
  pi = paste0(names(pi), ' ', round(pi * 100, 2), '%', collapse = ', ')
  
  hist_pl = hist_pl +
    labs(
      subtitle = paste0("N = ", x$N, "; ", pi)
    )
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Custom coloring
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # TODO
  
  
  # # Annotation of input entries
  # if(!is.null(annotate) & all(c("VAF", 'label') %in% colnames(annotate)))
  # {
  #   # Position the point at coord x = VAF and y = density
  #   m = max(densities$y, na.rm = TRUE)
  #
  #   annotate$y = round(annotate$VAF/binwidth)
  #   annotate$y = densities$y[annotate$y] + m * 0.02
  #
  #   hist_pl = hist_pl +
  #     geom_label_repel(data = annotate,
  #                      aes(
  #                        x = VAF,
  #                        y = y,
  #                        label = label,
  #                        color = factor(cluster)
  #                        # fill = factor(cluster)
  #                      ),
  #                      size = 1.5 * cex,
  #                      inherit.aes = FALSE,
  #                      box.padding = 0.95,
  #                      segment.size = .2 * cex, force = 1) +
  #     geom_point(data = annotate, aes(x = VAF, y = y,  color = factor(cluster)),
  #                size = .3 * cex, alpha = 1,
  #                inherit.aes = FALSE)
  # }
  #
  
  pl = add_fill_color_pl(x, hist_pl, colors)
  
  if(!is.null(annotation_extras))
  {
    points_df = mobster::Clusters_denovo(x, annotation_extras)
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
    
    nudge = max(points_df$density)/5
    
    pl = pl + 
      geom_point(
        data = points_df, 
        aes(x = VAF, y = density), 
        show.legend = F
      ) +
      ggrepel::geom_label_repel(
        data = points_df, 
        aes(x = VAF, y = density, label = label),
        nudge_y = nudge,
        nudge_x = nudge,
        direction = 'x', 
        angle = 45,
        vjust = 0,
        segment.size = 0.2,   
        show.legend = F
      ) 
  }

  return(pl)
}


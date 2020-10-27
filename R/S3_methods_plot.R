#' Plot a MOBSTER fit.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param beta_colors A vector of colors that are used to colour the Beta clusters.
#' Colors are used by order for \code{"C1"}, \code{"C2"}, \code{"C3"}, etc. If too
#' few colours are used the \code{rainbow} palette is used instead of \code{beta_colors}.
#' In all cases the colour of the tail is set by the other parameter \code{tail_color}.
#' @param tail_color The colour of the tail cluster, if any.
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
#'
#' @import sads
#' @import ggplot2
#'
#' @exportS3Method plot dbpmm
#' @export plot.dbpmm
#'
#' @examples
#' data(fit_example)
#' plot(fit_example$best)
#' plot(fit_example$best, secondary_axis = 'SSE')
#' plot(fit_example$best, cutoff_assignment = .7)
#' plot(fit_example$best, beta_colors = c("indianred3", "orange"))
#' plot(fit_example$best, tail_color = 'black')
plot.dbpmm = function(x,
                      cutoff_assignment = 0,
                      beta_colors = RColorBrewer::brewer.pal(n = 9, 'Set1'),
                      tail_color = "gainsboro",
                      na_color = 'gray',
                      annotation_extras = NULL,
                      secondary_axis = NULL,
                      ...)
{
  mobster:::is_mobster_fit(x)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting variables
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  binwidth = 0.01
  domain = seq(0, 1, binwidth)
  
  labels = names(mobster:::.params_Pi(x))
  labels.betas = mobster:::.params_Beta(x)$cluster
  
  pi = mobster:::.params_Pi(x)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting data
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Clustering assignments and clusters
  plot_data = mobster:::Clusters(x, cutoff_assignment) %>% dplyr::arrange(cluster)
  
  clusters = sort(unique(plot_data$cluster), na.last = TRUE)
  
  # Beta peaks
  Beta_peaks = x$Clusters %>%
    dplyr::filter(type == 'Mean', cluster != 'Tail')
  
  # Component mixture density and total as \sum_i pi_i f_i
  densities = suppressWarnings(mobster:::template_density(
    x,
    x.axis = domain,
    binwidth = binwidth,
    reduce = TRUE
  ))
  
  total_densities = tibble::as_tibble(densities) %>%
    dplyr::group_by(x) %>%
    dplyr::summarise(y = sum(y), cluster = 'f(x)')
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plotting aestetics
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  # Colours per cluster - tail is reserved
  Beta_colors = pio:::nmfy(paste0('C', 1:length(beta_colors)), beta_colors)
  MOBSTER_CLUSTER_COLORS = c(Beta_colors, `Tail` = tail_color)
  
  if (all(clusters %in% names(MOBSTER_CLUSTER_COLORS)))
  {
    MOBSTER_CLUSTER_COLORS = MOBSTER_CLUSTER_COLORS[names(MOBSTER_CLUSTER_COLORS) %in% clusters]
  }
  else
  {
    available =  clusters[clusters %in% names(MOBSTER_CLUSTER_COLORS)]
    missing =  clusters[!(clusters %in% names(MOBSTER_CLUSTER_COLORS))]
    
    warning(
      "You did not pass enough input colours, adding a gray colour\n",
      'Available: ', paste0(available, collapse = ', '), '\n',
      'Missing: ', paste0(missing, collapse = ', ')
    )
    
    MOBSTER_CLUSTER_COLORS = c(MOBSTER_CLUSTER_COLORS,
                               pio:::nmfy(paste(missing), rep("gray", length(missing)))
                               )
    # MOBSTER_CLUSTER_COLORS = rainbow(length(clusters))
  }
  
  # Plot title
  plot_title = ifelse(all(is.null(x$description)), 'mobster fit', x$description)
  
  # Plot substitle
  plot_subtitle = NULL
  
  pi = sort(pi)
  pi = paste0(names(pi), ' ', round(pi * 100, 1), '%', collapse = ', ')
  
  plot_subtitle = paste0("n = ", x$N, "; ", pi)
  
  # Plot caption
  plot_caption = NULL
  
  conv.steps = length(x$all.NLL)
  conv.epsilon = 0
  if (conv.steps >= 2)
    conv.epsilon = abs(rev(x$all.NLL)[1] - rev(x$all.NLL)[2])
  conv.epsilon = formatC(conv.epsilon, format = "e", digits = 0)
  
  plot_caption = bquote(
    .(x$fit.type) *
      " (" * .(conv.steps) ~ 'steps; ' * epsilon ~ '=' ~ .(conv.epsilon) *
      '), ' * z["nk"] * ' > ' * .(cutoff_assignment)
  )
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Main ggplot object is the histogram
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  hist_pl = ggplot(plot_data,
                   aes(VAF,
                       fill = factor(cluster),
                       y = ..count.. / sum(..count..))) +
    geom_histogram(
      alpha = .6,
      color = NA,
      position = 'identity',
      binwidth = binwidth
    ) +
    scale_fill_manual(values = MOBSTER_CLUSTER_COLORS) +
    geom_vline(
      xintercept = min(x$data$VAF),
      colour = 'black',
      linetype = "longdash",
      size = .3
    ) +
    guides(fill = guide_legend(title = "Cluster")) +
    geom_line(data = densities,
              aes(
                y = y,
                x = x,
                color = factor(cluster)
              ),
              size = 1) +
    geom_line(
      data = total_densities %>% mutate(y = y + max(total_densities$y, na.rm = TRUE) * 0.02),
      aes(y = y, x = x),
      color = 'black',
      alpha = .8,
      size = .5,
      linetype = 'dashed',
      inherit.aes = FALSE
    ) +
    geom_vline(data = Beta_peaks,
               aes(xintercept = fit.value,
                   color = factor(cluster)),
               linetype = "longdash") +
    scale_color_manual(values = MOBSTER_CLUSTER_COLORS) +
    guides(color = FALSE) +
    labs(
      title = bquote(.(plot_title)),
      subtitle = plot_subtitle,
      caption = plot_caption,
      x = "Observed Frequency",
      y = "Density"
    ) +
    mobster:::my_ggplot_theme() +
    theme(
      panel.background = element_rect(fill = 'white'),
      plot.caption = element_text(color = ifelse(x$status, "darkgreen",  "red"))
    ) +
    coord_cartesian(clip = 'off')
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Dashed SSE behind the overall density (optional), scaled to percentage
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (!is.null(secondary_axis))
  {
    if (secondary_axis == "SSE")
    {
      error = suppressWarnings(mobster:::.compute_fit_sqerr(x, binning = binwidth))
      
      m = max(ggplot_build(hist_pl)$data[[1]]$y)
      
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
        scale_y_continuous(sec.axis = sec_axis(~ . / m * 100, name = "SSE [cumulative %]"))
    }
    
    # if(secondary_axis == "N")
    # {
    #   xbins = hist(x$data$VAF, breaks = domain, plot = F)$counts
    #   xcumc = cumsum(xbins)
    #
    #   brk = seq(0, x$N, by = x$N/4)
    #
    #   # hist_pl =
    #     hist_pl +
    #     geom_line(
    #       data =  data.frame(x = domain[-1], y = (xcumc/x$N) * m),
    #       aes(y = y, x = x),
    #       color = 'darkgray',
    #       alpha = 1,
    #       size = .2,
    #       linetype = 'dashed',
    #       inherit.aes = FALSE
    #     ) +
    #     scale_y_continuous(
    #       sec.axis = sec_axis( ~ . ,
    #                            name = "Counts (n)",
    #                            breaks = brk / x$N * m,
    #                            labels = brk
    #                            )
    #       )
    # }
    
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Extra annotations (param) plus drivers
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  hist_pl = suppressWarnings(
    mobster:::add_extra_plot_annotations(x, annotation_extras, base_plot = hist_pl)
  )
  
  return(hist_pl)
}

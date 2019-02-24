#' Plot the initial density of a fit.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param cex Cex of the plot
#' @param colors If provided, these colours will be used for each cluster.
#' If not all clusters have available colours, errors are thrown.
#'
#' @return A ggplot object for the plot.
#' @export
#'
#' @examples
#' TODO
plot_init = function(x, cex = 1, colors = NA)
{
  stopifnot(inherits(x, "dbpmm"))
  
  # Simple plot.
  #
  # We just get the desnity with the usual functions
  # after assinging the initial parameters appropriately.
  binwidth = 0.01
  domain = seq(0, 1, binwidth)
  
  n = x
  n$Clusters$fit.value = n$Clusters$init.value
  
  initial.densities = template_density(n,
                                       x.axis = domain[2:(length(domain) - 1)],
                                       # Restricted for numerical errors
                                       binwidth = 0.01,
                                       reduce = TRUE)
  
  den_init_pl = ggplot() +
    labs(title = bquote(bold("Initialization")),
         x = "Adjusted VAF",
         y = "Density") +
    guides(fill = FALSE) +
    ylim(0, max(initial.densities$y)) +
    theme_light(base_size = 8 * cex) +
    geom_line(data = initial.densities, aes(y = y, x = x, color = cluster)) +
    guides(color = FALSE)
  
  
  return(add_color_pl(x, den_init_pl, colors))
}


#' Plot the mixing proportions of the mixture.
#'
#' @param x An object of class 'dbpmm'.
#' @param cex Cex of the plot.
#' @param colors If provided, these colours will be used for each cluster.
#' If not all clusters have available colours, errors are thrown.
#'
#' @return A plot of the mixing proportions of the mixture.
#' @export
#'
#' @examples
#' data(fit_example)
#' plot_mixing_proportions(fit_example$best)
plot_mixing_proportions = function(x, cex = 1, colors = NA)
{
  stopifnot(inherits(x, "dbpmm"))
  
  Proportions = x$Clusters %>%
    dplyr::filter(type == 'Mixing proportion')
  
  Proportions$fit.value = round(Proportions$fit.value, 2)
  
  if (!x$fit.tail)
    Proportions = Proportions %>% filter(cluster != 'Tail')
  
  pl = ggplot(data = Proportions, aes(x = cluster, y = fit.value, fill = cluster)) +
    geom_bar(stat = "identity", width = 0.3 * cex) +
    geom_hline(
      aes(yintercept = 0.02),
      colour = 'red',
      linetype = "longdash",
      size = 0.3
    ) +
    geom_text(
      data = NULL,
      aes(label = '2%', x = 0.1, y = 0.04),
      inherit.aes = FALSE,
      hjust = 0,
      colour = 'red',
      size = 2.5 * cex
    ) +
    labs(title  = bquote(bold('Mixing Proportions'))) +
    xlab("") +
    ylab("") +
    guides(fill = FALSE) +
    theme_light(base_size = 8 * cex) +
    ylim(c(0, 1))
  
  add_fill_pl(x, pl, colors)
}


#' Plot the scores for model selection.
#' 
#' @description Plots the scores via ICL, reICL, BIC and
#' AIC which can be used for model selection. It allows to
#' easily see if the model selected as best is consistently
#' better for all scores.
#'
#' @param x A list of fits computed via \code{mobster_fit}.
#' @param cex Cex of the plot
#'
#' @return A ggplot figure with the scores for model selection.
#' @export
#'
#' @examples
#' TODO
plot_fit_scores = function(x, cex = 1)
{
  stopifnot(is_list_mobster_fits(x))
  
  model.selection = 'ICL'
  if (!is.null(x$model.selection))
    model.selection = x$model.selection
  
  scores = c('ICL', 'reICL', 'BIC', 'AIC')
  
  scores = x$fits.table[, c(scores, 'K', 'tail')]
  scores = scores[complete.cases(scores),]
  
  scores = scores %>% mutate(tail = ifelse(tail, 'With Tail', 'Without Tail'))
  
  ranks = order(scores[, model.selection])
  scores = scores[ranks, , drop = FALSE]
  scores$rank = 1:nrow(scores)
  
  mscores = reshape2::melt(scores, c('K', 'tail', 'rank'))
  
  K_vals = unique(mscores$K)
  
  opt_scores = mscores %>%
    group_by(variable) %>%
    arrange(value) %>%
    filter(row_number() == 1)
  
  ggplot(data = mscores,
         aes(x = rank,
             y = value,
             color = variable)) +
    geom_line() +
    geom_point(aes(
      x = rank,
      y = value,
      shape = tail,
      size = K
    ),
    inherit.aes = F)  +
    geom_point(
      data = opt_scores,
      aes(
        x = rank,
        y = value,
        shape = tail,
        size = K
      ),
      inherit.aes = F,
      show.legend = FALSE,
      color = 'red'
    )  +
    scale_size(range = c(min(K_vals), max(K_vals)) * 0.7 * cex,
               breaks = min(K_vals):max(K_vals)) +
    facet_wrap( ~ variable) +
    theme_light(base_size =  8 * cex) +
    labs(
      title  = bquote(bold('Model selection:') ~ .(model.selection)),
      subtitle = bquote(.(nrow(scores)) ~ ' runs'),
      x = 'Models rank',
      y = 'Score'
    ) +
    guides(
      fill = FALSE,
      shape = guide_legend(title = 'Tail'),
      colour = guide_legend(title = 'Score')
    ) +
    theme(legend.position = 'bottom')
  
}


#' Plot summary for model selection.
#'
#' @description This plot is usefull to understand the different
#' fits and their rank with respect to some scoring used to
#' select the best model. The plot shows alternative solutions
#' and their rank as well.
#'
#' @param x A list of fits computed via \code{mobster_fit}.
#' @param TOP The top set of fits to use.
#' @param cex Cex of the plot.
#'
#' @return A complex figure with all plots arranged.
#'
#' @import ggpubr
#'
#' @export
#'
#' @examples
#' TODO
plot_model_selection = function(x,
                                TOP = 5,
                                cex = 1,
                                ...)
{
  stopifnot(is_list_mobster_fits(x))
  
  # =-=-=-=-=-=-=-=-
  # Best model fit
  # =-=-=-=-=-=-=-=-
  best = plot.dbpmm(x$best,
                    cex = cex)
  
  # =-=-=-=-=-=-=-=-
  # Scores rank best
  # =-=-=-=-=-=-=-=-
  rank = plot_fit_scores(x, cex)
  
  # =-=-=-=-=-=-=-=-
  # Goodness of fit
  # =-=-=-=-=-=-=-=-
  
  # SSE plot
  GOFP = plot_gofit(x, TOP = TOP, cex = cex)
  
  # =-=-=-=-=-=-=-=-=-=-
  # Low-rank solutions
  # =-=-=-=-=-=-=-=-=-=-
  
  # Other solutions ranked below top best -- maximum TOP 4 fixed
  TOP.plot = min(4, nrow(x$fits.table))
  
  other.best = NULL
  
  if (TOP > 1)
    other.best = lapply(2:TOP,
                        function(w)
                          plot.dbpmm(x$runs[[w]],
                                     cex = .6 * cex) +
                          labs(title = paste("Solution rank:", w)))
  
  
  # =-=-=-=-=-=-=-=-
  # Final layout
  # =-=-=-=-=-=-=-=-
  bottom = ggarrange(rank, GOFP, ncol = 2, nrow = 1)
  
  left_panel = ggarrange(best,
                         bottom,
                         nrow = 2,
                         ncol = 1)
  
  
  k = ceiling(sqrt(length(other.best)))
  
  right_panel = ggarrange(plotlist = other.best,
                          ncol = k,
                          nrow = k)
  
  figure = ggarrange(
    left_panel,
    right_panel,
    widths = c(1.5, 1),
    nrow = 1,
    ncol = 2
  )
  
  figure = annotate_figure(figure,
                           top = text_grob(
                             "MOBSTER model selection",
                             color = "black",
                             face = "bold",
                             size = 18 * cex
                           ))
  
  return(figure)
}


#' Plot the goodness of fit.
#'
#' @description Plot the SSE (sum of squared error) as a proxy for the
#' goodness of fit. For the \code{TOP} solutions the SSE trace is shown.
#'
#' @param x A list of fits computed via \code{mobster_fit}.
#' @param TOP The top set of fits to use.
#' @param cex Cex of the plot.
#'
#' @return A plot of the  goodness of fit.
#' @export
#'
#' @examples
#' TODO
plot_gofit = function(x, TOP = 5, cex = 1)
{
  stopifnot(is_list_mobster_fits(x))
  
  binning = 1e-2
  
  points = lapply(x$runs, function(w)
    .compute_fit_sqerr(w, binning = binning))
  points = lapply(seq_along(points), function(w) {
    points[[w]]$K = x$fits.table$K[w]
    points[[w]]$tail = x$fits.table$tail[w]
    points[[w]]$run = w
    
    
    points[[w]]
  })
  
  points = points[1:TOP]
  points = Reduce(bind_rows, points)
  
  
  ggplot(points, aes(
    x = x,
    y = cum.y,
    fill = factor(run),
    color = factor(run)
  )) +
    geom_line(show.legend = TRUE) +
    theme_light(base_size =  8 * cex) +
    # guides(colour = FALSE) +
    xlab('VAF') +
    ylab("SSE") +
    labs(title = bquote(bold('Goodness of fit (SSE)')),
         subtitle = paste0('Binwidth = ', binning)) +
    guides(fill = FALSE, color = guide_legend(title = "Solution")) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      panel.background = element_rect(fill = 'white')
    )
}

#' Plot the latent variables of the mixture.
#'
#' @param x An object of class 'dbpmm'.
#' @param cex Cex of the plot.
#'
#' @return A plot of the latent variables of the mixture.
#' @export
#'
#' @examples
#' TODO
plot_latent_variables = function(x, cex = 1)
{
  stopifnot(inherits(x, "dbpmm"))
  
  # Reshape and cut
  lv = reshape2::melt(x$z_nk[order(x$data$cluster),])
  
  lv$value = cut(lv$value,
                 breaks = c(-Inf, seq(0, 1, 0.05), Inf))
  
  colnames(lv) = c('Point', "Cluster", "Value")
  
  colors = colorRampPalette(RColorBrewer::brewer.pal(8, 'YlGnBu'))(20)
  
  ggplot(lv, aes(x = Cluster, y = Point, fill = Value)) +
    geom_raster() +
    scale_fill_manual(values = colors) +
    theme_light(base_size = 8 * cex) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      legend.text = element_text(size = 8 * cex),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    labs(title = bquote(bold("Latent variables")))
}





is_list_mobster_fits = function(x)
{
  nOK = all(c('fits.table', 'runs', 'best', 'model.selection') %in% names(x))
  lOK = is.list(x)
  mF = all(sapply(x$runs, class) == 'dbpmm')
  mbF = (class(x$best) == 'dbpmm')
  
  all(nOK, lOK, mF, mbF)
}

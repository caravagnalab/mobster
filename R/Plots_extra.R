#' Plot the initial density of a fit.
#'
#' @param x An object of class \code{"dbpmm"}.
#' @param colors If provided, these colours will be used for each cluster.
#' If a subset of colours is provided, palette Set1 from \code{RColorBrewer} is used.
#' By default the tail colour is provided as 'gainsboro'.
#'
#' @return A ggplot object for the plot.
#' @export
#'
#' @examples
#' data(fit_example)
#' plot_init(fit_example$best)
plot_init = function(x,
                     colors = c(`Tail` = 'gainsboro')
                     )
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
    labs(title = bquote("Initialization"),
         x = "Observed Frequency",
         y = "Density") +
    guides(fill = FALSE) +
    ylim(0, max(initial.densities$y)) +
    my_ggplot_theme() +
    geom_line(data = initial.densities, aes(y = y, x = x, color = cluster)) +
    guides(color = FALSE)
  
  
  return(add_color_pl(x, den_init_pl, colors))
}


#' Plot the mixing proportions of the mixture.
#'
#' @param x An object of class 'dbpmm'.
#' @param colors If provided, these colours will be used for each cluster.
#' If a subset of colours is provided, palette Set1 from \code{RColorBrewer} is used.
#' By default the tail colour is provided as 'gainsboro'.
#'
#' @return A plot of the mixing proportions of the mixture.
#' @export
#'
#' @examples
#' data(fit_example)
#' plot_mixing_proportions(fit_example$best)
plot_mixing_proportions = function(x,                       
                                   colors = c(`Tail` = 'gainsboro')
                                   )
{
  stopifnot(inherits(x, "dbpmm"))
  
  Proportions = x$Clusters %>%
    dplyr::filter(type == 'Mixing proportion')
  
  Proportions$fit.value = round(Proportions$fit.value, 2)
  
  if (!x$fit.tail)
    Proportions = Proportions %>% filter(cluster != 'Tail')
  
  pl = ggplot(data = Proportions, aes(x = cluster, y = fit.value, fill = cluster)) +
    geom_bar(stat = "identity", width = 0.3) +
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
      size = 2.5
    ) +
    labs(title  = bquote('Mixing Proportions')) +
    xlab("") +
    ylab(bquote('Proportions (' * pi * ')')) +
    guides(fill = FALSE) +
    my_ggplot_theme() +
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
#'
#' @return A ggplot figure with the scores for model selection.
#' @export
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' plot_fit_scores(fit_example)
plot_fit_scores = function(x)
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
             y = value)) +
    geom_line(show.legend = FALSE, size = .3) +
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
    scale_size(range = c(min(K_vals), max(K_vals)) * 0.7,
               breaks = min(K_vals):max(K_vals)) +
    facet_wrap( ~ variable, ncol = 2) +
    labs(
      title  = bquote('Scores for model selection'),
      subtitle = bquote(.(nrow(scores)) ~ 'runs, '~ .(model.selection) ~ "used"),
      x = 'Model rank',
      y = 'Score'
    ) +
    guides(
      fill = FALSE,
      shape = guide_legend(title = 'Tail'),
      colour = guide_legend(title = 'Score')
    ) +
    my_ggplot_theme()
  
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
                                ...)
{
  stopifnot(is_list_mobster_fits(x))
  
  # =-=-=-=-=-=-=-=-
  # Best model fit
  # =-=-=-=-=-=-=-=-
  best = plot.dbpmm(x$best)
  
  # =-=-=-=-=-=-=-=-
  # Scores rank best
  # =-=-=-=-=-=-=-=-
  rank = plot_fit_scores(x)
  
  # =-=-=-=-=-=-=-=-
  # Goodness of fit
  # =-=-=-=-=-=-=-=-
  
  # SSE plot
  GOFP = plot_gofit(x, TOP = TOP)
  
  # =-=-=-=-=-=-=-=-=-=-
  # Low-rank solutions
  # =-=-=-=-=-=-=-=-=-=-
  
  # Other solutions ranked below top best -- maximum TOP 4 fixed
  TOP.plot = min(4, nrow(x$fits.table))
  
  other.best = NULL
  
  if (TOP > 1)
    other.best = lapply(2:TOP,
                        function(w)
                          plot.dbpmm(x$runs[[w]]) +
                          labs(title = paste("Solution rank:", w)) +
                          my_ggplot_theme(cex = .6)
                        )
  
  
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
                             size = 18
                           ))
  
  return(figure)
}


#' Plot the goodness of fit.
#'
#' @description Plot the SSE (sum of squared error) as a proxy for the
#' goodness of fit. For the \code{TOP} available solutions the SSE trace is shown.
#'
#' @param x A list of fits computed via \code{mobster_fit}.
#' @param TOP The top set of fits to use, if more than the one available
#' only the ones in \code{x} are used.
#'
#' @return A plot of the goodness of fit.
#' @export
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' plot_gofit(fit_example),
#' 
#' # This will subset the call to the one available
#' plot_gofit(fit_example, TOP = 100)
plot_gofit = function(x, TOP = 5)
{
  stopifnot(is_list_mobster_fits(x))
  
  binning = 1e-2
  
  nsol = length(x$runs)
  
  if(nsol < TOP) message("Required TOP-", TOP, " solutions, but only ", nsol, " are available.")
  TOP = min(TOP, nsol)
  
  # Compute the SSE
  points = lapply(x$runs, function(w)
    suppressWarnings(mobster:::.compute_fit_sqerr(w, binning = binning)))
  
  points = lapply(seq_along(points), function(w) {
    points[[w]]$K = x$fits.table$K[w]
    points[[w]]$tail = x$fits.table$tail[w]
    points[[w]]$run = w
    
    
    points[[w]]
  })
  
  points = points[1:TOP]
  points = Reduce(bind_rows, points)
  
  points$run = paste0('Solution #', points$run)
  ggplot(points, aes(
    x = x,
    y = y,
    fill = factor(run)
    # color = factor(run)
  )) +
    geom_line(show.legend = TRUE, size = .3) +
    my_ggplot_theme() +
    # guides(colour = FALSE) +
    ylab("SSE") +
    labs(title = bquote('Goodness of fit'),
         x = "Observed Frequency",
         subtitle = paste0('TOP-', TOP, ' solutions')) +
    guides(fill = FALSE, color = FALSE) +
    facet_wrap(~run, ncol = 1)
}

#' Plot the latent variables of the mixture.
#' 
#' @description It renders a heatmap where the latent variables
#' (reponsibilities) are shown and colured according to their value.
#' This function also calls function \code{Clusters}, using a parameter
#' that determines if a point is not to be assigned its best cluster
#' based on a cutoff.
#' 
#'
#' @param x A MOBSTER fit.
#' @param cutoff_assignment The parameter used to call function
#' \code{Clusters}, which does not assign a point to its best cluster
#' if the value of the corresponding latent variable is not above the cutoff.
#'
#' @return A plot of the latent variables of the mixture.
#' 
#' @export
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' plot_latent_variables(fit_example)
#' plot_latent_variables(fit_example, cutoff_assignment = .9)
plot_latent_variables = function(x, cutoff_assignment = 0)
{
  stopifnot(inherits(x, "dbpmm"))
  
  # assignments
  assignments = Clusters(x, cutoff_assignment) %>%
    dplyr::select(VAF, cluster) %>%
    data.frame(stringsAsFactors = FALSE)
  
  not_assign = is.na(assignments$cluster)
  n = sum(not_assign)
  p = (n/nrow(assignments)) * 100 
  ordering = order(assignments$cluster, na.last = TRUE)
  
  # Reshape and cut
  lv = reshape2::melt(x$z_nk[ordering, , drop = FALSE])
  
  # lv$value = cut(lv$value,
  #                breaks = c(-Inf, seq(0, 1, 0.05), Inf))
  # 
  colnames(lv) = c('Point', "Cluster", "Value")
  
  # 
  # colors = colorRampPalette(RColorBrewer::brewer.pal(8, 'YlGnBu'))(20)
  cuts_below = cuts_below = c()
  
  if(cutoff_assignment - 0.05 > 0) cuts_below  = seq(0, cutoff_assignment - 0.05, 0.05)
  if(cutoff_assignment - 0.05 < 1) cuts_above  = seq(cutoff_assignment, 1, 0.05)
  
  lv$Value = cut(lv$Value, breaks = c(-Inf, cuts_below, cuts_above, Inf))
  
  lblues = RColorBrewer::brewer.pal(5, 'Blues')
  lreds = RColorBrewer::brewer.pal(5, 'Reds')
  lreds = c('darkorange2', 'darkred')
  
  colors_below = colorRampPalette(lblues[1:3])(length(cuts_below))
  colors_above = colorRampPalette(lreds)(length(cuts_above))
  
  ggplot(lv, aes(x = Cluster, y = Point, fill = Value)) +
    geom_raster() +
    scale_fill_manual(values = c(colors_below, colors_above)) +
    my_ggplot_theme() +
    theme(
      legend.text = element_text(size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    guides(fill = guide_legend('')) +
    labs(
      title = bquote("Latent variables"),
      subtitle = bquote(
        .(n) ~" non assignable ("* .(p) *'%) with cutoff ' * z['nk']  ~' > ' * .(cutoff_assignment) ),
      y = paste0("Points (n =", x$N, ')')
      )
}

is_list_mobster_fits = function(x)
{
  nOK = all(c('fits.table', 'runs', 'best', 'model.selection') %in% names(x))
  lOK = is.list(x)
  mF = all(sapply(x$runs, class) == 'dbpmm')
  mbF = (class(x$best) == 'dbpmm')
  
  all(nOK, lOK, mF, mbF)
}

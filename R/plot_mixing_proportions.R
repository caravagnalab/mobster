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
  is_mobster_fit(x)

  if(is_mobsterL(x)){
    Proportions = x$Clusters %>%
      dplyr::filter(type == 'Mixing proportion')

    Proportions$fit.value = Proportions$fit.value
    Proportions$karyotype <- ""

  } else if (is_mobsterhL(x)) {
    Proportions = clonality_interpreter(x) %>%  rename(fit.value = "mixing")

  }


  if (!has_tail(x))
    Proportions = Proportions %>% filter(cluster != 'Tail')

  pl = ggplot(data = Proportions, aes(x = cluster, y = round(fit.value, 2), fill = cluster)) +
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
    labs(title  = bquote('Mixing proportions')) +
    xlab("") +
    ylab(bquote('Proportions (' * pi * ')')) +
    guides(fill = FALSE) +
    mobster:::my_ggplot_theme() +
    ylim(c(0, 1)) + facet_wrap(~ karyotype, scales = 'free_y')

  mobster:::add_fill_color_pl(x$data$cluster %>%  unique(), pl, colors)
}

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
#' @importFrom reshape2 melt
#' 
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' plot_fit_scores(fit_example)
plot_fit_scores = function(x)
{
  is_list_mobster_fits(x)
  
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
      shape = guide_legend(title = ''),
      colour = guide_legend(title = 'Score')
    ) +
    my_ggplot_theme()
  
}
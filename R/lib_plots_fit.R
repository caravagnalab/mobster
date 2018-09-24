getColors_model = function(x, alpha = 1, palette = 'Set1', tail.color = 'darkgray')
{
  labels.betas = x$Clusters %>%
    dplyr::filter(type == 'Mixing proportion', cluster != 'Tail') %>%
    dplyr::select(cluster) %>%
    dplyr::pull()
  
  # brewer palettes
    colors = c(
      `Tail` = tail.color,
      scols(labels.betas, palette = palette)
    )

  colors = ggplot2::alpha(colors, alpha)
  names(colors) = c('Tail', labels.betas)
  
  colors
}

getColors_model = function(x, alpha = 1, palette = 'Set1', tail.color = 'darkgray')
{
  labels.betas = x$Clusters %>%
    dplyr::filter(type == 'Mixing proportion', cluster != 'Tail') %>%
    dplyr::select(cluster) %>%
    dplyr::pull()
  
  # brewer palettes, plus tail color
  colors = scols(labels.betas, palette = palette)
  if(x$fit.tail) colors = c(`Tail` = tail.color, colors)
   
  # alpha 
  colors = ggplot2::alpha(colors, alpha)

  # name the vector 
  if(x$fit.tail) names(colors) = c('Tail', labels.betas)
  else names(colors) = labels.betas

  colors
}

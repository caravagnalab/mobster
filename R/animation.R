#' Title
#'
#' @param x 
#' @param density 
#' @param palette 
#' @param tail.color 
#'
#' @return
#' @export
#'
#' @import plotly 
#'
#' @examples
MOBSTER_animate = function(x, density = 0.05,  
                           transition_length = 1,
                           state_length = 1,
                           lib = 'plotly',
                           palette = 'Set1', tail.color = 'gray')
{
  if(any(is.null(x$trace))) stop("Trace is null -- run fit with trace = TRUE.")
  
  # Load colors
  colors = getColors_model(x, alpha = 1, palette = palette, tail.color = tail.color)
  
  # Prepare trace
  tr = split(x$trace, f = x$trace$step)
  
  # subset steps
  steps = seq(1, length(tr), round(density * length(tr)))
  
  # per point density
  tr.points = lapply(steps, function(w, x) {
    
    new.x = x
    new.x$Clusters = tr[[w]]
    
    points = template_density(new.x, x.axis = seq(0, 1, 0.01), binwidth = 0.01, reduce = TRUE)
    points$step = w
    
    points
  },
  x = x)
  
  tr.points = Reduce(rbind, tr.points)
  
  if(lib == 'plotly')
  {
    require(plotly)
    
    p <- tr.points %>%
      plot_ly(
        x = ~x,
        y = ~y,
        frame = ~step,
        color = ~cluster,
        type = 'scatter',
        mode = 'markers',
        showlegend = TRUE
      )
    
    return(p)
  }
  
  # 
  require(ggplot2)
  require(gganimate)
  require(magrittr)
  
  movie = tr.points
  
  movie.render = ggplot(movie, aes(x, y, color = cluster)) + 
    # geom_histogram(data = x$data, aes(VAF, y = ..count../sum(..count..), color = cluster, fill = cluster), binwidth = 0.01) + 
    geom_histogram(data = x$data, aes(VAF, y = ..count../sum(..count..), color = NULL, fill = NULL), binwidth = 0.01) + 
    geom_point() +
    labs(
      title = "MOBSTER fit animation",
      x = "Observed Frequency", y = "Density") +
    scale_color_manual(values = colors, labels = names(colors)) +
    scale_fill_manual(values = colors, labels = names(colors)) +
    transition_states(
      step,
      transition_length = 1,
      state_length = 1
    ) +
    enter_fade() + 
    exit_shrink() +
    ease_aes('sine-in-out')
  
  movie.render
}
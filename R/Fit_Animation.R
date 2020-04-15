#' Render the animation of a MOBSTER fit
#'
#' @param x An object of class \code{"dbpmm"} computed with 
#' \code{trace = TRUE} in \code{mobster_fit}.
#' @param density The density of steps in the trace that should be used
#' to render the animation. For instance, 5% would mean that if there are
#' n points in the trace, one every n*0.05 is used to animate the fit.
#' @param lib Library used for renderining; this can be either 
#' \code{'plotly'} or \code{'gganimate'}.
#'
#' @return The animation object, according to the chosen  \code{lib}.
#' 
#' @export
#'
#' @import plotly
#' @import gganimate
#'
#' @examples
#' data(fit_example)
#' 
#' model = mobster_fit(
#' fit_example$best$data, 
#' parallel =  FALSE, 
#' samples = 3,
#' init = 'random',
#' trace = TRUE,
#' K = 2,
#' tail = TRUE)
#' 
#' fit_animation(model$best, lib = 'plotly')
fit_animation = function(x,
                         density = 0.05,
                         lib = 'plotly')
{
  stopifnot(inherits(x, "dbpmm"))
  if (any(is.null(x$trace)))
    stop("Trace is null -- run fit with trace = TRUE.")
  
  # Prepare trace
  tr = split(x$trace, f = x$trace$step)
  
  # subset steps
  steps = seq(1, length(tr), round(density * length(tr)))
  
  # per point density
  tr.points = lapply(steps, function(w, x) {
    new.x = x
    new.x$Clusters = tr[[w]]
    
    points = template_density(
      new.x,
      x.axis = seq(0, 1, 0.01),
      binwidth = 0.01,
      reduce = TRUE
    )
    points$step = w
    
    points
  },
  x = x)
  
  tr.points = Reduce(rbind, tr.points)
  
  if (lib == 'plotly')
  {
    require(plotly)
    
    p <- tr.points %>%
      plot_ly(
        x = ~ x,
        y = ~ y,
        frame = ~ step,
        color = ~ cluster,
        type = 'scatter',
        mode = 'markers',
        showlegend = TRUE
      )
    
    return(p)
  }
  
  movie = tr.points
  
  require(gganimate)
  
  movie.render = ggplot(movie, aes(x, y, color = cluster)) +
    # geom_histogram(data = x$data, aes(VAF, y = ..count../sum(..count..), color = cluster, fill = cluster), binwidth = 0.01) +
    geom_histogram(
      data = x$data,
      aes(
        VAF,
        y = ..count.. / sum(..count..),
        color = NULL,
        fill = NULL
      ),
      binwidth = 0.01
    ) +
    geom_point() +
    labs(title = "MOBSTER fit animation",
         x = "Observed Frequency", 
         y = "Density") +
    transition_states(step,
                      transition_length = 1,
                      state_length = 1) +
    enter_fade() +
    exit_shrink() +
    ease_aes('sine-in-out')
  
  movie.render
}

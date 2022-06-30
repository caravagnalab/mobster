#' Plot the NLL trace.
#' 
#' @description It plots the Negative log-Likelihood (NLL) trace of a fit.
#' 
#' @param x A MOBSTER fit.
#'
#' @return A plot of the the NLL trace for the input fit.
#' 
#' @export
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' plot_NLL(fit_example$best)
plot_NLL = function(x)
{
  is_mobster_fit(x)
  
  NLL = data.frame(
    x = 1:length(x$all.NLL), 
    NLL = x$all.NLL
  )
  
  steps = nrow(NLL)
  
  status = paste(
    ifelse(x$fit.type == 'MM', 'Moments Matching', 'Maximum Likelihood'),
    ifelse(x$status, 'converged', 'NOT converged'),
    'in', steps, 'steps'
  )
  
  ggplot(NLL[-1,], aes(x = x, y = NLL)) +
    geom_line(color = 'steelblue') +
    geom_point(size = 1.5) +
    my_ggplot_theme() +
    labs(
      title = bquote("Negative log-Likelihood"),
      subtitle = bquote(.(status)),
      x = 'Step'
    ) 
  
}
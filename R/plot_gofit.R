#' Plot the goodness of fit.
#'
#' @description Plot the SSE (sum of squared error) as a proxy for the
#' goodness of fit. For the \code{TOP} available solutions the SSE trace is shown.
#'
#' @param x A list of fits computed via \code{mobster_fit}.
#' @param TOP The top set of fits to use, if more than the one available
#' only the ones in \code{x} are used.
#' @param nx Columns in the matrix layout.
#' @param ny Rows in the matrix layout.
#' @param TOP The first \code{TOP} fits are used.
#' @return A plot of the goodness of fit.
#' @export
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' plot_gofit(fit_example, TOP = 3)
plot_gofit = function(x, TOP = 6, nx = 3, ny = 2)
{
  is_list_mobster_fits(x)
  
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
    ylab("SSE") +
    labs(title = bquote('Goodness of fit'),
         x = "Observed Frequency",
         subtitle = paste0('TOP-', TOP, ' solutions')) +
    guides(fill = FALSE, color = FALSE) +
    facet_wrap(~run, ncol = nx, nrow = ny)
}
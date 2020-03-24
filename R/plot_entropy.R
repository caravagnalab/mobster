#' Plot the entropy of a MOBSTER mixture.
#' 
#' @description Returns a plot of the entropy and the reduced entropy
#' for this mixture, binning the domain with bins of size 1e-3. The
#' full entropy is coloured with a gradient, and the reduced is dashed.
#'
#' @param x An object of class \code{"dbpmm"}.
#'
#' @return A ggplot object for the plot.
#' 
#' @export
#'
#' @examples
#' data(fit_example)
#' plot_entropy(fit_example$best)
plot_entropy = function(x)
{
  domain = seq(0, 1, 0.001)

  # LV
  N = length(domain)
  z_nk  = matrix(0, nrow = N, ncol = x$K)
  colnames(z_nk) = names(x$pi)
  
  pdf.w = z_nk
  
  # Get log density per component
  for (k in 1:x$K)
    pdf.w[, k] = ddbpmm(x, data = domain, components = paste(k), log = TRUE)
  
  # Calculate probabilities using the logSumExp trick for numerical stability
  Z = apply(pdf.w, 1, .log_sum_exp)
  z_nk   = pdf.w - Z
  z_nk   = apply(z_nk, 2, exp)
  
  # Full entropy
  entropy_trace = data.frame(
    x = domain,
    y = - apply(z_nk * log(z_nk), 1, sum, na.rm = TRUE)
  )
  
  # Reduced entropy - remove 1 LV column
  cz_nk = z_nk[, 2:ncol(z_nk), drop = FALSE]
  
  # This is un-normalized -- we compute the empirical normalizing constant (C)
  C = rowSums(cz_nk)
  for (i in 1:nrow(cz_nk))
    cz_nk [i, ] = cz_nk [i, ] / C[i]
  
  # Same as above
  reduced_trace = data.frame(
    x = domain,
    y = - apply(cz_nk * log(cz_nk), 1, sum, na.rm = TRUE)
  )
  
  ggplot(entropy_trace,  aes(x = x, y = y, color = y)) +
    geom_line(size = 1) +
    geom_line(data = reduced_trace, color = 'darkred', linetype = 'dashed', size = .3) +
    mobster:::my_ggplot_theme() +
    # scale_color_distiller(palette = 'Spectral', direction = -1) +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = 'Entropy of the latent variables',
      subtitle = 'Dashed line: reduced entropy',
      y = 'Value',
      x = 'Observed Frequency'
    ) +
    guides(color = guide_colourbar("Entropy  ", barwidth = unit(3, 'cm')))
}
  


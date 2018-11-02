
plot_2D_cluster_density = function(x, s.x, s.y, b.c, col = 'red', cov = 100, cutoff = 0.05)
{
  binfit = x$fit.Binomial
  rownames(binfit$theta_k) = x$samples

  # mixing proportion
  pi = binfit$pi_k[b.c]

  # Binomial parameters
  p.x = binfit$theta_k[s.x, b.c]
  p.y = binfit$theta_k[s.y, b.c]

  # less than 1% do not plot anything
  eplot = function() {
    ggplot() + annotation_custom(text_grob(paste0('< ', cutoff), size = 10)) + theme(element_blank(), axis.text = element_blank())
  }

  if(all(p.x < cutoff, p.y < cutoff)) return(eplot())

  # Otherwise plot
  domain = seq(0, 1, 0.01)

  grid = expand.grid(domain, domain, stringsAsFactors = FALSE)
  colnames(grid) = c('x', 'y')

  grid$nv.x = round(grid$x * cov)
  grid$nv.y = round(grid$y * cov)

  mpdf = function(x, y, b.c) {
    p.x = binfit$theta_k[s.x, b.c]
    p.y = binfit$theta_k[s.y, b.c]

    pi * dbinom(x, size = cov, prob = p.x) *
      dbinom(y, size = cov, prob = p.y)
  }

  all.points = lapply(
    b.c,
    # colnames(binfit$theta_k),
    function(cl)
    {
      cbind(grid,
            group = cl,
            pdf =  apply(grid, 1, function(w)
              mpdf(w['nv.x'], w['nv.y'], b.c = cl)))
    })

  all.points = Reduce(rbind, all.points)

  ggplot(all.points, aes(x = x, y = y)) +
    geom_raster(aes(fill = pdf), interpolate = TRUE) +
    theme_light(base_size = 10) +
    facet_wrap(~group) +
    scale_fill_gradient (low = 'gainsboro', high = col) +
    theme(legend.position="none") +
    labs(x = s.x, y = s.y) +
    geom_hline(yintercept = cutoff, size = .3, color = 'red',linetype = "longdash") +
    geom_vline(xintercept = cutoff, size = .3, color = 'red',linetype = "longdash")
}


binomial2D_cluster_template_density = function(x, s.x, s.y, col = 'red', cov = 100, cutoff = 0.05)
{
  if (is.null(x$fit.Binomial))
    stop("Binomial clusters are not available!")
  
  binfit = x$fit.Binomial
  rownames(binfit$theta_k) = x$samples
  
  # mixing proportion
  pi = binfit$pi_k
  
  # purity (we have to adjust the binomial parameters for purity if we want to plto with adjusted VAF)
  purity = x$purity
  
  # Otherwise plot 
  domain = seq(0, 1, 0.02)
  
  grid = expand.grid(domain, domain, stringsAsFactors = FALSE)
  colnames(grid) = c('x', 'y')
  
  grid$nv.x = round(grid$x * cov)
  grid$nv.y = round(grid$y * cov)
  
  mpdf = function(x, y, b.c) {
    p.x = binfit$theta_k[s.x, b.c] / purity[s.x]
    p.y = binfit$theta_k[s.y, b.c] / purity[s.y]
    
    pi[b.c] * dbinom(x, size = cov, prob = p.x) *
      dbinom(y, size = cov, prob = p.y) 
  }
  
  all.points = lapply(
    colnames(binfit$theta_k),
    function(cl)
    {
      cbind(grid,
            group = cl,
            pdf =  apply(grid, 1, function(w)
              mpdf(w['nv.x'], w['nv.y'], b.c = cl)))
    })
  
  Reduce(rbind, all.points)
}


#' Title
#'
#' @param x 
#' @param palette 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
mobster_plt_2D_cluster_density = function(x, cov = 100, palette = "Set1", ...) {
  stopifnot(!is.null(x$fit.Binomial))
  
  # Get colors
  colors = mobster:::scols(names(x$fit.Binomial$pi_k), palette = palette)
  
  # Clusters names
  cls = x$fit.Binomial$pi_k
  
  # paired of samples
  pairs = combn(x$samples, 2)
  
  pl_list = NULL
  
  for(i in 1:ncol(pairs))
  {
    # Get a plot for each cluster
    pl = lapply(
      seq_along(cls), 
      function(c) plot_2D_cluster_density(
        x = x, 
        b.c = names(cls)[c], 
        s.x = pairs[1, i], 
        s.y = pairs[2, i], 
        col = colors[c], 
        ...)
    )
    
    header_pl = ggplot() + 
      annotation_custom(text_grob(paste0(pairs[1, i], ' vs ', pairs[2, i]), size = 20)) + 
      theme(element_blank(), axis.text = element_blank()) 
    
    
    # Display it as a row
    pl_fig = ggarrange(plotlist = append(list(header_pl), pl), ncol = length(cls) + 1, nrow = 1)
    pl_list = append(pl_list, list(pl_fig))
  }
  
  fg = ggarrange(plotlist = pl_list, ncol = 1, nrow = ncol(pairs))
  annotate_figure(fg, top = paste("Multivariate Binomial density for", x$mobster_global$description), 
                  bottom = paste0("Density at coverage", cov, "x"))
}

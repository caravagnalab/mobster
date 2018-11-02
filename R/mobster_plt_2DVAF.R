#' Title
#'
#' @param data
#' @param x
#' @param y
#' @param cluster
#' @param marginal
#'
#' @return
#' @export
#'
#' @examples
mobster_plt_2DVAF = function(obj,
                             x,
                             y,
                             cluster = NULL,
                             cluster.label = 'cluster',
                             marginal = FALSE,
                             density = TRUE,
                             cex = 1,
                             alpha = 0.3,
                             downsample = 1000,
                             scales.fixed = TRUE,
                             # VAF.range = c(0.05, 0.6),
                             marginal.remove.zeroes = T,
                             palette = 'Set1') 
{
  caption = ""
  
  ####
  # Prepare input VAF data. Match clusters if required
  # and downsample if there are too many entries
  data = VAF_table(obj, samples = c(x, y), suffix = '') 
  
  if (!is.null(cluster)) {
    
    stopifnot(cluster.label %in% colnames(cluster))
    stopifnot(is.data.frame(cluster))
    
    data = full_join(data, cluster, by = 'id')    
  }

  # downsample data if too many
  if(!is.null(downsample) & nrow(data) > downsample)
  {
    stopifnot(is.numeric(downsample))
    data = data[sample(1:nrow(data), downsample), , drop = FALSE]
    
    caption = paste0(caption, "Downsampled (N = ", downsample, ')')
  }
  
  data = data  %>% select(-id)
  
  ####
  # Ggplot object, depends on clustering or not.  With clustering, we set colour
  # of the points according to the labels. If it is not, we use a density-depenednt coloring
  if (!is.null(cluster))
  { 
    labels = unique(data[, cluster.label]) %>% pull(!!cluster.label)
    colors = mobster:::scols(sort(labels), palette)
    
    p = ggplot(data = data,
               aes(
                 x = eval(parse(text = x)),
                 y = eval(parse(text = y)),
                 color = factor(eval(parse(text = cluster.label)))
               )) +
      scale_color_manual(values = colors) +
      geom_point(alpha = alpha, size = 1 * cex)
  }
  else {
    # Get density of points in 2 dimensions.
    # @param x A numeric vector.
    # @param y A numeric vector.
    # @param n Create a square n by n grid to compute density.
    # @return The density within each square.
    get_density <- function(x, y, n = 100) {
      dens <- MASS::kde2d(x = x, y = y, n = n)
      ix <- findInterval(x, dens$x)
      iy <- findInterval(y, dens$y)
      ii <- cbind(ix, iy)
      return(dens$z[ii])
    }
    
    require(viridis)
    
    # annotated density
    data$density =  get_density(data %>% pull(!!x), data %>% pull(!!y))
    
    p = ggplot(data,
               aes(x = eval(parse(text = x)),
                   y = eval(parse(text = y)))) +
      scale_color_viridis() +
      geom_point(alpha = alpha, aes(color = density), size = 2 * cex)
  }
  
  ####
  # Add density points if required. Density plot can only go with clustering  
  if(density)
  {
    if(is.null(cluster)) stop("Density requires also Binomial clusters")
  
    # Get median coverage of the data
    den.coverage = DP(obj, samples = c(x,y)) %>% pull(value)
    den.coverage = round(median(den.coverage, na.rm = TRUE))
    
    # Compute the density
    density.points = binomial2D_cluster_template_density(obj, x, y)

    p = p + 
      geom_contour(data = density.points, 
                   aes(x = x, y = y, z = pdf, color = group),
                   inherit.aes = F) 
    
    caption = paste0(caption, 
                     '\n', den.coverage, "x density adjusted for purity (",
                     x, ' = ', obj$purity[x], '; ',
                     y, ' = ', obj$purity[y], 
                     ")")
  }
  
  ####
  # Adjust plot ranges if required 
  if (scales.fixed & max(data[, x], na.rm = T) < 1) p = p + xlim(0, 1)
  if (scales.fixed & max(data[, y], na.rm = T) < 1) p = p + ylim(0, 1)
  
  ####
  # Plot style
  p = p +
    geom_vline(xintercept = 0, colour = "darkgray", size = .3) +
    geom_hline(yintercept = 0, colour = "darkgray", size = .3) +
    labs(
      title = bquote(bold(.(x)) ~ "vs" ~ bold(.(y))),
      caption = caption,
      x = x,
      y = y
    ) +
    guides(color = guide_legend(title = cluster.label)) +
    theme_light(base_size = 8 * cex) +
    theme(
      panel.background = element_rect(fill = alpha("gray95", 1)),
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      legend.text = element_text(size = 8 * cex)
    ) 
    # geom_vline(xintercept = max(VAF(obj, samples = x) %>% pull(value)),
    #            colour = "black",
    #            linetype = "longdash", size = .3) +
    # geom_hline(yintercept = max(VAF(obj, samples = y) %>% pull(value)),
    #            colour = "black",
    #            linetype = "longdash", size = .3) +
    # geom_vline(xintercept = min(VAF(obj, samples = x) %>% pull(value)),
    #            colour = "black",
    #            linetype = "longdash", size = .3) +
    # geom_hline(yintercept = min(VAF(obj, samples = y) %>% pull(value)),
    #            colour = "black",
    #            linetype = "longdash", size = .3)
    # 
  
  
  ####
  # Return the object or add the marginal histogram
  if (!marginal) return(p)
  
  data = data[data[, x] > 0,]
  data = data[data[, y] > 0,]
  
  require(ggExtra)
  
  ggMarginal(
    p,
    fill = "black",
    binwidth = 0.01,
    alpha = 1,
    aes = aes(
      x = eval(parse(text = x)),
      colour = eval(parse(text = cluster)),
      y = eval(parse(text = y))
    ),
    data = data,
    type = "histogram",
    xparams = list(colour = "black", size = 0.1),
    yparams = list(colour = "black", size = 0.1)
  )
  
}

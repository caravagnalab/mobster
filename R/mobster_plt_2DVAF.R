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
                             cex = 1,
                             # VAF.range = c(0.05, 0.6),
                             marginal.remove.zeroes = T,
                             palette = 'Set1') 
{
  # pio::pioHdr(
  #   header = 'MOBSTER -- plot 2D VAF',
  #   prefix = '\t-',
  #   toPrint = c(
  #     `Sample IDs` = paste(x, y, sep = ' -vs- '),
  #     `VAF range` = paste(VAF.range, collapse = ' -- ')
  #   )
  # )
  
  data = VAF_table(obj, samples = c(x, y), suffix = '') 
  
  if (!is.null(cluster)) {
    
    stopifnot(cluster.label %in% colnames(cluster))
    stopifnot(is.data.frame(cluster))
    
    data = full_join(data, cluster, by = 'id')    
  }

  data = data  %>% select(-id)

  require(ggplot2)
  
  # data = data.frame(data)
  
  if (!is.null(cluster))
    p = ggplot(data, aes(
      x = eval(parse(text = x)),
      colour = factor(eval(parse(text = cluster.label))),
      y = eval(parse(text = y))
    )) +
    scale_color_brewer(palette = palette, drop = FALSE) +
    geom_point(alpha = 0.1, size = .8)
  
  if (is.null(cluster))
    p = ggplot(data, aes(
      x = eval(parse(text = x)),
      colour = 'Unclustered',
      y = eval(parse(text = y))
    )) +
    geom_point(alpha = 0.1, colour = 'black', size = .8)
  
  if (max(data[, x], na.rm = T) < 1)
    p = p + xlim(0, 1)

  if (max(data[, y], na.rm = T) < 1)
    p = p + ylim(0, 1)
  
  
  p = p +
    # theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      # panel.grid.minor = element_blank(),
      axis.line = element_line(
        size = 0.5,
        linetype = "solid",
        colour = "white"
      )
    ) +
    geom_vline(xintercept = 0, colour = "darkgray", size = .3) +
    geom_hline(yintercept = 0, colour = "darkgray", size = .3) +
    # geom_density_2d(aes(fill = ..level..), geom = "polygon", alpha = .5) +
    labs(
      title = paste(x, "vs", y),
      # subtitle = paste0("",
      # caption = paste0("Dashed cutoffs (", VAF.range[1], ', ', VAF.range[2], ')'),
      x = x,
      y = y
    ) +
    guides(colour = guide_legend(title = cluster.label)) +
    theme_light(base_size = 8 * cex) +
    theme(
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
  
  
  if (!marginal)
    return(p)
  
  data = data[data[, x] > 0,]
  data = data[data[, y] > 0,]
  
  require(ggExtra)
  ggMarginal(
    p,
    fill = "gainsboro",
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

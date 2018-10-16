# Plot
#' Title
#'
#' @param data
#' @param cluster
#' @param title
#' @param cex
#' @param cutoff
#' @param palette
#' @param free
#'
#' @return
#' @export
#'
#' @examples
MOBSTER_mplot_mixing = function(data, cluster, title = "Mixing proportions", cex = 1, cutoff = 0.05, palette = 'Set1', free = FALSE)
{
  pi.df = get_mixing_proportions(data, cluster)
  pi.df = reshape2::melt(pi.df, id = 'Cluster')
  pi.df$value = round(pi.df$value, 2)

  pi.plot = ggplot(pi.df, aes(value, x = Cluster, fill = factor(Cluster))) +
    geom_bar(stat = "identity", position="dodge") +
    scale_fill_brewer(palette = palette) +
    theme_light(base_size = 8 * cex) +
    guides(fill=guide_legend(title="Cluster"))  +
    facet_wrap(~variable) +
    labs(
      title = bquote(bold(.(title))),
      y = bquote("Mixing proportions (" ~pi~")"),
      subtitle = 'Multivariate clustering of read counts'
    ) +
    geom_hline(yintercept = cutoff, colour = 'red', linetype = "longdash") +
    theme(legend.position="bottom",
          legend.text = element_text(size = 8 * cex),
          legend.key.size = unit(.3 * cex, "cm")
    ) +
    geom_text(data = pi.df[pi.df$value > 0.05, ], aes(label=value), vjust=1.6, color="white", size = 2)

  if(!free) pi.plot = pi.plot + coord_cartesian(ylim = c(0,1))

  pi.plot
}

#' Title
#'
#' @param data
#' @param cluster
#' @param VAF.samples
#' @param title
#' @param cex
#' @param cutoff
#' @param palette
#' @param free
#'
#' @return
#' @export
#'
#' @examples
MOBSTER_mplot_peaks = function(data, cluster, VAF.samples, title = "Clone peaks", cex = 1, cutoff = 0.05, palette = 'Set1', free = FALSE)
{
  labels = data[, cluster]

  # Get summary from posterior clustering assignments
  values = lapply(VAF.samples,
                  function(s){
                    sp = split(data[, s], f = labels)
                    sapply(sp, summary)
                  })

  # Get median
  values = Reduce(
    rbind,
    lapply(values, function(w)w['Median', ]))

  # Create data frame for plot
  values = as.data.frame(values)
  values = round(values, 2)

  rownames(values) = NULL
  values$Sample = samples

  # Melt
  df.scores = reshape2::melt(values, id = 'Sample')

  # Colors are taken from pi
  colors = RColorBrewer::brewer.pal(n = ncol(values) - 1, palette)

  # scale alpha accordingly
  pi.df = get_mixing_proportions(data, cluster)
  pi.df[, cluster] = pi.df[, cluster] * 2
  pi.df[pi.df[, cluster] > 1, cluster] = 1
  pi.df[, cluster] = cut(pi.df[, cluster], breaks = seq(0, 1, 0.1), labels = seq(0, 1, 0.1)[-1])
  pi.df[, cluster] = as.numeric(pi.df[, cluster])
  # levels(pi.df[, cluster])

  for(i in 1:length(colors)) colors[i] = alpha(colors[i], pi.df[i, cluster] * 0.1)

  ggplot(df.scores,
         aes(value, x = variable, fill = variable)) +
    geom_hline(yintercept = 0.5, colour = 'black', linetype = "longdash") +
    geom_bar(stat = "identity", position="dodge") +
    scale_fill_manual(values = colors) +
    facet_wrap(~Sample, nrow = 1) +
    labs(
      title = bquote(bold(.(title))),
      y = bquote("Median VAF"),
      x = bquote('Cluster (transparency based on proportions)'),
      subtitle = 'Multivariate clustering of read counts'
    ) +
    geom_hline(yintercept = cutoff, colour = 'red', linetype = "longdash") +
    theme_light(base_size = 8) +
    theme(
      # legend.position="bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      legend.text = element_text(size = 8 * cex)) +
    geom_text(data = df.scores[df.scores$value > 0, ], aes(label=value), vjust=1.6, color="black", size = 2 * cex) +
    guides(fill = guide_legend(title = "Cluster")) +
    coord_cartesian(ylim = c(0,1))

}


#' Title
#'
#' @param data
#' @param samples
#' @param VAF.range
#'
#' @return
#' @export
#'
#' @examples
MOBSTER_mplot_fits_MOBSTER_sciClone = function(
  data,
  samples,
  cluster1 = 'sciClone.cluster',
  cluster2 = 'MOBSTER.sciClone.cluster',
  palette.cluster1 = 'Set1',
  palette.cluster2 = 'Spectral',
  cex = 1,
  VAF.range = c(0.05, 0.6))
{

  pio::pioHdr("MOBSTER - Multivariate fits plot",
              c(`Cluster 1` = cluster1,
                `Cluster 2` = cluster2,
                `Samples` = paste(samples, collapse = ',')
                )
  )


  require(ggplot2)

  #  Diagonal
  MB.figure = ggpubr::ggarrange(
    plotlist = plot_diagonal_MOBSTER(data$best.MOBSTER, samples, cex = cex),
    nrow = 1,
    ncol = length(data$best.MOBSTER),
    labels = LETTERS[seq(data$best.MOBSTER)]
  )

  labels = LETTERS[-seq(data$best.MOBSTER)]

  plots = NULL
  id = 1


  for (s in seq(samples)) {
    for (w in s:length(samples)) {
      if (s != w) {

        pl.1 = plot_2DVAF(
          data$data,
          x = paste0('VAF.', samples[s]),
          y = paste0('VAF.', samples[w]),
          cluster = cluster1,
          palette = palette.cluster1,
          VAF.range = VAF.range,
          cex = cex
        )

        pl.2 = plot_2DVAF(
          data$data,
          x = paste0('VAF.projected.', samples[s]),
          y = paste0('VAF.projected.', samples[w]),
          cluster = cluster2,
          palette = palette.cluster2,
          VAF.range = VAF.range,
          cex = cex
        )

        fig = ggpubr::ggarrange(pl.1,
                                pl.2,
                                nrow = 1,
                                ncol = 2,
                                labels = labels[id])

        plots = append(plots, list(fig))
        id = id + 1
      }
    }
  }

  twoBtwo = ggpubr::ggarrange(
    plotlist = plots,
    ncol = 2,
    nrow = length(plots)/2
  )

  figure = ggpubr::ggarrange(
    MB.figure,
    twoBtwo,
    ncol = 1,
    nrow = 2,
    heights = c(.15, 1)
  )

  figure

  # plotlist = append(list(MB.figure), plots)
  # figure = ggpubr::ggarrange(
  #   plotlist = plotlist,
  #   nrow = length(plotlist),
  #   ncol = 1,
  #   common.legend = T,
  #   heights = c(.75, rep(1, length(plotlist)))
  # )
}




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
plot_2DVAF = function(data, x, y, cluster = NULL,
                      marginal = FALSE,
                      cex = 1,
                      VAF.range = c(0.05, 0.6), marginal.remove.zeroes = T,
                      palette = 'Set1') {

  pio::pioHdr(header = 'MOBSTER -- plot 2D VAF', prefix = '\t-',
              toPrint = c(
                `Sample IDs` = paste(x, y, sep = ' -vs- '),
                `VAF range` = paste(VAF.range, collapse = ' -- ')
              ))

  # data = data[!is.na(data[, cluster]), ]

  # print(cluster)
  if(!is.null(cluster)) data[, cluster] = paste(data[, cluster])

  # there is a small bug here in ggplot, gotta render rows different
  # one another to avoid a
  #
  # Error: `data` must be uniquely named but has duplicate elements
  #
  # n = sum(data[, x] == 0)
  # r = abs(runif(n)/1e9) # epsilon noise ~ O(1e-10)
  # data[data[, x] == 0, x] = r
  #
  # n = sum(data[, y] == 0)
  # r = abs(runif(n)/1e9) # epsilon noise ~ O(1e-10)
  # data[data[, y] == 0, y] = r

  require(ggplot2)

  # data = data.frame(data)

  if(!is.null(cluster))
    p = ggplot(data, aes(x = eval(parse(text = x)), colour = factor(eval(parse(text = cluster))), y = eval(parse(text = y)))) +
    scale_color_brewer(palette = palette, drop=FALSE) +
    geom_point(alpha = 0.6)

  if(is.null(cluster))
    p = ggplot(data, aes(x = eval(parse(text = x)), colour = 'Unclustered', y = eval(parse(text = y)))) +
    geom_point(alpha = 0.6, colour = 'black')

  maxX = max(data[, x], na.rm = T)
  if(maxX < 1)  p = p + xlim(0,1)


  p = p +
    # theme_minimal() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid",
                                   colour = "white")) +
    geom_vline(xintercept = 0, colour="darkgray") +
    geom_hline(yintercept = 0, colour="darkgray") +
    # geom_density_2d(aes(fill = ..level..), geom = "polygon", alpha = .5) +
    labs(
      title = paste(x, "vs", y),
      subtitle = "",
      caption = paste0("Dashed cutoffs (", VAF.range[1], ', ', VAF.range[2], ')'),
      x = x, y = y) +
    xlim(0, 1) +
    ylim(0, 1) +
    guides(colour = guide_legend(title = cluster)) +
    theme_light(base_size = 8 * cex) +
    theme(
      legend.position="bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      legend.text = element_text(size = 8 * cex)) +
    geom_vline(xintercept = VAF.range[1], colour="red", linetype = "longdash") +
    geom_hline(yintercept = VAF.range[1], colour="red", linetype = "longdash") +
    geom_vline(xintercept = VAF.range[2], colour="steelblue", linetype = "longdash") +
    geom_hline(yintercept = VAF.range[2], colour="steelblue", linetype = "longdash")


  if(!marginal) return(p)

  if(marginal.remove.zeroes) {
    data = data[data[, x] > 0, ]
    data = data[data[, y] > 0, ]
  }

  require(ggExtra)
  ggMarginal(
    p,
    fill = "gainsboro", binwidth = 0.01, alpha = 1,
    aes = aes(
      x = eval(parse(text = x)),
      colour = eval(parse(text = cluster)),
      y = eval(parse(text = y))
    ),
    data = data,
    # type = 'density',
    type = "histogram",
    xparams = list(colour = "black", size = 0.1),
    yparams = list(colour = "black", size = 0.1))

}


#' Title
#'
#' @param data
#' @param cluster
#' @param samples
#' @param palette
#' @param cex
#' @param title
#'
#' @return
#' @export
#'
#' @examples
MOBSTER_mplot_mapping2MOBSTER_clusters = function(data,
                                                  MOBSTER.fits,
                                                  cluster = 'sciClone.cluster',
                                                  samples,
                                                  palette = 'Set2',
                                                  cex = 1,
                                                  tail.color = 'gray',
                                                  title = paste0('Mapping of ', cluster, ' clusters to MOBSTER ones')
)
{
  labels = data[, cluster]
  groups = split(data, f = labels)

  values = lapply(
    groups,
    function(g) {
      lapply(samples, function(s) {

        w = which(g[, paste0('VAF.', s)] > 0)

        t = table(g[w, paste0('MOBSTER.', s, '.cluster')])
        t = as.data.frame(t)

        if(nrow(t)>0) t$Sample = s
        else t = NULL
        t
      })
    })

  values = lapply(values, Reduce, f = rbind)
  values = lapply(seq(values), function(g) {
    cbind(values[[g]], Cluster = names(values)[g])
  })

  values = Reduce(rbind, values)

  # get Palette
  maxBeta = max(sapply(MOBSTER.fits, function(w) w$Kbeta))

  col = scols(
    paste0('C', 1:maxBeta), palette = palette)
  col = c(col, `Tail` = tail.color)


  p1 = ggplot(values,
         aes(Freq, x = Sample, fill = factor(Var1))) +
    geom_bar(stat = "identity",  position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = col, labels = names(col)) +
    theme_light(base_size = 8 * cex) +
    guides(fill=guide_legend(title="MOBSTER Cluster"))  +
    # facet_grid(~Cluster, scales = "free_y", space = "free_y") +
    facet_wrap(~Cluster, scales="free_x", nrow = 1)+
    labs(
      title = bquote(bold(.(title))),
      y = bquote("Number of observations"),
      subtitle = 'Multivariate clustering of read counts'
    ) +
    # geom_hline(yintercept = cutoff, colour = 'red', linetype = "longdash") +
    theme(legend.position="bottom",
          legend.text = element_text(size = 8 * cex),
          legend.key.size = unit(.3 * cex, "cm")
    ) +
    coord_flip()

  dummy <- ggplot(values,
                  aes(Freq, x = Sample)) +
    facet_wrap(~Cluster, scales="free_x", nrow = 1) +
    geom_rect(aes(fill=factor(Cluster)), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    scale_fill_brewer(palette = 'Set2')
  theme_minimal()

  library(gtable)

  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(dummy)

  gtable_select <- function (x, ...)
  {
    matches <- c(...)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    x
  }

  panels <- grepl(pattern="panel", g2$layout$name)
  strips <- grepl(pattern="strip_t", g2$layout$name)
  g2$layout$t[panels] <- g2$layout$t[panels] - 1
  g2$layout$b[panels] <- g2$layout$b[panels] - 1

  require(grid)
  new_strips <- gtable_select(g2, panels | strips)
  # grid.newpage()
  # grid.draw(new_strips)

  gtable_stack <- function(g1, g2){
    g1$grobs <- c(g1$grobs, g2$grobs)
    g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
    g1$layout <- rbind(g1$layout, g2$layout)
    g1
  }

  ## ideally you'd remove the old strips, for now they're just covered
  new_plot <- gtable_stack(g1, new_strips)
  ggpubr::as_ggplot(new_plot)
}

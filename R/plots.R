#' Plot an object of class dbpmm. This is a set of multiple frames plot via ggplot2.
#'
#' @param x an object of class dbpmm
#' @param palette a palette for colours from RColorBrewer
#' @param tail.color a 2-components vector for the color of the tail, and the estimated density
#' @param alpha alpha channel
#'
#' @return if assigned, it return the list of ggplot2 objects.
#' @export
#'
#' @import sads
#' @import RColorBrewer
#' @import ggplot2
#' @import ggthemes
#'
#' @examples something..
plot.dbpmm = function(x, annotation = NULL, palette = 'Spectral', tail.color = c('gainsboro', 'darkgray'), alpha = .5, max.height.hist = TRUE, cex = 1,
                      plaw.ranges = c(0.1, 0.3), bg.color = 'ivory2', main = 'DBPMM fit', ...)
{
  domain = seq(0, 1, 0.01)
  labels = names(x$pi)
  labels.betas = colnames(x$beta)

  pi = x$pi
  ICL = round(as.numeric(x$scores$ICL), 2)
  NLL = round(as.numeric(x$scores$NLL), 2)
  K = as.numeric(x$K)

  col = RColorBrewer::brewer.pal(
    RColorBrewer::brewer.pal.info[palette, 'maxcolors'], palette)
  col = colorRampPalette(col)(x$Kbeta)

  col.histogram = c(tail.color[1], col)
  col.lines = c(tail.color[2], col)

  names(col.histogram) = names(col.lines) = labels



  # Plot 1 -- main histogram
  df = data.frame(X = x$X, Cluster = x$labels, Color = col[x$labels])
  vvv = lapply(1:x$K,
               function(w)
                 data.frame(X = domain,
                            Cluster = labels[w],
                            y = ddbpmm(x, data = domain, components = w, log = FALSE)))
  names(vvv) = labels

  # Histogram coloured according to clustering assignments
  # tit = bquote(bold(main) ~ K['Beta']~'=' ~ .(x$Kbeta) ~ '    Tail: OFF')
  # if(!all(is.na(x$tail)))
  #   tit = bquote(bold(main) ~
  #                       K['Beta']~'=' ~ .(x$Kbeta)  ~ '   '~
  #                       tau['tail'] ~'='~.(round(x$shape, 2)) ~ ' '~ gamma['tail'] ~'='~.(round(x$scale, 2)))
  tit = bquote(bold(.(main)))

  extented.labels = levels(df$Cluster)
  if(!all(is.na(x$tail))) extented.labels['Tail'] = paste('Tail : ', round(x$shape, 2), '; x >', round(x$scale, 2), sep ='')
  else extented.labels['Tail'] = bquote('Tail: OFF')

  betavals = x$beta[c('mean', 'var'), , drop = FALSE]
  betavals['mean', ] = round(betavals['mean', ], 3)
  betavals['var', ]  = format(betavals['var', ] , scientific = T, digits = 3)

  for(clus in paste('C', 1:x$Kbeta, sep = ''))
    extented.labels[clus] = paste(clus, ' : ', betavals['mean', clus], ' (', betavals['var', clus], ')', sep = '')


  p = ggplot(df, aes(X, fill = Cluster, y = ..count../sum(..count..))) +
    geom_histogram(alpha = alpha, position = 'identity', binwidth = 0.01) +
    scale_fill_manual(values = col.histogram, labels = extented.labels) +
    labs(
      title = tit,
      subtitle = annotation,
      x = "Observed Frequency", y = "") +
    theme_classic(base_size = 10 * cex)


  if(max.height.hist) {
    yMax = max(ggplot_build(p)$data[[1]]$y)
    p = p + ylim(0, yMax)
  }


  df.m = data.frame(variable = labels.betas, Mean = x$beta['mean', labels.betas])
  rownames(df.m) = df.m$variable

  # Add densities to the plot
  # scaling.ggplot = ggplot_build(p)$data[[1]]
  # scaling.denfun = max(unlist(vvv))
  # scaling = max(max(scaling.ggplot$y), scaling.denfun)

  # prop.to = max(scaling.ggplot$y)/scaling.denfun
  #
  # vvv = lapply(vvv, function(w) {w$y = w$y * prop.to; w})
  vvv = lapply(vvv, function(w) {w$y = w$y * 0.01; w}) # Scale density wrt binwidth

  for(i in seq(labels)) p = p + geom_line(data = vvv[[i]], aes(y = y, x = X), colour = col.lines[i])

  # Add means for Beta components
  for(i in labels.betas) p = p +  geom_vline(data = df.m[i ,], aes(xintercept = Mean), colour = col.lines[i], linetype = "longdash")

  # Annotate convergency value
  if(x$status) p = p + annotate("text", x = .9, y = yMax , label = paste(x$fit.type, ': CONVERGED'), size = 3, colour = 'darkgreen')
  else   p = p + annotate("text", x = .9, y = yMax, label = paste(x$fit.type, ': NOT CONVERGED'), size = 3, colour = 'red')


  tab.annotations = data.frame(round(x$beta[c('mean', 'var'), ], 2))
  colnames(tab.annotations) = colnames(x$beta)

  # p = p + annotation_custom(gridExtra::tableGrob(tab.annotations),  xmin = .9,  ymax= Inf)


  # Initial Condition -- density plot
  subdomain = domain[2:(length(domain) - 1)]
  vvvi = lapply(1:x$K,
               function(w)
                 data.frame(X = subdomain,
                            Cluster = labels[w],
                            y = ddbpmm(x, data = subdomain, components = w, init = TRUE, log = FALSE)))
  names(vvvi) = labels

  maxBeta = max(unlist(vvvi[2:x$K]))

  pinit = ggplot() +
    labs(
      title = bquote(italic("Initialization")),
      x = "", y = "") +
    guides(fill = FALSE) +
    ylim(0, maxBeta) +
    theme(
      panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

    # theme_pander(base_size = 7 * cex)


  for(i in seq(labels)) pinit = pinit + geom_line(data = vvvi[[i]], aes(y = y, x = X), colour = col.lines[i])

  df.mi = data.frame(variable = labels.betas,
                     Mean = .MeanVarBeta(x$beta['a.init', labels.betas], x$beta['b.init', labels.betas]))
  rownames(df.mi) = df.mi$variable

  for(i in labels.betas) pinit = pinit +  geom_vline(data = df.mi[i ,], aes(xintercept = Mean.mean), colour = col.lines[i], linetype = "longdash")

  # Mixing proportions
  df.pi = data.frame(Cluster = names(pi), Proportions = pi)

  h = ggplot(data = df.pi, aes(x = Cluster, y = Proportions, fill = Cluster)) +
    geom_bar(stat = "identity", alpha = alpha, width = 0.3 * cex) +
    scale_fill_manual(values = col.histogram) +
    coord_flip() +
    geom_hline(aes(yintercept = 0.02), colour = 'red', linetype = "longdash") +
    labs(title  = bquote(italic('Proportions'))) +
    xlab("") + ylab("") +
    guides(fill = FALSE) +
    theme(
      panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
    #
    # theme_pander(base_size = 7 * cex)

  ###### Mean parameter (Beta)
  df.pars = NULL
  for(i in seq(labels.betas))
    df.pars = rbind(
      df.pars,
      data.frame(Cluster = labels.betas[i],
                 Sample = rbeta(1000, x$beta['a', i], x$beta['b', i])))

  df.tail = NULL
  if(!all(is.na(x$tail))) {
    df.tail = data.frame(Cluster = labels[1],
                       Sample = sads::rpareto(1000, shape = x$tail$shape, scale = x$tail$scale))
    df.tail = df.tail[df.tail$Sample < 1, ]
  }

  pa = ggplot(data = rbind(df.tail, df.pars), aes(x = Cluster, y = Sample, fill = Cluster)) +
    geom_boxplot(alpha = alpha) +
    scale_fill_manual(values = col.histogram) +
    guides(fill=FALSE) +
    labs(title  = bquote(italic('Means'))) +
    xlab('') + ylab("") +
    theme(
      panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

  # +
  #   theme_pander(base_size = 7 * cex)


  ###### Scores
  df.scores =  data.frame(value = t(x$scores))
  df.scores$score = rownames(df.scores)

  ps = ggplot(data = df.scores, aes(x = score, y = value)) +
    # geom_boxplot(alpha = alpha) +
    geom_point(shape = 21, colour = "black", fill = "red", size = 1 * cex, stroke = 1) +
    # scale_colour_gradient(low = "blue")+
    guides(fill = FALSE) +
    labs(title  = bquote(italic('Scores'))) +
    xlab('') + ylab('') +
    theme(
      panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    coord_flip()

    # +
    # theme_pander(base_size = 7 * cex)



  # GET EQUATION AND R-SQUARED AS STRING
  # SOURCE: http://goo.gl/K4yh

  lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                     list(a = format(coef(m)[1], digits = 2),
                          b = format(coef(m)[2], digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
    format(summary(m)$r.squared, digits = 2)
  }

  if(all(!is.na(x$shape)))
  {
    x.tail = seq(plaw.ranges[1], plaw.ranges[2], 0.01)

    ## Cumdensity for the tail
    tail.X = x$X[x$labels == 'Tail']
    tail.X = tail.X[tail.X >= plaw.ranges[1] & tail.X <= plaw.ranges[2]]
    ht = hist(tail.X, breaks = x.tail, plot = FALSE)
    cms = cumsum(ht$counts)
    # cms = cms[1:which.max(cms)]
    # cms = cms[cms > 0]


    # cumdf = sort(1/tail.X)
    cumdf = data.frame(x = x.tail[1:length(cms)], y = cms)
    rs = lm_eqn(cumdf)

    r2 = ggplot(data = cumdf, aes(x = x, y = y)) +
      geom_point(shape = 21, colour = "black", size = 1 * cex, stroke = 1) +
      geom_smooth(method='lm', colour = "red") +
      guides(fill=FALSE) +
      labs(title  = bquote(italic('Power-law'))) +
      xlab('f ~ VAF') + ylab('Cumulative M(f)') +
      # scale_x_continuous(labels = NULL) +
      theme(
        panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
      annotate("text", x = Inf, y = min(cumdf$y) + 3, label = paste("R^2 == ", rs), hjust = 1, parse = T, size = 2 * cex)
  }
  else
  {
    r2 = ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100)
  }

  top.layout = matrix(1, nrow = 2, ncol = 5)
  bottom.layout = 2:6

  .multiplot(p, pinit, pa, h, r2, ps,
            layout = rbind(top.layout, bottom.layout))

  invisible(list(p, pa, h))
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
.multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

.plot.fit.summary = function(x, alpha = .5) {

  x = x[complete.cases(x), ]

  print(x)

  x$tail[x$tail] = 'With Tail'
  x$tail[x$tail == 'FALSE'] = 'Without Tail'

  x$K = as.factor(x$K)
  x$tail = as.factor(x$tail)

  outL = function(x) {
    # compute lower and upper whiskers
    boxplot.stats(x$ICL)$stats[5]
  }

  outliers = sapply(split(x, f = x$K), outL)
  outliers = max(outliers)

  range = c(min(x$ICL), outliers)

  pa = ggplot(data = x, aes(x = K, y = ICL, fill = K)) +
    geom_boxplot(alpha = alpha) +
    scale_fill_brewer(palette = 'Set1') +
    labs(title  = bquote(bold('ICL Score: ')~ .(nrow(x)) ~' runs' ~ '- no upper outliers')) +
    xlab('K') + ylab(bquote(italic('ICL'))) +
    theme_light() +
    facet_wrap(~tail, nrow = 1) +
    coord_cartesian(ylim = range * 1.05)

  pa
}

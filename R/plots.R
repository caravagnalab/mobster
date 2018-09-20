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
plot.dbpmm = function(x, annotation = NULL, palette = 'Spectral', tail.color = c('darkgray', 'dimgrey'), alpha = .8, max.height.hist = TRUE, cex = 1,
                      plaw.ranges = c(0.1, 0.3), bg.color = 'ivory2', histogram.main = 'Fit', title = 'MOBSTER',
                      plot.rsquare = FALSE, plot.scores = FALSE,
                      silent = FALSE, ...)
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
  col = scols(1:x$Kbeta, palette = palette)

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
  tit = bquote(bold(.(histogram.main)))

  extented.labels = levels(df$Cluster)
  if(!all(is.na(x$tail))) extented.labels['Tail'] = paste('Tail : ', round(x$shape, 2), '; x >', round(x$scale, 2), sep ='')
  else extented.labels['Tail'] = bquote('Tail: OFF')

  betavals = x$beta[c('mean', 'var'), , drop = FALSE]
  betavals['mean', ] = round(betavals['mean', ], 3)
  betavals['var', ]  = format(betavals['var', ] , scientific = T, digits = 3)

  for(clus in paste('C', 1:x$Kbeta, sep = ''))
    extented.labels[clus] = paste(clus, ' : ', betavals['mean', clus], ' (', betavals['var', clus], ')', sep = '')


  conv.steps = length(x$all.NLL)
  conv.epsilon = 0
  if(conv.steps >= 2) conv.epsilon = abs(x$all.NLL[length(x$all.NLL)] - x$all.NLL[length(x$all.NLL) - 1])

  label.fit = ifelse(x$fit.type == 'MM', 'Moments Matching', "Maximum Likelihood")
  label.fit = paste0(label.fit)

  p = ggplot(df, aes(X, fill = Cluster, y = ..count../sum(..count..))) +
    geom_histogram(alpha = alpha, position = 'identity', binwidth = 0.01) +
    scale_fill_manual(values = col.histogram, labels = extented.labels) +
    labs(
      title = tit,
      subtitle = annotation,
      caption = bquote(.(label.fit) ~ '-' ~ .(conv.steps) ~ 'steps,' ~ epsilon ~ '=' ~ .(conv.epsilon)),
      x = "Observed Frequency", y = "Density") +
    theme_classic(base_size = 10 * cex) +
    geom_vline(xintercept = min(x$X), colour = 'black', linetype = "longdash") +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      plot.caption = element_text(color = ifelse(x$status, "darkgreen",  "red"))
    )



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

  # # Annotate convergency value
  # if(x$status) p = p + annotate("text", x = .9, y = yMax , label = paste(x$fit.type, ': CONVERGED'), size = 3, colour = 'darkgreen')
  # else   p = p + annotate("text", x = .9, y = yMax, label = paste(x$fit.type, ': NOT CONVERGED'), size = 3, colour = 'red')


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
    theme_light(base_size = 8 * cex) +
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
    theme_light(base_size = 8 * cex) +
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
    geom_violin(alpha = alpha, trim = TRUE) +
    geom_boxplot(width = .1) +
    scale_fill_manual(values = col.histogram) +
    guides(fill=FALSE) +
    labs(title  = bquote(italic('Means'))) +
    xlab('') + ylab("") +
    theme_light(base_size = 8 * cex) +
    theme(
      panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

  ###### Scores
  ps = NULL
  if(plot.scores)
  {
    df.scores =  data.frame(value = t(x$scores))
    df.scores$score = rownames(df.scores)

    ps = ggplot(data = df.scores, aes(x = score, y = value)) +
      # geom_boxplot(alpha = alpha) +
      geom_point(shape = 21, colour = "black", fill = "red", size = 1 * cex, stroke = 1) +
      # scale_colour_gradient(low = "blue")+
      guides(fill = FALSE) +
      labs(title  = bquote(italic('Scores'))) +
      xlab('') + ylab('') +
      theme_light(base_size = 8 * cex) +
      theme(
        panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
      coord_flip()
  }



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

  r2 = NULL
  if(plot.rsquare)
  {
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
        theme_light(base_size = 8 * cex) +
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
  }

  ncols = 3
  if(plot.rsquare & plot.scores) ncols = ncols + 2
  if(!plot.rsquare & plot.scores) ncols = ncols + 1
  if(plot.rsquare & !plot.scores) ncols = ncols + 1

  top.layout = matrix(1, nrow = 2, ncol = ncols)
  bottom.layout = 2:(ncols + 1)

  cat("Layout: ")

  if(!silent & plot.rsquare & plot.scores)
    .multiplot(p, pinit, pa, h, r2, ps,
               layout = rbind(top.layout, bottom.layout), title = bquote(bold(.(title))), fontsize = 16 * cex)

  if(!silent & !plot.rsquare & plot.scores)
    .multiplot(p, pinit, pa, h, ps,
               layout = rbind(top.layout, bottom.layout), title = bquote(bold(.(title))), fontsize = 16 * cex)

  if(!silent & plot.rsquare & !plot.scores)
    .multiplot(p, pinit, pa, h, ps,
               layout = rbind(top.layout, bottom.layout), title = bquote(bold(.(title))), fontsize = 16 * cex)

  if(!silent & !plot.rsquare & !plot.scores)
    .multiplot(p, pinit, pa, h,
               layout = rbind(top.layout, bottom.layout), title = bquote(bold(.(title))), fontsize = 16 * cex)

  invisible(list(mainHist = p, Parameters = pa, Proportions = h, Initialization = pinit))
}


scols = function (v, palette = "Spectral")
{
  colors = NULL

  pmax = RColorBrewer::brewer.pal.info[palette, "maxcolors"]

  if (length(v) <= pmax)
  {
    colors = RColorBrewer::brewer.pal(n = length(v), palette)
    colors = colors[1:length(v)]
  }
  else
  {
    colors = RColorBrewer::brewer.pal(n = pmax, palette)
    colors = colorRampPalette(colors)(length(v))
  }

  names(colors) = v

  return(colors)
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
.multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL, title="MOBSTER",
                       fontsize = 18, fontfamily = "Helvetica") {
  require(grid)

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

  if (nchar(title)>0){
    layout <- rbind(rep(0, ncol(layout)), layout)
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout),
                                               ncol(layout),
                                               heights = if (nchar(title)>0) {unit(c(0.5, rep(5,nrow(layout)-1)), "null")}
                                               else {unit(c(rep(5, nrow(layout))), "null")})))

    # Make each plot, in the correct location
    if (nchar(title)>0) {
      grid.text(title,
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol(layout)),
                gp = gpar(fontsize = fontsize, fontfamily = fontfamily))
    }

    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

.plot.fit.summary = function(x, alpha = .5, silent = FALSE, range = NULL, cex = 1) {

  model.selection = 'ICL'
  if(!is.null(x$model.selection)) model.selection = x$model.selection

  x = x$fits.table

  pio::pioTit("Creating boxplot of scores for model selection.")
  x = x[complete.cases(x), ]

  # It's in the table
  model.selection = x$model.selection[1]

  pio::pioDisp(x)

  x$tail[x$tail] = 'With Tail'
  x$tail[x$tail == 'FALSE'] = 'Without Tail'

  x$K = as.factor(x$K)
  x$tail = as.factor(x$tail)
  values = x[, model.selection]

  # outL = function(x) {
  #   # compute lower and upper whiskers
  #   boxplot.stats(values)$stats[5]
  # }
  #
  # if(all(is.null(range)))
  # {
  #   message('range = NULL : removing score outliers')
  #   outliers = sapply(split(x, f = x$K), outL)
  #   outliers = max(outliers)
  #
  #   range = c(min(values), outliers) * 1.05
  # }
  # else message("Using custom input range for ICL")


  # pa = ggplot(data = x, aes(x = K, y = ICL, fill = K)) +
  #   geom_jitter(position = position_jitter(width = 0.2, height = 0.005),
  #               alpha = 0.05, aes(color = K),
  #               size = 2) +
  #   geom_violin(alpha = alpha, trim = FALSE) +
  #   # geom_boxplot(alpha = alpha) +
  #   scale_fill_brewer(palette = 'Set1') +
  #   labs(title  = bquote(bold('ICL Score: ')~ .(nrow(x)) ~' runs' ~ '- no upper outliers')) +
  #   xlab('K') + ylab(bquote(italic('ICL'))) +
  #   theme_light() +
  #   facet_wrap(~tail, nrow = 1) +
  #   coord_cartesian(ylim = range * 1.05)

  pa = ggplot(data = x, aes(x = K, y = ICL, fill = K)) +
    geom_jitter(position = position_jitter(width = 0.2, height = 0.005),
                alpha = alpha, aes(colour = K),
                size = 2) +
    labs(title  = bquote(bold('Model selection')),
         subtitle = bquote(.(nrow(x)) ~' runs')) +
    xlab('Beta components (K)') + ylab(bquote(italic(.(model.selection)))) +
    theme_light(base_size =  8 * cex) +
    facet_wrap(~tail, nrow = 1) +
    # coord_cartesian(ylim = range) +
    guides(fill = FALSE, colour = FALSE)


  # require(ggrepel)

  if(!silent) print(pa)

  invisible(list(boxPlot = pa))
}


plot.fits = function(fits)
{
  .getGGplotHistogram(fits[[1]])


}


.getGGplotHistogram = function(x, annotation = NULL, palette = 'Spectral',
                               tail.color = c('gainsboro', 'darkgray'),
                               alpha = .5, max.height.hist = TRUE, cex = 1,
                               bg.color = 'ivory2', main = 'MOBSTER fit', ...)
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
  vvv = lapply(vvv, function(w) {w$y = w$y * 0.01; w}) # Scale density wrt binwidth

  for(i in seq(labels)) p = p + geom_line(data = vvv[[i]], aes(y = y, x = X), colour = col.lines[i])

  # Add means for Beta components
  for(i in labels.betas) p = p +  geom_vline(data = df.m[i ,], aes(xintercept = Mean), colour = col.lines[i], linetype = "longdash")

  # Annotate convergency value
  if(x$status) p = p + annotate("text", x = .9, y = yMax , label = paste(x$fit.type, ': CONVERGED'), size = 3, colour = 'darkgreen')
  else   p = p + annotate("text", x = .9, y = yMax, label = paste(x$fit.type, ': NOT CONVERGED'), size = 3, colour = 'red')

  p
}



.getIconGGplotHistogram = function(x, palette = 'Spectral',
                                   tail.color = c('gainsboro', 'darkgray'),
                                   alpha = .5, max.height.hist = TRUE, cex = 1,
                                   bg.color = 'ivory2')
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

  p = ggplot(df, aes(X, fill = Cluster, y = ..count../sum(..count..))) +
    geom_histogram(alpha = alpha, position = 'identity', binwidth = 0.01) +
    scale_fill_manual(values = col.histogram, labels = labels) +
    theme_void(base_size = 10 * cex) +
    guides(fill = FALSE)


  if(max.height.hist) {
    yMax = max(ggplot_build(p)$data[[1]]$y)
    p = p + ylim(0, yMax)
  }


  df.m = data.frame(variable = labels.betas, Mean = x$beta['mean', labels.betas])
  rownames(df.m) = df.m$variable

  # Add densities to the plot
  vvv = lapply(vvv, function(w) {w$y = w$y * 0.01; w}) # Scale density wrt binwidth

  for(i in seq(labels)) p = p + geom_line(data = vvv[[i]], aes(y = y, x = X), colour = col.lines[i])

  p
}


# mobster_oncoplot = function(data, samples, col.clonal = 'steelblue', col.subclonal = 'orange', col.private = 'red', col.below.cutoff = 'darkgray', col.background = 'gainsboro', cutoff = 0.05)
# {
#   scoreCol = function(x) {
#     score = 0
#     for(i in 1:length(x)) {
#       if(x[i]) {
#         score = score + 2^(length(x)-i*1/x[i])
#       }
#     }
#     return(score)
#   }
#
#   # Values
#   inp = data[, samples]
#
#   # NA for everything below the cutoff
#   inp[inp < cutoff] = NA
#
#   # Clonal, subclonal and private values
#   which.clonal = apply(inp, 1,
#                        function(x) all(!is.na(x) & x > cutoff))
#
#   which.subclonal = apply(inp, 1,
#                           function(x) {
#                             k = length(x)
#                             x = x[!is.na(x)]
#                             sum(as.numeric(x > cutoff)) < k & sum(as.numeric(x > cutoff)) > 1
#                             })
#
#   which.private = apply(inp, 1,
#                         function(x) {
#                           x = x[!is.na(x)]
#                           sum(as.numeric(x > cutoff)) == 1
#                           })
#
#   # Clonal sorted by frequency
#   clonal = inp[which.clonal, ]
#   clonal[T] = 1
#   # clonal = clonal[order(rowSums(clonal, na.rm = T), decreasing = T), ]
#
#   head(clonal)
#
#   # Subclonal sorted by mutual exclusivity
#   subclonal = inp[which.subclonal, ]
#   subclonal[!is.na(subclonal)] = 1
#   subclonal[is.na(subclonal)] = 0
#
#   scores   = apply(subclonal, 1, scoreCol)
#   subclonal = subclonal[order(scores, decreasing=TRUE), ]
#
#
#   head(subclonal)
#
#
#
#
#   # Pheatmap
#   plot.data = rbind(clonal, subclonal)
#   head(plot.data)
#   tail(plot.data)
#
#   colors = c(
#     `0` = col.background,
#     `1` = col.clonal,
#     `2` = col.subclonal,
#     `3` = col.private)
#
#
#   require(pheatmap)
#   pheatmap(
#     plot.data,
#     cluster_cols = F,
#     cluster_rows = F,
#     color = colors,
#     border_color = NA,
#     na_col = col.below.cutoff)
#
#
#   Order = function(data) {
#     data[is.na(data)] = 0
#
#     scoreCol = function(x) {
#       score = 0
#       for(i in 1:length(x)) {
#         if(x[i]) {
#           score = score + 2^(length(x)-i*1/x[i])
#         }
#       }
#       return(score)
#     }
#
#     sharedData  = data[which(apply(data, 1, sum) > 1),]
#     sharedIndex = which(apply(data, 1, sum) > 1)
#
#     scores   = apply(sharedData, 1, scoreCol)
#     topOrder = sharedIndex[order(scores, decreasing=TRUE)]
#
#     privateData  = data[which(apply(data, 1, sum) <= 1),]
#     privateIndex = which(apply(data, 1, sum) <= 1)
#
#     scores   = apply(privateData, 1, scoreCol)
#     bottomOrder = privateIndex[order(scores, decreasing=TRUE)]
#
#     return(c(topOrder, bottomOrder))
#   }
#
#
#
#
#
#
#
#
# }


#' Title
#'
#' @param res
#' @param fig.lab
#' @param title
#' @param palette
#' @param TOP
#' @param cex
#' @param boxplot.ICL.range
#'
#' @return
#' @export
#'
#' @examples
plot_report_MOBSTER = function(res, fig.lab = "Figure 1", title = 'MOBSTER top fit', palette = 'Set1', TOP = 5, cex = 1, boxplot.ICL.range = NULL)
{
  # require(dbpmm)
  require(ggplot2)
  require(ggpubr)

  # Best model fit -- get plots
  best = plot.dbpmm(
    res$best,
    plot.rsquare = F,
    plot.scores = F,
    alpha = .9,
    tail.color = c('darkgray', 'dimgrey'),
    bg.color = 'gainsboro',
    cex = 1 * cex,
    title = "",
    histogram.main = title,
    palette = palette,
    silent = TRUE,
    annotation = paste0('Top-fit from N = ', length(res$best$X))
  )

  # Extract scores table
  scores = lapply(res$runs,
                  function(x)
                  {
                    r = x$scores
                    r$K = x$Kbeta
                    r$Tail = ifelse(all(is.na(x$tail)), 'NO', 'YES')

                    r
                  }
  )
  scores = Reduce(rbind, scores)

  # Set format for columns with numbers
  scores.col = colnames(res$best$scores)
  scores[, scores.col] = apply(scores[, scores.col], 2, round, digit = 0)

  # Subset TOP scores and create a dataframe for plot
  scores = scores[!duplicated(scores[, c('K', 'Tail')]), ]
  if(nrow(scores) > TOP) scores = scores[1:TOP, ]
  pio::pioTit("Table with TOP scores")
  pio::pioDisp(scores)

  library(ggpubr)
  df.scores = ggtexttable(scores, rows = NULL,
                          theme = ttheme(base_size = 6 * cex,
                                         tbody.style = tbody_style(
                                           fill = rev(get_palette("RdBu", TOP)),
                                           size = 6 * cex)))

  for(s in scores.col){
    idx = which.min(scores[, s]) + 1
    df.scores = table_cell_bg(df.scores, row = idx, column = which(colnames(scores) == s),
                              linewidth = 2.5,
                              fill=alpha("darkgreen", .5))
  }

  # Model selection boxplot
  boxP = .plot.fit.summary(
    res, alpha = .6, silent = TRUE, range = boxplot.ICL.range)$boxPlot

  # Other solutions ranked below top best -- maximum TOP 4 fixed
  TOP.plot = min(4, nrow(scores))

  solutions = as.integer(rownames(scores)[2:(TOP.plot + 1)])
  solutions = solutions[!is.na(solutions)]

  other.best = lapply(
    seq(solutions),
    function(w)
      plot(
        res$runs[[solutions[w]]],
        histogram.main = paste0("Solution #", 1+w),
        annotation = paste0("Overall rank : ", solutions[w]),
        silent = T,
        palette = palette,
        cex = .6 * cex)$mainHist
  )

  # Final layout
  top_panel_best = ggarrange(
    best$mainHist,
    ggarrange(
      best$Initialization,
      best$Parameters,
      best$Proportions,
      ncol = 3,
      nrow = 1,
      labels = c("B", "C", "D")
    ),
    ggarrange(
      boxP,
      df.scores,
      ncol = 2,
      nrow = 1,
      labels = c("E", "F"),
      widths = c(1, 1.5)
    ),
    heights = c(2.5, 1),
    nrow = 3,
    ncol = 1,
    common.legend = T,
    labels = 'A'
  )

  right_panel = ggarrange(
    plotlist = other.best,
    ncol = 1,
    nrow = length(other.best),
    labels = LETTERS[7:(7 + TOP.plot - 1)]
  )

  figure = ggarrange(
    top_panel_best,
    right_panel,
    widths = c(1.5,1),
    nrow = 1, ncol = 2
  )

  figure = annotate_figure(
    figure,
    top = text_grob("MOBSTER (model selection)", color = "black", face = "bold", size = 22 * cex),
    bottom = text_grob(
      " Panels: (A-D) Best fit; (E,F) Scores from model selection; (G-*) Lower-scoring fits.",
      color = "black",
      face = "bold",
      hjust = 0,
      x = 0,
      size = 8 * cex),
    fig.lab = fig.lab,
    fig.lab.face = "bold"
  )

  return(figure)
}

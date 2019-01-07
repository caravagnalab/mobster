#' Plot an object of class dbpmm. This is a set of multiple frames plot via ggplot2.
#'
#' @param x an object of class dbpmm
#' @param palette a palette for colours from RColorBrewer
#' @param tail.color a 2-components vector for the color of the tail, and the estimated density
#' @param alpha alpha channel
#' @param annotation 
#' @param binwidth 
#' @param max.height.hist 
#' @param cex 
#' @param histogram.main 
#' @param bg.color 
#' @param title 
#' @param silent 
#' @param ... 
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
plot.dbpmm = function(x, 
                      annotation = NULL, 
                      palette = 'Set1', 
                      tail.color = 'darkgray', 
                      alpha = .8, 
                      binwidth = 0.01,
                      max.height.hist = TRUE, 
                      annotate = NULL,
                      cex = 1,
                      histogram.main = 'Fit', # Title of the histogram
                      bg.color = 'ivory2',    # Background color
                      title = 'MOBSTER',      # Plot title
                      silent = FALSE,         # Just return plots
                      ...)
{
  # Prepare variables
  domain = seq(0, 1, binwidth)
  
  labels = names(mobster:::.params_Pi(x))
  
  labels.betas = mobster:::.params_Beta(x)$cluster

  pi = mobster:::.params_Pi(x)
  
  ICL = round(as.numeric(x$scores$ICL), 2)
  NLL = round(as.numeric(x$scores$NLL), 2)
  K = as.numeric(x$K)

  # Load colors
  colors = mobster:::getColors_model(x, alpha = .9, palette = palette, tail.color = tail.color)

  ############### Plot 1 -- main histogram
  
  # Main plotting data
  x$data$cluster = factor(x$data$cluster, levels = names(colors))
  
  # Text for the plot -- convergence
  conv.steps = length(x$all.NLL)
  conv.epsilon = 0
  if(conv.steps >= 2) conv.epsilon = abs(rev(x$all.NLL)[1] - rev(x$all.NLL)[2])
  conv.epsilon = formatC(conv.epsilon, format = "e", digits = 0)
  
  sse = max(mobster:::.compute_fit_sqerr(x, binning = binwidth)$cum.y)
  sse = formatC(sse, format = "e", digits = 3)

  label.fit = bquote(
    .(x$fit.type) *
      " ("*omega*" = "* .(conv.steps) ~ 'steps; ' * epsilon ~ '=' ~ .(conv.epsilon) *
      "; SSE"~.(sse)*')')
  
  # label.fit = label.fit ~ ifelse(x$fit.tail, 'With tail', "Without tail")

  
  # Main ggplot object
  hist_pl = ggplot(x$data, aes(VAF, fill = factor(cluster), y = ..count../sum(..count..))) +
    geom_histogram(alpha = alpha, position = 'identity', binwidth = binwidth) +
    scale_fill_manual(values = colors, labels = names(colors)) +
    guides(fill = guide_legend(title = "Cluster")) + 
    labs(
      title = bquote(bold(.(histogram.main))),
      subtitle = annotation,
      caption = label.fit,
      x = "Observed Frequency", y = "Density") +
    theme_classic(base_size = 10 * cex) +
    geom_vline(xintercept = min(x$data$VAF), colour = 'black', linetype = "longdash") +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      panel.background = element_rect(fill = 'white'),
      plot.caption = element_text(color = ifelse(x$status, "darkgreen",  "red"))
    )
  
  #### Prepare plot for density values
  densities = mobster:::template_density(
    x, 
    x.axis = domain, 
    binwidth = binwidth, 
    reduce = TRUE)
    
  # Add the trace and the mean of each component
  hist_pl = hist_pl + 
    geom_line(data = densities, aes(y = y, x = x, color = factor(cluster)), size = 1 * cex) +
    scale_color_manual(values = colors, labels = names(colors)) +
    guides(color = FALSE)
    
  Beta_peaks = x$Clusters %>%
    dplyr::filter(type == 'Mean', cluster != 'Tail')
  
  hist_pl = hist_pl + 
    geom_vline(data = Beta_peaks, aes(xintercept = fit.value, color = cluster), linetype = "longdash")
  
  # Rescale y in case the density goes too high
  if(max.height.hist) 
  {
    yMax = max(ggplot_build(hist_pl)$data[[1]]$y)
    hist_pl = hist_pl + ylim(0, yMax)
  }
  
  #### Overall density and SSE
  densities = tibble::as_tibble(densities)
  densities = densities %>% group_by(x) %>% summarise(y = sum(y), cluster = 'f(x)')
  
  m = max(densities$y, na.rm = TRUE)

  hist_pl = hist_pl +
    geom_line(data = densities %>% mutate(y = y + m * 0.02), 
              aes(y = y, x = x), 
              color = 'black',
              alpha = .8,
              size = .5 * cex,
              linetype = 'dashed',
              inherit.aes = FALSE)
  
  # SSE
  error = .compute_fit_sqerr(x, binning = binwidth)

  # scale to percentage for plotting
  me = max(error$cum.y)
  error = error %>% mutate(cum.y = (cum.y / me) * m )
  
  hist_pl = hist_pl +
    geom_line(data = error, aes(y = cum.y, x = x), 
              color = 'darkgray',
              alpha = 1,
              size = .2 * cex,
              linetype = 'dashed',
              inherit.aes = FALSE) +
    scale_y_continuous(sec.axis = sec_axis(~./ m * 100, name = "SSE [cumulative %]"))  

  # Annotation of input entries
  if(!is.null(annotate) & all(c("VAF", 'label') %in% colnames(annotate)))
  {
    # Position the point at coord x = VAF and y = density
    m = max(densities$y, na.rm = TRUE)

    annotate$y = round(annotate$VAF/binwidth)
    annotate$y = densities$y[annotate$y] + m * 0.02
    
    hist_pl = hist_pl +
      geom_label_repel(data = annotate, 
                       aes(
                         x = VAF, 
                         y = y,
                         label = label, 
                         color = factor(cluster)
                         # fill = factor(cluster)
                       ),
                       size = 1.5 * cex,
                       inherit.aes = FALSE,
                       box.padding = 0.95, 
                       segment.size = .2 * cex, force = 1) +
      geom_point(data = annotate, aes(x = VAF, y = y,  color = factor(cluster)), 
                 size = .3 * cex, alpha = 1, 
                 inherit.aes = FALSE)
  }

  ############### Plot 2 -- initial Condition, density plot
  n = x
  n$Clusters$fit.value = n$Clusters$init.value
  
  initial.densities = mobster:::template_density(
    n, 
    x.axis = domain[2:(length(domain) - 1)], # Restricted for numerical errors
    binwidth = 0.01,
    reduce = TRUE)
  
  den_init_pl = ggplot() +
    labs(
      title = bquote("Initialization"),
      x = "", y = "") +
    guides(fill = FALSE) +
    ylim(0, max(initial.densities$y)) +
    theme_light(base_size = 8 * cex) +
    theme(
      panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
    
  den_init_pl = den_init_pl + 
    geom_line(data = initial.densities, aes(y = y, x = x, color = cluster)) +
    scale_color_manual(values = colors, labels = names(colors)) +
    guides(color = FALSE)
  
  ############### Plot 2 -- Mixing proportions
  Proportions = x$Clusters %>%
    dplyr::filter(type == 'Mixing proportion')
  
  Proportions$fit.value = round(Proportions$fit.value, 2)
  
  if(!x$fit.tail) Proportions = Proportions %>% filter(cluster != 'Tail')
  
  prop_hist_pl = ggplot(data = Proportions, aes(x = cluster, y = fit.value, fill = cluster)) +
    geom_bar(stat = "identity", alpha = alpha, width = 0.3 * cex) +
    scale_fill_manual(values = colors, labels = names(colors)) +
    coord_flip() +
    geom_hline(aes(yintercept = 0.02), colour = 'red', linetype = "longdash") +
    labs(title  = bquote('Proportions')) +
    xlab("") + ylab("") +
    guides(fill = FALSE) +
    theme_light(base_size = 8 * cex) +
    theme(
      panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    ylim(c(0,1))
    
  

  ###### Mean parameter (Beta)
  Beta.a = x$Clusters %>%
    dplyr::filter(cluster != 'Tail', type == 'a') %>%
    dplyr::select(fit.value) %>%
    dplyr::pull()
  
  Beta.b = x$Clusters %>%
    dplyr::filter(cluster != 'Tail', type == 'b') %>%
    dplyr::select(fit.value) %>%
    dplyr::pull()
  
  df.pars = NULL
  for(i in seq(labels.betas))
    df.pars = rbind(
      df.pars,
      data.frame(Cluster = labels.betas[i],
                 Sample = rbeta(1000, Beta.a[i], Beta.b)))

  ###### Mean parameter (Pareto)
  shape = x$Clusters %>%
    dplyr::filter(cluster == 'Tail', type == 'Shape') %>%
    dplyr::pull(fit.value)
  
  scale = x$Clusters %>%
    dplyr::filter(cluster == 'Tail', type == 'Scale') %>%
    dplyr::pull(fit.value)
  
  
  df.tail = NULL
  if(x$fit.tail) 
  {
    df.tail = data.frame(Cluster = labels[1],
                         Sample = sads::rpareto(1000, shape = shape, scale = scale))
    
    df.tail = df.tail[df.tail$Sample < 1, ]
  }

  box_mean_pl = ggplot(data = rbind(df.tail, df.pars), aes(x = Cluster, y = Sample, fill = Cluster)) +
    geom_violin(aes(color = NULL), alpha = alpha, trim = TRUE) +
    geom_boxplot(width = .1, outlier.size = .5) +
    scale_fill_manual(values = colors, labels = names(colors)) +
    guides(fill=FALSE) +
    labs(title  = bquote('Means')) +
    xlab('') + ylab("") +
    theme_light(base_size = 8 * cex) +
    theme(
      panel.background = element_rect(fill = alpha(bg.color, 0.4), colour = "white", size = 0.5, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())


  ######################### Figure layout 
  figure = ggpubr::ggarrange(
    hist_pl,
    ggpubr::ggarrange(
      den_init_pl, prop_hist_pl, box_mean_pl,
      ncol = 3, nrow  = 1
    ),
    nrow = 2,
    ncol = 1, 
    heights = c(2,1)
  )
  
  figure = ggpubr::annotate_figure(figure, top = title)


  if(!silent) print(figure)

  invisible(list(mainHist = hist_pl, 
                 Parameters = box_mean_pl, Proportions = prop_hist_pl, Initialization = den_init_pl))
}


scols = function (v, palette = "Spectral")
{
  colors = NULL

  pmax = RColorBrewer::brewer.pal.info[palette, "maxcolors"]

  if (length(v) <= pmax)
  {
    colors = suppressWarnings(RColorBrewer::brewer.pal(n = length(v), palette))
    colors = colors[1:length(v)]
  }
  else
  {
    colors = suppressWarnings(RColorBrewer::brewer.pal(n = pmax, palette))
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


.plot.fit.summary = function(x, TOP = 5, alpha = .5, silent = FALSE, range = NULL, cex = 1) {

  model.selection = 'ICL'
   if(!is.null(x$model.selection)) model.selection = x$model.selection
  # 
  model.selection = x$model.selection[1]
  x = x$fits.table
  
  pio::pioTit("Creating boxplot of scores for model selection.")
  x = x[complete.cases(x), ]
  
  # It's in the table
  pio::pioDisp(x)
  
  x$tail[x$tail] = 'With Tail'
  x$tail[x$tail == 'FALSE'] = 'Without Tail'
  
  x$K = as.factor(x$K)
  x$tail = as.factor(x$tail)
  values = x[, model.selection]
  
  x$rank = 1:nrow(x)
  
  pos <- position_jitter(width = 0.2, seed = 1)
  
  pa = ggplot(data = x, aes(x = K, y = eval(paste(text = model.selection)), fill = interaction(K,tail), colour = interaction(K,tail))) +
    geom_point(alpha = alpha, position = pos, size = 2) + 
    # geom_jitter(position = position_jitter(width = 0.2, height = 0.005),
    #             alpha = alpha, aes(colour = K),
    #             size = 2) +
    labs(title  = bquote(bold('Model selection')),
         subtitle = bquote(.(nrow(x)) ~' runs')) +
    xlab('Beta components (K)') + ylab(bquote(italic(.(model.selection)))) +
    theme_light(base_size =  8 * cex) +
    facet_grid(K~tail) +
    guides(fill = FALSE, colour = FALSE)
  
  
  
  library(ggrepel)
  pa = pa + geom_label_repel(data = x[1:TOP, ], aes(label = rank), 
                        position = pos,
                        fontface = 'bold', 
                        color = 'black', 
                        angle = 45,
                        box.padding = 0.95, 
                        size = 0.5 * cex,
                        segment.size = 0.2 * cex,
                        segment.color = 'gainsboro')
  

  if(!silent) print(pa)

  invisible(list(boxPlot = pa))
}




# plot.fits = function(fits)
# {
#   .getGGplotHistogram(fits[[1]])
# 
# 
# }


# .getGGplotHistogram = function(x, annotation = NULL, palette = 'Spectral',
#                                tail.color = c('gainsboro', 'darkgray'),
#                                alpha = .5, max.height.hist = TRUE, cex = 1,
#                                bg.color = 'ivory2', main = 'MOBSTER fit', ...)
# {
#   domain = seq(0, 1, 0.01)
# 
#   labels = names(x$pi)
#   labels.betas = colnames(x$beta)
# 
#   pi = x$pi
#   ICL = round(as.numeric(x$scores$ICL), 2)
#   NLL = round(as.numeric(x$scores$NLL), 2)
#   K = as.numeric(x$K)
# 
#   col = RColorBrewer::brewer.pal(
#     RColorBrewer::brewer.pal.info[palette, 'maxcolors'], palette)
#   col = colorRampPalette(col)(x$Kbeta)
# 
#   col.histogram = c(tail.color[1], col)
#   col.lines = c(tail.color[2], col)
# 
#   names(col.histogram) = names(col.lines) = labels
# 
# 
#   # Plot 1 -- main histogram
#   df = data.frame(X = x$X, Cluster = x$labels, Color = col[x$labels])
#   vvv = lapply(1:x$K,
#                function(w)
#                  data.frame(X = domain,
#                             Cluster = labels[w],
#                             y = ddbpmm(x, data = domain, components = w, log = FALSE)))
#   names(vvv) = labels
# 
#   # Histogram coloured according to clustering assignments
#   tit = bquote(bold(.(main)))
# 
#   extented.labels = levels(df$Cluster)
#   if(!all(is.na(x$tail))) extented.labels['Tail'] = paste('Tail : ', round(x$shape, 2), '; x >', round(x$scale, 2), sep ='')
#   else extented.labels['Tail'] = bquote('Tail: OFF')
# 
#   betavals = x$beta[c('mean', 'var'), , drop = FALSE]
#   betavals['mean', ] = round(betavals['mean', ], 3)
#   betavals['var', ]  = format(betavals['var', ] , scientific = T, digits = 3)
# 
#   for(clus in paste('C', 1:x$Kbeta, sep = ''))
#     extented.labels[clus] = paste(clus, ' : ', betavals['mean', clus], ' (', betavals['var', clus], ')', sep = '')
# 
# 
#   p = ggplot(df, aes(X, fill = Cluster, y = ..count../sum(..count..))) +
#     geom_histogram(alpha = alpha, position = 'identity', binwidth = 0.01) +
#     scale_fill_manual(values = col.histogram, labels = extented.labels) +
#     labs(
#       title = tit,
#       subtitle = annotation,
#       x = "Observed Frequency", y = "") +
#     theme_classic(base_size = 10 * cex)
# 
# 
#   if(max.height.hist) {
#     yMax = max(ggplot_build(p)$data[[1]]$y)
#     p = p + ylim(0, yMax)
#   }
# 
# 
#   df.m = data.frame(variable = labels.betas, Mean = x$beta['mean', labels.betas])
#   rownames(df.m) = df.m$variable
# 
#   # Add densities to the plot
#   vvv = lapply(vvv, function(w) {w$y = w$y * 0.01; w}) # Scale density wrt binwidth
# 
#   for(i in seq(labels)) p = p + geom_line(data = vvv[[i]], aes(y = y, x = X), colour = col.lines[i])
# 
#   # Add means for Beta components
#   for(i in labels.betas) p = p +  geom_vline(data = df.m[i ,], aes(xintercept = Mean), colour = col.lines[i], linetype = "longdash")
# 
#   # Annotate convergency value
#   if(x$status) p = p + annotate("text", x = .9, y = yMax , label = paste(x$fit.type, ': CONVERGED'), size = 3, colour = 'darkgreen')
#   else  p = p + annotate("text", x = .9, y = yMax, label = paste(x$fit.type, ': NOT CONVERGED'), size = 3, colour = 'red')
# 
#   p
# }
# 
# 
# 
# .getIconGGplotHistogram = function(x, palette = 'Spectral',
#                                    tail.color = c('gainsboro', 'darkgray'),
#                                    alpha = .5, max.height.hist = TRUE, cex = 1,
#                                    bg.color = 'ivory2')
# {
#   domain = seq(0, 1, 0.01)
# 
#   labels = names(x$pi)
#   labels.betas = colnames(x$beta)
# 
#   pi = x$pi
#   ICL = round(as.numeric(x$scores$ICL), 2)
#   NLL = round(as.numeric(x$scores$NLL), 2)
#   K = as.numeric(x$K)
# 
#   col = RColorBrewer::brewer.pal(
#     RColorBrewer::brewer.pal.info[palette, 'maxcolors'], palette)
#   col = colorRampPalette(col)(x$Kbeta)
# 
#   col.histogram = c(tail.color[1], col)
#   col.lines = c(tail.color[2], col)
# 
#   names(col.histogram) = names(col.lines) = labels
# 
# 
#   # Plot 1 -- main histogram
#   df = data.frame(X = x$X, Cluster = x$labels, Color = col[x$labels])
#   vvv = lapply(1:x$K,
#                function(w)
#                  data.frame(X = domain,
#                             Cluster = labels[w],
#                             y = ddbpmm(x, data = domain, components = w, log = FALSE)))
#   names(vvv) = labels
# 
#   # Histogram coloured according to clustering assignments
# 
#   p = ggplot(df, aes(X, fill = Cluster, y = ..count../sum(..count..))) +
#     geom_histogram(alpha = alpha, position = 'identity', binwidth = 0.01) +
#     scale_fill_manual(values = col.histogram, labels = labels) +
#     theme_void(base_size = 10 * cex) +
#     guides(fill = FALSE)
# 
# 
#   if(max.height.hist) {
#     yMax = max(ggplot_build(p)$data[[1]]$y)
#     p = p + ylim(0, yMax)
#   }
# 
# 
#   df.m = data.frame(variable = labels.betas, Mean = x$beta['mean', labels.betas])
#   rownames(df.m) = df.m$variable
# 
#   # Add densities to the plot
#   vvv = lapply(vvv, function(w) {w$y = w$y * 0.01; w}) # Scale density wrt binwidth
# 
#   for(i in seq(labels)) p = p + geom_line(data = vvv[[i]], aes(y = y, x = X), colour = col.lines[i])
# 
#   p
# }


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

.plot.goodness_of_fit = function(x, binning = 1e-2, TOP = 5, cex = 1) 
{
  points = lapply(x$runs, function(w) mobster:::.compute_fit_sqerr(w, binning = binning))
  points = lapply(seq_along(points), function(w) { 
    points[[w]]$K = x$fits.table$K[w]
    points[[w]]$tail = x$fits.table$tail[w]
    points[[w]]$run = w
    
    
    points[[w]] })
  
  points = points[1:TOP]
  points = Reduce(bind_rows, points)
  
  
  ggplot(points, aes(x = x, y = cum.y, fill = factor(run), color = factor(run))) +
    geom_line(show.legend = TRUE) +
    theme_light(base_size =  8 * cex) +
    # guides(colour = FALSE) +
    xlab('VAF') + 
    ylab("SSE") +
    labs(
      title = bquote(bold('Goodnees of fit')), 
      subtitle = paste0('Binwidth = ', binning)
    ) +
    guides(fill = FALSE, color = guide_legend(title = "Solution")) + 
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      panel.background = element_rect(fill = 'white')
    )
}



############### Plot 1 -- main histogram for a MOBSTER fit
.mobster_fit_histogram = function(
  x,
  alpha,
  binwidth,
  palette, 
  tail.color,
  cex,
  annotate
) 
{
  # Fit variables
  domain = seq(0, 1, binwidth)
  
  labels = names(.params_Pi(x))
  labels.betas = mobster:::.params_Beta(x)$cluster
  
  pi = mobster:::.params_Pi(x)
  K = as.numeric(x$K)
  
  # Load colors
  colors = mobster:::getColors_model(x, alpha = alpha, palette = palette, tail.color = tail.color)
  
  # Main plotting data
  x$data$cluster = factor(x$data$cluster, levels = names(colors))
  
  ##### Text for the plot -- convergence status below the plot
  # - number of steps
  # - epsilon
  # - SSE
  conv.steps = length(x$all.NLL)
  
  conv.epsilon = ifelese(
    conv.steps < 2,
    0,
    formatC(
      abs(rev(x$all.NLL)[1] - rev(x$all.NLL)[2]),
      format = "e", 
      digits = 0
    )
  )

  SSE_fit = .compute_fit_sqerr(x, binning = binwidth)$cum.y
  conv.sse = formatC(max(SSE_fit), format = "e", digits = 3)
  
  label.fit = bquote(
    .(x$fit.type) *
      " ("*omega*" = "* .(conv.steps) ~ 'steps; ' * epsilon ~ '=' ~ .(conv.epsilon) *
      "; SSE"~.(conv.sse)*')')
  
  ##### Main ggplot object -- histogram
  hist_pl = ggplot(x$data, 
                   aes(
                     VAF, 
                     fill = factor(cluster), 
                     y = ..count.. /sum(..count..))) +
    geom_histogram(alpha = alpha,
                   position = 'identity',
                   binwidth = binwidth) +
    scale_fill_manual(values = colors, labels = names(colors)) +
    guides(fill = guide_legend(title = "MOBSTER\nClusters")) +
    labs(
      title = bquote(bold(.(histogram.main))),
      subtitle = annotation,
      caption = label.fit,
      x = "Observed Frequency",
      y = "Density [a.u.]"
    ) +
    theme_classic(base_size = 10 * cex) +
    geom_vline(
      xintercept = min(x$data$VAF),
      colour = 'black',
      linetype = "longdash",
      size = 05 * cex
    ) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex, "cm"),
      panel.background = element_rect(fill = 'white'),
      plot.caption = element_text(color = ifelse(x$status, "darkgreen",  "red"))
    )
  
  ##### Mixture compoenent density on top of the histogram
  densities = mobster:::template_density(
    x, 
    x.axis = domain, 
    binwidth = binwidth, 
    reduce = TRUE)
  
  # Add the trace and the mean of each component
  hist_pl = hist_pl + 
    geom_line(data = densities, aes(y = y, x = x, color = factor(cluster)), size = 1 * cex) +
    scale_color_manual(values = colors, labels = names(colors)) +
    guides(color = FALSE)
  
  Beta_peaks = x$Clusters %>%
    dplyr::filter(type == 'Mean', cluster != 'Tail')
  
  hist_pl = hist_pl + 
    geom_vline(data = Beta_peaks, aes(xintercept = fit.value, color = cluster), linetype = "longdash")
  
  # Rescale y in case the density goes too high (which could happen with the power law)
  if(max.height.hist) 
  {
    yMax = max(ggplot_build(hist_pl)$data[[1]]$y)
    hist_pl = hist_pl + ylim(0, yMax)
  }
  
  ##### Overall density and SSE
  # - cumulative density
  densities = tibble::as_tibble(densities)
  densities = densities %>% group_by(x) %>% summarise(y = sum(y), cluster = 'f(x)')
  
  m = max(densities$y, na.rm = TRUE)
  
  hist_pl = hist_pl +
    geom_line(data = densities %>% mutate(y = y + m * 0.02), # offset by 2% of the max to avoid overlays
              aes(y = y, x = x), 
              color = 'black',
              alpha = .8,
              size = .5 * cex,
              linetype = 'dashed',
              inherit.aes = FALSE)
  
  # SSE
  error = .compute_fit_sqerr(x, binning = binwidth)
  
  # scale to percentage for plotting
  me = max(SSE_fit$cum.y)
  error = error %>% mutate(cum.y = (cum.y / me) * m )
  
  hist_pl = hist_pl +
    geom_line(data = SSE_fit, aes(y = cum.y, x = x), 
              color = 'darkgray',
              alpha = 1,
              size = .2 * cex,
              linetype = 'dashed',
              inherit.aes = FALSE) +
    scale_y_continuous(sec.axis = sec_axis(~./ m * 100, name = "SSE [cumulative %]"))  
  
  # Annotation of input entries
  if(!is.null(annotate) & all(c("VAF", 'label') %in% colnames(annotate)))
  {
    # Position the point at coord x = VAF and y = density
    m = max(densities$y, na.rm = TRUE)
    
    annotate$y = round(annotate$VAF/binwidth)
    annotate$y = densities$y[annotate$y] + m * 0.02
    
    hist_pl = hist_pl +
      geom_label_repel(data = annotate, 
                       aes(
                         x = VAF, 
                         y = y,
                         label = label, 
                         color = factor(cluster)
                         # fill = factor(cluster)
                       ),
                       size = 1.5 * cex,
                       inherit.aes = FALSE,
                       box.padding = 0.95, 
                       segment.size = .2 * cex, force = 1) +
      geom_point(data = annotate, aes(x = VAF, y = y,  color = factor(cluster)), 
                 size = .3 * cex, alpha = 1, 
                 inherit.aes = FALSE)
  }
  
  hist_pl
}

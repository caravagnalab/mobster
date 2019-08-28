#' Co-clustering probability from non-parametric bootstrap 
#' 
#' @description Compute co-clustering probability from non-parametric bootstrap resamples
#' 
#' This function can be applied only to nonparametric bootstrap resamples.
#'
#' @param resamples Data resampled
#' @param fits Fits (one per resample)
#'
#' @return Nothing
#' @export
#'
#' @examples
mobster_plt_bootstrap_cocluster = function(x, resamples, fits, background.color = 'white') 
{
 
  
  # normalize values
  co.clustering = co.clustering/n
  
  colors = mobster:::scols(seq(0.1, 1, 0.05), palette = 'YlGnBu')
  colors = c(`0` = background.color, colors)
  
  ann.colors = mobster:::getColors_model(x)
  ann.VAF =  mobster:::scols(seq(0, 1, 0.05), palette = 'Oranges')
  
  # pheatmap::pheatmap(co.clustering[1:500, 1:500],
  figure = pheatmap::pheatmap(co.clustering,
                              main = 'Co-clustering (non parametric bootstrap)',
                              cluster_rows = FALSE,
                              cluster_cols = FALSE, 
                              show_rownames = FALSE,
                              show_colnames = FALSE,
                              border = NA,
                              color = colors,
                              annotation_row = ordered.data,
                              annotation_col = ordered.data, 
                              silent = TRUE,
                              annotation_colors = list(cluster = ann.colors, VAF = ann.VAF)
  )
  
  return(list(figure = figure, co.clustering = co.clustering))
}


#' Plots bootstrap fits of model and parameters.  
#' 
#' @description Plots model's stability (percentage of model
#' fits across resamples) and bootstrap estimates of the parameters' boxplots. 
#' 
#' This function can be applied also to both parametric and nonparametric bootstrap.
#'
#' @param resamples Data resampled
#' @param fits Fits (one per resample)
#' @param tail.color 
#'
#' @return A figure
#' @export
#'
#' @examples
mobster_plt_bootstrap_params = function(resamples, fits, tail.color = 'darkgray', cex = 1, palette = "Set1") 
{
  
  
  
  n = length(resamples)
  
  res = NULL
  
  ########################## Analyze outputs
  # -- model matching
  tab = lapply(seq(fits),
               function(w)
                 data.frame(
                   `resample` = w,
                   `tail` = fits[[w]]$fit.tail,
                   `K` =  fits[[w]]$Kbeta,
                   `Model` = paste0(
                     'K = ',
                     fits[[w]]$Kbeta,
                     ifelse(fits[[w]]$fit.tail,
                            " with tail",
                            " without tail")
                   ),
                   stringsAsFactors = FALSE
                 ))
  
  tab = Reduce(rbind, tab)
  
  # Extract model frequency
  tab.counts = as.data.frame(table(tab$Model))
  colnames(tab.counts) = c("Model", "Freq")
  tab.counts$Freq = tab.counts$Freq / n
  
  tab = tab[, c('tail', 'K', 'Model')]
  tab = tab[!duplicated(tab),]
  tab = merge(tab, tab.counts, by = 'Model')
  
  tab$tail = ifelse(tab$tail, 'With tail', 'Without tail')
  
  pio::pioTit("Model frequency")
  pio::pioDisp(tab)
  res = append(res, list(tab))
  
  # plot it
  pl_modelFreq = ggplot(data = tab, aes(x = factor(K), y = Freq, fill = Model)) +
    geom_bar(stat = "identity", position = 'dodge') +
    ylim(0, 1) +
    theme_light(base_size = 8) +
    facet_grid(K ~ tail) +
    labs(
      title = bquote(bold("Model frequency")),
      subtitle = paste0("Bootstrap resamples (n = ", n, ")"),
      x = 'Number of Beta components',
      y = 'Frequency'
    ) +
    guides(fill = FALSE)
  
  ########################## Analyze outputs
  # -- tail statistics
  tab = lapply(seq(fits),
               function(w)
                 cbind(fits[[w]]$Clusters, `resample` = w))
  
  tab = tibble::as_tibble(Reduce(rbind, tab))
  
  # empty plot -- in case there is no tail
  pl_tailParam = ggplot(data.frame()) +  geom_blank()
  
  # get outputs from the fits
  shape.stat = tab %>%
    filter(cluster == 'Tail',
           type == 'Shape' | type == 'Scale' | type == 'Mixing proportion')
  
  shape.stat = reshape2::melt(shape.stat[, c('type', 'fit.value')], id = 'type')
  shape.stat$value = as.numeric(shape.stat$value)
  
  pl_tailParam = ggplot(data = shape.stat, aes(
    x = factor(type),
    y = value,
  )) +
    geom_violin(alpha = .8, color = NA, fill=tail.color) + 
    geom_boxplot(alpha = 1, width = 0.2) +
    theme_light(base_size = 8 * cex) +
    facet_wrap( ~ type, scale = 'free') +
    labs(
      title = bquote(bold("Tail parameters")),
      subtitle = paste0("Bootstrap resamples (n = ", n, ")"),
      y = "Fit value",
      x = "Parameter"
    ) + 
    guides(fill = FALSE)
  
  shape.stat = as_tibble(shape.stat)
  
  out = shape.stat %>%
    filter(variable == 'fit.value') %>%
    group_by(type) %>%
    summarise(
      quantile_0.05 = quantile(value, .05),
      quantile_0.95 = quantile(value, .95)
    )
  
  pio::pioTit("Tail fits : [.05, .95] quantiles from bootstrap replicates")
  print(out)
  res = append(res, list(out))
  
  
  ########################## Analyze outputs
  # -- Beta statistics
  #
  # Note: we reorder the clusters IDs to match labels across resamples
  tab = lapply(seq(fits),
               function(w) {
                 # get clusters
                 cluster.table = fits[[w]]$Clusters %>%
                   filter(cluster != 'Tail')
                 
                 # sort the means -- might lead to a different naming
                 sorted.means = cluster.table %>%
                   filter(type == 'Mean')
                 
                 sorted.means = sorted.means[order(sorted.means$fit.value, decreasing = TRUE),]
                 
                 # create a map to the new names
                 new.names = paste0("C", 1:nrow(sorted.means))
                 names(new.names) = sorted.means$cluster
                 
                 # push new names in
                 new.cluster = NULL
                 for (s in cluster.table$cluster)
                   new.cluster = c(new.cluster, new.names[s])
                 
                 cluster.table$cluster = new.cluster
                 
                 cbind(cluster.table, `resample` = w)
               })
  
  Beta.params = tibble::as_tibble(Reduce(rbind, tab))
  Beta.params$fit.value = as.numeric(Beta.params$fit.value)
  
  upq = Beta.params %>% group_by(type) %>% summarise(upper = quantile(fit.value, probs = .9))
  lq = Beta.params %>% group_by(type) %>% summarise(lower = quantile(fit.value, probs = .1))
  
  m = Beta.params %>% group_by(type, cluster) %>% summarise(fit.value = 
                                                              round(median(fit.value), digits = 4))
  
  out = Beta.params %>%
    group_by(type, cluster) %>%
    summarise(
      quantile_0.05 = quantile(fit.value, .05),
      quantile_0.95 = quantile(fit.value, .95)
    )
  pio::pioTit("Tail fits : [.05, .95] quantiles from bootstrap replicates")
  print(out)
  res = append(res, list(out))
  
  
  # Beta.params %>% filter(a < lq %>% filter(type == 'a') %>% pull(lower) )
  
  
  # plot
  pl_BetaParam = ggplot(data = Beta.params, aes(
    x = (cluster),
    y = fit.value,
    fill = factor(cluster)
  )) +
    geom_violin(alpha = .6, color = NA) + 
    geom_boxplot(alpha = 1, width = 0.2) +
    # scale_y_log10() + 
    theme_light(base_size = 8 * cex) +
    facet_wrap( ~ type, scale = 'free', nrow = 1) +
    labs(
      title = bquote(bold("Beta fits (cluster ids sorted by mean value)")),
      subtitle = paste0("Bootstrap resamples (n = ", n, ")"),
      y = "Fit value",
      x = "Cluster"
    ) +
    scale_fill_brewer(palette = palette) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(title = "Cluster"), color = FALSE)
  
  pl_BetaParam = pl_BetaParam + geom_label_repel(data = m, aes(label = format(fit.value, scientific = T, digits = 2)), size = 2 * cex)
  
  figure = ggpubr::ggarrange(
    ggpubr::ggarrange(pl_modelFreq, pl_tailParam, labels = c("A", "B")),
    pl_BetaParam,
    labels = c('', 'C'),
    nrow = 2
  )
  
  list(figure = figure, table = res)
}


#' Plots bootstrap fits of model and parameters.  
#' 
#' @description Plots model's stability (percentage of model
#' fits across resamples) and bootstrap estimates of the parameters' boxplots. 
#' 
#' This function can be applied also to both parametric and nonparametric bootstrap.
#'
#' @param resamples Data resampled
#' @param fits Fits (one per resample)
#' @param tail.color 
#'
#' @return A figure
#' @export
#'
#' @examples
mobster_plt_bootstrap = function(resamples, fits,
                                 
                                 fit, 
                                 bootstrap.fits,  
                                 alpha = 0.025,
                                 statistics = 'model_frequency',
                                 tail.color = 'darkgray', cex = 1, palette = "Set1") 
{
  
  boot.values = bootstrapped_statistics(fit = fit, bootstrap.fits = bootstrap.fits, alpha = alpha)
  n = length(bootstrap.fits)
  
  
  
  res = NULL
  
  # ########################## Analyze outputs
  # # -- model matching
  # tab = lapply(seq(fits),
  #              function(w)
  #                data.frame(
  #                  `resample` = w,
  #                  `tail` = fits[[w]]$fit.tail,
  #                  `K` =  fits[[w]]$Kbeta,
  #                  `Model` = paste0(
  #                    'K = ',
  #                    fits[[w]]$Kbeta,
  #                    ifelse(fits[[w]]$fit.tail,
  #                           " with tail",
  #                           " without tail")
  #                  ),
  #                  stringsAsFactors = FALSE
  #                ))
  # 
  # tab = Reduce(rbind, tab)
  # 
  # # Extract model frequency
  # tab.counts = as.data.frame(table(tab$Model))
  # colnames(tab.counts) = c("Model", "Freq")
  # tab.counts$Freq = tab.counts$Freq / n
  # 
  # tab = tab[, c('tail', 'K', 'Model')]
  # tab = tab[!duplicated(tab),]
  # tab = merge(tab, tab.counts, by = 'Model')
  # 
  # tab$tail = ifelse(tab$tail, 'With tail', 'Without tail')
  # 
  # pio::pioTit("Model frequency")
  # pio::pioDisp(tab)
  # res = append(res, list(tab))
  # 
  # # plot it
  # pl_modelFreq = ggplot(data = tab, aes(x = factor(K), y = Freq, fill = Model)) +
  #   geom_bar(stat = "identity", position = 'dodge') +
  #   ylim(0, 1) +
  #   theme_light(base_size = 8) +
  #   facet_grid(K ~ tail) +
  #   labs(
  #     title = bquote(bold("Model frequency")),
  #     subtitle = paste0("Bootstrap resamples (n = ", n, ")"),
  #     x = 'Number of Beta components',
  #     y = 'Frequency'
  #   ) +
  #   guides(fill = FALSE)
  
  ########################## Analyze outputs
  # -- tail statistics
  tab = lapply(seq(fits),
               function(w)
                 cbind(fits[[w]]$Clusters, `resample` = w))
  
  tab = tibble::as_tibble(Reduce(rbind, tab))
  
  # empty plot -- in case there is no tail
  pl_tailParam = ggplot(data.frame()) +  geom_blank()
  
  # get outputs from the fits
  shape.stat = tab %>%
    filter(cluster == 'Tail',
           type == 'Shape' | type == 'Scale' | type == 'Mixing proportion')
  
  shape.stat = reshape2::melt(shape.stat[, c('type', 'fit.value')], id = 'type')
  shape.stat$value = as.numeric(shape.stat$value)
  
  pl_tailParam = ggplot(data = shape.stat, aes(
    x = factor(type),
    y = value,
  )) +
    geom_violin(alpha = .8, color = NA, fill=tail.color) + 
    geom_boxplot(alpha = 1, width = 0.2) +
    theme_light(base_size = 8 * cex) +
    facet_wrap( ~ type, scale = 'free') +
    labs(
      title = bquote(bold("Tail parameters")),
      subtitle = paste0("Bootstrap resamples (n = ", n, ")"),
      y = "Fit value",
      x = "Parameter"
    ) + 
    guides(fill = FALSE)
  
  shape.stat = as_tibble(shape.stat)
  
  out = shape.stat %>%
    filter(variable == 'fit.value') %>%
    group_by(type) %>%
    summarise(
      quantile_0.05 = quantile(value, .05),
      quantile_0.95 = quantile(value, .95)
    )
  
  pio::pioTit("Tail fits : [.05, .95] quantiles from bootstrap replicates")
  print(out)
  res = append(res, list(out))
  
  
  ########################## Analyze outputs
  # -- Beta statistics
  #
  # Note: we reorder the clusters IDs to match labels across resamples
  tab = lapply(seq(fits),
               function(w) {
                 # get clusters
                 cluster.table = fits[[w]]$Clusters %>%
                   filter(cluster != 'Tail')
                 
                 # sort the means -- might lead to a different naming
                 sorted.means = cluster.table %>%
                   filter(type == 'Mean')
                 
                 sorted.means = sorted.means[order(sorted.means$fit.value, decreasing = TRUE),]
                 
                 # create a map to the new names
                 new.names = paste0("C", 1:nrow(sorted.means))
                 names(new.names) = sorted.means$cluster
                 
                 # push new names in
                 new.cluster = NULL
                 for (s in cluster.table$cluster)
                   new.cluster = c(new.cluster, new.names[s])
                 
                 cluster.table$cluster = new.cluster
                 
                 cbind(cluster.table, `resample` = w)
               })
  
  Beta.params = tibble::as_tibble(Reduce(rbind, tab))
  Beta.params$fit.value = as.numeric(Beta.params$fit.value)
  
  upq = Beta.params %>% group_by(type) %>% summarise(upper = quantile(fit.value, probs = .9))
  lq = Beta.params %>% group_by(type) %>% summarise(lower = quantile(fit.value, probs = .1))
  
  m = Beta.params %>% group_by(type, cluster) %>% summarise(fit.value = 
                                                              round(median(fit.value), digits = 4))
  
  out = Beta.params %>%
    group_by(type, cluster) %>%
    summarise(
      quantile_0.05 = quantile(fit.value, .05),
      quantile_0.95 = quantile(fit.value, .95)
    )
  pio::pioTit("Tail fits : [.05, .95] quantiles from bootstrap replicates")
  print(out)
  res = append(res, list(out))
  
  
  # Beta.params %>% filter(a < lq %>% filter(type == 'a') %>% pull(lower) )
  
  
  # plot
  pl_BetaParam = ggplot(data = Beta.params, aes(
    x = (cluster),
    y = fit.value,
    fill = factor(cluster)
  )) +
    geom_violin(alpha = .6, color = NA) + 
    geom_boxplot(alpha = 1, width = 0.2) +
    # scale_y_log10() + 
    theme_light(base_size = 8 * cex) +
    facet_wrap( ~ type, scale = 'free', nrow = 1) +
    labs(
      title = bquote(bold("Beta fits (cluster ids sorted by mean value)")),
      subtitle = paste0("Bootstrap resamples (n = ", n, ")"),
      y = "Fit value",
      x = "Cluster"
    ) +
    scale_fill_brewer(palette = palette) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(title = "Cluster"), color = FALSE)
  
  pl_BetaParam = pl_BetaParam + geom_label_repel(data = m, aes(label = format(fit.value, scientific = T, digits = 2)), size = 2 * cex)
  
  figure = ggpubr::ggarrange(
    ggpubr::ggarrange(pl_modelFreq, pl_tailParam, labels = c("A", "B")),
    pl_BetaParam,
    labels = c('', 'C'),
    nrow = 2
  )
  
  list(figure = figure, table = res)
}

# 
# plot_2D_cluster_density = function(x, s.x, s.y, b.c, col = 'red', cov = 100, cutoff = 0.05)
# {
#   binfit = x$fit.Binomial
#   rownames(binfit$theta_k) = x$samples
#   
#   # mixing proportion
#   pi = binfit$pi_k[b.c]
#   
#   # Binomial parameters
#   p.x = binfit$theta_k[s.x, b.c]
#   p.y = binfit$theta_k[s.y, b.c]
#   
#   # less than 1% do not plot anything
#   eplot = function() {
#     ggplot() + annotation_custom(text_grob(paste0('< ', cutoff), size = 10)) + theme(element_blank(), axis.text = element_blank())
#   }
#   
#   if(all(p.x < cutoff, p.y < cutoff)) return(eplot())
#   
#   # Otherwise plot
#   domain = seq(0, 1, 0.01)
#   
#   grid = expand.grid(domain, domain, stringsAsFactors = FALSE)
#   colnames(grid) = c('x', 'y')
#   
#   grid$nv.x = round(grid$x * cov)
#   grid$nv.y = round(grid$y * cov)
#   
#   mpdf = function(x, y, b.c) {
#     p.x = binfit$theta_k[s.x, b.c]
#     p.y = binfit$theta_k[s.y, b.c]
#     
#     pi * dbinom(x, size = cov, prob = p.x) *
#       dbinom(y, size = cov, prob = p.y)
#   }
#   
#   all.points = lapply(
#     b.c,
#     # colnames(binfit$theta_k),
#     function(cl)
#     {
#       cbind(grid,
#             group = cl,
#             pdf =  apply(grid, 1, function(w)
#               mpdf(w['nv.x'], w['nv.y'], b.c = cl)))
#     })
#   
#   all.points = Reduce(rbind, all.points)
#   
#   ggplot(all.points, aes(x = x, y = y)) +
#     geom_raster(aes(fill = pdf), interpolate = TRUE) +
#     theme_light(base_size = 10) +
#     facet_wrap(~group) +
#     scale_fill_gradient (low = 'gainsboro', high = col) +
#     theme(legend.position="none") +
#     labs(x = s.x, y = s.y) +
#     geom_hline(yintercept = cutoff, size = .3, color = 'red',linetype = "longdash") +
#     geom_vline(xintercept = cutoff, size = .3, color = 'red',linetype = "longdash")
# }

#' 
#' binomial2D_cluster_template_density = function(x, s.x, s.y, col = 'red', cov = 100, cutoff = 0.05)
#' {
#'   if (is.null(x$fit.Binomial))
#'     stop("Binomial clusters are not available!")
#'   
#'   binfit = x$fit.Binomial
#'   rownames(binfit$theta_k) = x$samples
#'   
#'   # mixing proportion
#'   pi = binfit$pi_k
#'   
#'   # purity (we have to adjust the binomial parameters for purity if we want to plto with adjusted VAF)
#'   purity = x$purity
#'   
#'   # Otherwise plot 
#'   domain = seq(0, 1, 0.01)
#'   
#'   grid = expand.grid(domain, domain, stringsAsFactors = FALSE)
#'   colnames(grid) = c('x', 'y')
#'   
#'   grid$nv.x = round(grid$x * cov)
#'   grid$nv.y = round(grid$y * cov)
#'   
#'   mpdf = function(x, y, b.c) {
#'     p.x = binfit$theta_k[s.x, b.c] / purity[s.x]
#'     p.y = binfit$theta_k[s.y, b.c] / purity[s.y]
#'     
#'     pi[b.c] * dbinom(x, size = cov, prob = p.x) *
#'       dbinom(y, size = cov, prob = p.y) 
#'   }
#'   
#'   all.points = lapply(
#'     colnames(binfit$theta_k),
#'     function(cl)
#'     {
#'       cbind(grid,
#'             group = cl,
#'             pdf =  apply(grid, 1, function(w)
#'               mpdf(w['nv.x'], w['nv.y'], b.c = cl)))
#'     })
#'   
#'   Reduce(rbind, all.points)
#' }
#' 
#' 
#' #' Title
#' #'
#' #' @param x 
#' #' @param palette 
#' #' @param ... 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' mobster_plt_2D_cluster_density = function(x, cov = 100, palette = "Set1", ...) {
#'   stopifnot(!is.null(x$fit.Binomial))
#'   
#'   # Get colors
#'   colors = mobster:::scols(names(x$fit.Binomial$pi_k), palette = palette)
#'   
#'   # Clusters names
#'   cls = x$fit.Binomial$pi_k
#'   
#'   # paired of samples
#'   pairs = combn(x$samples, 2)
#'   
#'   pl_list = NULL
#'   
#'   for(i in 1:ncol(pairs))
#'   {
#'     # Get a plot for each cluster
#'     pl = lapply(
#'       seq_along(cls), 
#'       function(c) plot_2D_cluster_density(
#'         x = x, 
#'         b.c = names(cls)[c], 
#'         s.x = pairs[1, i], 
#'         s.y = pairs[2, i], 
#'         col = colors[c], 
#'         ...)
#'     )
#'     
#'     header_pl = ggplot() + 
#'       annotation_custom(text_grob(paste0(pairs[1, i], ' vs ', pairs[2, i]), size = 20)) + 
#'       theme(element_blank(), axis.text = element_blank()) 
#'     
#'     
#'     # Display it as a row
#'     pl_fig = ggarrange(plotlist = append(list(header_pl), pl), ncol = length(cls) + 1, nrow = 1)
#'     pl_list = append(pl_list, list(pl_fig))
#'   }
#'   
#'   fg = ggarrange(plotlist = pl_list, ncol = 1, nrow = ncol(pairs))
#'   annotate_figure(fg, top = paste("Multivariate Binomial density for", x$mobster_global$description), 
#'                   bottom = paste0("Density at coverage", cov, "x"))
#' }

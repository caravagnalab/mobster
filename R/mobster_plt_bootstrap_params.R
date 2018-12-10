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
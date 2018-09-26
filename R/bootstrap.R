# lapply(.parametric_bootstrap_resamples(x, n = 2), function(w) hist(w$VAF, breaks = seq(0, 1, 0.01)))
.parametric_bootstrap_resamples = function(x, n = 100)
{
  stopifnot(n >= 1)
  
  data.size = nrow(x$data)
  
  lapply(1:n,
         function(w) {
           resample = rdbpmm(x, n = data.size)
           data.frame(id = 1:data.size, VAF = resample)
         })
}

.nonparametric_bootstrap_resamples = function(x, n = 100)
{
  stopifnot(n >= 1)
  
  data.size = nrow(x$data)
  
  lapply(1:n,
         function(w) {
           ids = sample(1:data.size, data.size, replace = TRUE)
           data.frame(
             id = 1:data.size,
             VAF = x$data$VAF[ids],
             original.id = ids
           )
         })
}


#' Parallel implementation of parametric/ nonparametric 
#' bootstrap for MOBSTER.
#'
#' @param x a fit by MOBSTER
#' @param n number of parametric resamples (drawn via \code{rdbpmm})
#' @param cores.ratio  ratio of cores to use
#' @param ... fit parameters for \code{mobster_fit}
#' @param palette plot palette used for Beta clusters
#' @param alpha transparency
#' @param tail.color for the tail's boxplot
#' @param bootstrap "parametric" or "nonparametric"
#'
#' @return Data from the fits, resamples and a plottable figure.
#' @export
#'
#' @examples
mobster_bootstrap = function(x,
                             n = 100,
                             bootstrap = 'parametric',
                             palette = 'Set1',
                             alpha = 1,
                             tail.color = 'darkgray',
                             cores.ratio = 0.8,
                             ...)
{
  pio::pioHdr(
    "MOBSTER bootstrap",
    toPrint = c(`Resamples` = paste(n),
                `Type of bootstrap` = bootstrap),
    prefix = '\t-'
  )
  
  stopifnot(bootstrap %in% c('parametric', 'nonparametric'))

  pio::pioTit(paste0("Bootstrapping for this MOBSTER model"))
  print(x)
  
  pio::pioTit(paste0("Creating ", bootstrap, " bootstrap resamples"))

  resamples = NULL
  
  # Get resamples -- n datasets with same size of x
  if (bootstrap == 'parametric')
    resamples = .parametric_bootstrap_resamples(x, n)
  
  if (bootstrap == 'nonparametric')
    resamples = .nonparametric_bootstrap_resamples(x, n)
  
  pio::pioTit("Running fits (might take some time)")
  
  # Setup clusters for parallel computing
  cl = .setup_parallel(cores.ratio = cores.ratio)
  
  # perform parallel inferences with custom parameters
  # -- disable internal parallelism
  # -- randomize seed
  fits = foreach(
    num = 1:n,
    .packages = c("crayon", "mobster"),
    .export = ls(globalenv(), all.names = TRUE)
  ) %dopar%
  {
    # best fit from resample
    mobster_fit(x = resamples[[num]],
                parallel = FALSE,
                seed = NULL,
                ...)$best
  }
  
  return(list(resamples = resamples, fits = fits))
}

#' Analyze fits by bootstrap. 
#' 
#' @description Computes model's stability (percentage of model
#' fits across resamples) and parameters' boxplots. 
#' 
#' This function can be applied also to both parametric and nonparametric bootstrap.
#'
#' @param resamples Data resampled
#' @param fits Fits (one per resample)
#'
#' @return A figure
#' @export
#'
#' @examples
mobster_bootstrap_params = function(resamples, fits) 
{
  n = length(resamples)
  
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
  
  # plot it
  pl_modelFreq = ggplot(data = tab, aes(x = factor(K), y = Freq, fill = Model)) +
    geom_bar(stat = "identity", position = 'dodge') +
    ylim(0, 1) +
    theme_light(base_size = 8) +
    facet_wrap( ~ tail, nrow = 1) +
    labs(
      title = "Model frequency",
      subtitle = paste0("Bootrstrap resamples n = ", n, " "),
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
      fill = factor(type)
    )) +
      geom_boxplot(alpha = .8, color = tail.color) +
      theme_light(base_size = 8) +
      facet_wrap( ~ type, scale = 'free') +
      labs(
        title = "Tail fits",
        subtitle = paste0("Bootrstrap resamples n = ", n, " "),
        y = "Fit value",
        x = "Parameter"
      ) + guides(fill = FALSE)
  
  
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
  
  # plot
  pl_BetaParam = ggplot(data = Beta.params, aes(
    x = (cluster),
    y = fit.value,
    fill = factor(cluster)
  )) +
    # geom_violin(alpha = .8) +
    geom_boxplot(alpha = .8) +
    theme_light(base_size = 8) +
    facet_wrap( ~ type, scale = 'free', nrow = 1) +
    labs(
      title = "Beta fits (cluster ids sorted by mean value)",
      subtitle = paste0("Bootrstrap resamples n = ", n, " "),
      y = "Fit value",
      x = "Parameter"
    ) +
    scale_fill_brewer(palette = palette) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(title = "Cluster"))
  
  figure = ggpubr::ggarrange(
    ggpubr::ggarrange(pl_modelFreq, pl_tailParam, labels = c("A", "B")),
    pl_BetaParam,
    labels = c('', 'C'),
    nrow = 2
  )
  
  figure
}



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
mobster_bootstrap_cocluster = function(x, resamples, fits, background.color = 'white') 
{
  # sum up occrences 
  .coocc = function(l, M) 
  {
    cluster.labels = unique(l)
    
    for (cl in cluster.labels) 
    {
      # A unique is for nonparametric bootstrap
      # where we resample the same samples twice
      cl.assignments = unique(names(l[l == cl]))
      if(is.null(cl.assignments) | length(cl.assignments) == 1) next;
      
      pairs = combn(cl.assignments, 2, simplify = F)
      
      for (p in 1:length(pairs)) 
      {
        M[pairs[[p]][1], pairs[[p]][2]] = M[pairs[[p]][1], 
                                            pairs[[p]][2]] + 1
        M[pairs[[p]][2], pairs[[p]][1]] = M[pairs[[p]][2], 
                                            pairs[[p]][1]] + 1
      }
    }
    M
  }
  
  ########################## Analyze outputs
  # -- Co-clustering probability
  #
  # Note: only with nonparametric bootstrap
  
  # number of resamples
  n = length(resamples)

  # number of mutations
  N = nrow(fits[[1]]$data)
  
  if(!('original.id' %in% colnames(fits[[1]]$data)))
    stop('Missing the original.id column in the bootstrap data! \n\nAre you sure this is
         a result from non parametric bootstrap?')
  
  # Still assuming that id is the row index of the VAF
  co.clustering = matrix(0, nrow = N, ncol = N)
  
  rownames(co.clustering) =
    colnames(co.clustering) = 1:N
  
  # Extract co-clustering labels
  for(w in seq(fits))
  {    
    cluster.results = fits[[w]]$data
    
    cluster.labels = cluster.results$cluster
    names(cluster.labels) = cluster.results$original.id
    
    co.clustering = .coocc(cluster.labels, co.clustering)
  }
  
  ordered.data = x$data[, c('VAF', 'cluster')]
  rownames(ordered.data) = 1:N
  
  # sort heatmap by cluster
  ordering = order(ordered.data$cluster)
  
  ordered.data = ordered.data[ordering, ]
  co.clustering = co.clustering[ordering, ordering]
  
  # normalize values
  co.clustering = co.clustering/n
  
  colors = mobster:::scols(seq(0.1, 1, 0.05), palette = 'YlGnBu')
  colors = c(`0` = background.color, colors)
  
  ann.colors = mobster:::getColors_model(x)
  ann.VAF =  mobster:::scols(seq(0, 1, 0.05), palette = 'Oranges')
  
  # pheatmap::pheatmap(co.clustering[1:500, 1:500],
  pheatmap::pheatmap(co.clustering,
                     main = 'Co-clustering (non parametric bootstrap)',
                     cluster_rows = FALSE,
                     cluster_cols = FALSE, 
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     border = NA,
                     color = colors,
                     annotation_row = ordered.data,
                     annotation_col = ordered.data,
                     annotation_colors = list(cluster = ann.colors, VAF = ann.VAF)
  )
}


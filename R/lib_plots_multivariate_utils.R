get_mixing_proportions = function(data, cluster)
{
  pi = as.vector(table(data[, cluster])) /
    sum(!is.na(data[, cluster]))

  pi.df = data.frame(sort(pi, decreasing = TRUE))
  colnames(pi.df) = cluster
  pi.df$Cluster = 1:nrow(pi.df)

  pi.df
}

show_VAF_per_segment =  function(segments, muts, cutoff, ...)
{
  plots = lapply(1:nrow(segments),
                 function(s) {
                   # map SNVs to this segment
                   map = map_SNV_2_CNA_segment(segments = segments,
                                               segment.id = s,
                                               muts = muts,
                                               ...)

                   map = map[, startsWith(colnames(map), 'VAF')]

                   map = reshape2::melt(map)
                   map = map[map$value > cutoff,]

                   mxlim = max(1, max(map$value, na.rm = T))

                   # lbl = paste0(
                   #
                   # )

                   ggplot(map, aes(value, fill = variable)) +
                     geom_histogram(binwidth = 0.01) +
                     facet_wrap( ~ variable, scale = 'free', nrow = 1) +
                     coord_cartesian(xlim = c(0, mxlim)) +
                     guides(fill = FALSE) +
                     labs(title = paste(segments[s, ], collapse = ',')) +
                     geom_vline(aes(xintercept = 0.5), colour = 'red', linetype = "longdash", size = .3) +
                     geom_vline(aes(xintercept = 0.25), colour = 'red', linetype = "longdash", size = .3) +
                     geom_vline(aes(xintercept = 0.33), colour = 'blue', linetype = "longdash", size = .3) +
                     geom_vline(aes(xintercept = 0.66), colour = 'blue', linetype = "longdash", size = .3)
                 })

  figure = ggpubr::ggarrange(
    plotlist = plots,
    nrow = nrow(segments),
    ncol = 1
  )

  figure
}

map_SNV_2_CNA_segment = function(segments,
                                 segment.id,
                                 muts,
                                 labels = list(
                                   Chromosome = c(`chromosome` = 'chr', `from` = 'from', `to` = 'to'),
                                   Mutations = c(`chromosome` = 'chr', `position` = 'pos')
                                 ))
{
  stopifnot(all(unlist(labels$Chromosome) %in% colnames(segments)))
  stopifnot(all(unlist(labels$Mutations) %in% colnames(muts)))

  map = NULL
  map = muts[muts[, labels$Mutations['chromosome']] == segments[segment.id, labels$Chromosome['chromosome']], ]
  map = map[map[, labels$Mutations['position']] >= segments[segment.id, labels$Chromosome['from']], ]
  map = map[map[, labels$Mutations['position']] <= segments[segment.id, labels$Chromosome['to']], ]

  map
}


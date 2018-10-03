#
#' @title Prepare multivariate inputs.
#'
#' @details
#'  Requirement for a multivariate analysis, a data frame with this format
#'
#' > head(data)
#' id  NV.B1  DP.B1       VAF.B1  NV.B2  DP.B2      VAF.B2       chr    pos
#' 241    62    123    0.5040650     54    112   0.4821429    pseudo      1
#' 253    46    100    0.4600000     57    108   0.5277778    pseudo      2
#' 264    65    119    0.5462185     63    120   0.5250000    pseudo      3
#' 275    48    112    0.4285714     68    129   0.5271318    pseudo      4
#' 286    45    106    0.4245283     66    117   0.5641026    pseudo      5
#' 297    59    116    0.5086207     69    136   0.5073529    pseudo      6
#
#' In the above dataset we have 2 samples, named "B1" and "B2". Input columns
#' "id", "chr" and "pos" are optional and are auto-generated in case they are
#' missing from the input dataframe.
#'
#' @param data input SNVs detected across all samples
#' @param CN Copy Number segments
#' @param samples samples ID
#' @param VAF.range VAF range for mutations to analyze
#' @param subsample VAF range for mutations to analyze
#'
#' @return data for input in sciClone
#' @export
#'
#' @examples
format_input = function(data, CN, samples, VAF.range = c(0.05, 0.6), subsample = NULL)
{
  cn = colnames(data)
  vaf.columns = paste0("VAF.", samples)

  if(! all(paste0("DP.", samples) %in% cn)) stop("Missing DP entries?")
  if(! all(paste0("NV.", samples) %in% cn)) stop("Missing DP entries?")
  if(! all(paste0("VAF.", samples) %in% cn)) stop("Missing DP entries?")

  # Fix missing data
  if(!("chr" %in% cn)) {
    message("Adding Chromosome 'chr' field to dataframe with value 'pseudo'")
    data$chr = 'pseudo'
  }

  if(!("pos" %in% cn)) {
    message("Adding Position 'pos' field to dataframe with value the row id")
    data$pos = as.numeric(1:nrow(data))
  }

  if(!("pos" %in% cn)) {
    message("Adding ID 'id' field to dataframe with value the row id")
    data$id = as.numeric(1:nrow(data))
  }

  if(all(is.null(CN))){
    message("Missing segments data, assuming everything is diploid and using start/ stops from row ids.")
    CN = data.frame(chr = 'pseudo',
                    start = as.numeric(seq_len(nrow(data))),
                    stop = as.numeric(seq_len(nrow(data))),
                    segment_mean = 2)
  }

  pio::pioTit(paste("Original input dataframe - n =", nrow(data)))
  print(tibble::as.tibble(data))




  if(!is.null(subsample))
  {
    pio::pioTit(paste0("Subsampling to ", subsample, " SNVs"))
    if(nrow(data) > subsample)
    {
      data = data[sample(1:nrow(data), subsample), ]
      message("Subsampling carried out succesfully")
    }
  }


  # Disassamble data
  input.SNVs = lapply(samples, function(w){

    # A single sample
    this.data = data[, paste0(c("NV.", "DP.", "VAF."), w)]

    # subset data in the format for sciClone
    this.data$ref_reads = as.numeric(this.data[, 2] - this.data[, 1]) # coverage - nv
    this.data$var_reads = as.numeric(this.data[, 1]) # nv
    this.data$vaf = this.data[, 3] * 100 # VAF scaled in [0, 100]

    # SNVs in sciClone format
    this.data = this.data[, c('ref_reads', 'var_reads', 'vaf')]

    # this.data = cbind(this.data, data[, c('id', 'chr', 'pos')])
    this.data = cbind(data[, c('chr', 'pos')], this.data)

    colnames(this.data) = c('chr',	'start', 'refCount', 'varCount', 'VAF')


    this.data
  })

  # CN are replicated across samples
  input.CNs = lapply(seq(input.SNVs), function(w) CN)

  names(input.SNVs) = names(input.CNs) = samples


  inputs = list(data = data, SNVs = input.SNVs, CNs = input.CNs, samples = samples)

  inputs
}




# data is obtained from the above function
#' Title
#'
#' @param data
#' @param samples
#' @param minimumDepth
#'
#' @return
#' @export
#'
#' @examples
sciClone_fit = function(data, samples, minimumDepth = 30, maximumClusters = 10, ...) {

  library(sciClone)
  sciClone.fit = sciClone(
    vafs = data$SNVs,
    copyNumberCalls = data$CNs,
    sampleNames = samples,
    clusterMethod = 'bmm',
    maximumClusters = maximumClusters,
    verbose = FALSE,
    minimumDepth = minimumDepth, ...)

  factor.values = sort(unique(sciClone.fit@vafs.merged$cluster))
  data$data$sciClone.cluster = factor(sciClone.fit@vafs.merged$cluster, levels = factor.values)

  data$sciClone.fit = sciClone.fit

  data
}

#' Title
#'
#' @param data
#' @param samples
#' @param minimumDepth
#'
#' @return
#' @export
#'
#' @examples
MOBSTER_sciClone_fit = function(data, samples, VAF.adjustment.max = TRUE, minimumDepth = 30, projection.tails = "Global", ...) {

  cnames = paste("MOBSTER", samples, "cluster", sep = '.')
  nc = ncol(data$data) + 1
  for(s in samples)
    data$data = cbind(data$data, NA)
  colnames(data$data)[nc:(nc+length(samples)-1)] = cnames

  pio::pioTit("Running MOBSTER on each input samples")

  MOBSTER.fits = lapply(
    samples,
    function(s)
    {
      pio::pioTit(paste("MOBSTER fit for", s))

      # Run MOBSTER on each sample
      inputs = data$SNVs[[s]]
      inputs = inputs$VAF/100

      # Gain "space" in the wider frequency spectrum
      if(VAF.adjustment.max) inputs = inputs/max(inputs)

      inputs = inputs[inputs > 0]

      if(length(inputs) == 0) stop("0 entries with VAF > 0 in one MOBSTER input??")

      library(dbpmm)

      dbpmm.fit(
        X = inputs,
        K = 1:3,
        samples = 3,
        init = 'peaks',
        tail = c(TRUE, FALSE),
        fit.type = 'MM',
        maxIter = 6000,
        epsilon = 1e-8,
        parallel = FALSE
      )

    })
  names(MOBSTER.fits) = samples
  data$MOBSTER.fits = MOBSTER.fits

  data$best.MOBSTER = lapply(MOBSTER.fits, function(w) w$best)
  names(data$best.MOBSTER) = samples

  pio::pioTit("MOBSTER fits")
  print(data$best.MOBSTER)

  # MOBSTER clustering assignments
  for(s in seq(samples))
  {
    data$data[
      data$SNVs[[s]]$VAF > 0,
      cnames[s]] = data$best.MOBSTER[[s]]$labels
  }
  head(data$data)

  pio::pioTit(paste("Tail detection -- projection type:", projection.tails))
  pio::pioDisp(data$data[, c(cnames)])

  # Projected SNVs -- 2 possible protocols
  data$projected.SNVs = data$SNVs

  prj.names = paste0('VAF.projected.', samples)

  if(projection.tails == 'Local')
  {
    for(s in seq(samples))
    {
      original = data$projected.SNVs[[s]]$VAF/100
      clusters = data$data[, cnames[s]]

      original[is.na(clusters)] = 0      # projection
      original[clusters == "Tail"] = 0   # projection

      data$projected.SNVs[[s]]$VAF = original
      data$data = cbind(data$data, original)
      colnames(data$data)[ncol(data$data)] = prj.names[s]
    }
  }

  if(projection.tails == 'Global')
  {
    any.tail = apply(data$data[, cnames], 1, function(x) any(x == 'Tail', na.rm = TRUE) )

    # Some had VAF < in the original data, because they did not pass QC or minimum VAF.
    # De facto these are considered tails, as they are not for sure a cluster. Thus one
    # entry with negative VAF is enough to discard the SNV

    samples.names = paste0('VAF.', samples)
    any.negVAF = apply(data$data[, samples.names], 1, function(x) any(x < 0, na.rm = TRUE) )

    data$data$any.negVAF = any.negVAF
    data$data$any.tail = any.tail | any.negVAF

    for(s in seq(samples))
    {
      data$projected.SNVs[[s]]$VAF[any.tail] = 0 # projection
      data$projected.SNVs[[s]]$refCount[any.tail] = 0 # projection


      data$data = cbind(data$data, data$projected.SNVs[[s]]$VAF/100)
      colnames(data$data)[ncol(data$data)] = prj.names[s]
    }
  }

  # we can further remove those that are all 0s everywhere
  idx.remove = apply(data$data[, prj.names], 1, function(w) all(w==0))

  for(s in seq(samples))
  {
    data$projected.SNVs[[s]] = data$projected.SNVs[[s]][!idx.remove, ]
  }

  # Clustering projected read counts
  pio::pioTit("sciClone fit on projected read counts")

  library(sciClone)
  MOBSTER.sciClone.fit = sciClone(
    vafs = data$projected.SNVs,
    copyNumberCalls = data$CNs,
    sampleNames = samples,
    clusterMethod = 'bmm',
    verbose = TRUE,
    minimumDepth = minimumDepth, ...)

  factor.values = sort(unique(MOBSTER.sciClone.fit@vafs.merged$cluster))
  data$data$MOBSTER.sciClone.cluster[!idx.remove] =
    factor(MOBSTER.sciClone.fit@vafs.merged$cluster, levels = factor.values)

  data$MOBSTER.sciClone.fit = MOBSTER.sciClone.fit

  data
}



#' Grid plot for multivariate analysis.
#'
#' @param data
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
plot_grid = function(data, samples, below.cluster = 'sciClone.cluster', top.cluster = 'MOBSTER.sciClone.cluster',
                     diagonal = 'MOBSTER') {

  pio::pioHdr("MOBSTER - Multivariate grid plot",
              c(`Bottom block` = below.cluster,
                `Top block` = top.cluster,
                `Diagonal` = paste(diagonal)
                )
   )

  require(ggplot2)

  MB = plots = NULL

  #  Diagonal
  if(!is.null(diagonal) & diagonal == 'MOBSTER')
    MB = plot_diagonal_MOBSTER(data$best.MOBSTER, samples)

  if(!is.null(below.cluster))
  {
    for(s in seq(samples)) {
      for(w in s:length(samples)) {
        if(s != w) {

          pl = plot_2DVAF(
            data$data,
            x = paste0('VAF.', samples[s]),
            y = paste0('VAF.', samples[w]),
            cluster = below.cluster)

          plots = append(plots, list(pl))
        }
      }
    }
  }

  if(!is.null(top.cluster))
  {
    for(s in seq(samples)) {
      for(w in s:length(samples)) {
        if(s != w) {

          pl = plot_2DVAF(
            data$data,
            x = paste0('VAF.projected.', samples[s]),
            y = paste0('VAF.projected.', samples[w]),
            cluster = top.cluster)

          plots = append(plots, list(pl))
        }
      }
    }
  }

  layout = matrix(0, ncol = length(samples), nrow = length(samples))
  diag(layout) = 1:length(samples)

  combs = length(samples) * (length(samples)-1) / 2
  layout[lower.tri(layout)] = (1:combs) + length(samples)
  layout[upper.tri(layout)] = (1:combs) + length(samples) + combs

  .multiplot(plotlist = plots, layout = layout)
}

#' Grid plot for multivariate analysis.
#'
#' @param data
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
plot_grid_plain = function(data, samples, VAF.range = c(0.05, 0.6), cex = 1, title = 'VAF pairwise plot') {

  pio::pioHdr('MOBSTER - plot 2D VAF',
              toPrint = c(
                `Sample IDs` = paste(samples, collapse = ', '),
                `VAF range` = paste(VAF.range, collapse = ' -- ')
              ))

  #  Diagonal is just a marginal freq historgram
  plots = lapply(
    samples,
    function(w) {

      w = paste0('VAF.', w)

      ggplot(data[data[, w] >= VAF.range[1] & data[, w] <= VAF.range[2], ],
             aes(eval(parse(text = w)), y = ..count../sum(..count..))) +
        geom_histogram(alpha = .8, position = 'identity', binwidth = 0.01) +
        labs(
          title = w,
          x = "Observed Frequency", y = "Density") +
        theme_classic(base_size = 10 * cex)
    })

  require(ggplot2)

  # 2D freq historgram
  for(s in seq(samples)) {
    for(w in s:length(samples)) {
      if(s != w) {

        pl = plot_2DVAF(
          data,
          x = paste0('VAF.', samples[s]),
          y = paste0('VAF.', samples[w]),
          cluster = NULL,
          VAF.range = VAF.range)

        plots = append(plots, list(pl))
      }
    }
  }

  layout = matrix(0, ncol = length(samples), nrow = length(samples))
  diag(layout) = 1:length(samples)

  combs = length(samples) * (length(samples)-1) / 2
  layout[lower.tri(layout)] = (1:combs) + length(samples)

  .multiplot(plotlist = plots, layout = layout, title = title)
}


#' Title
#'
#' @param x
#' @param title
#'
#' @return
#' @export
#'
#' @examples
plot_data = function(x, VAF.columns, DP.columns, NV.columns, title = 'Raw data: VAF (full and low-freq.), DP, NV')
{
  x.orig = nrow(x)
  x = x[complete.cases(x), ]
  x.comp = nrow(x)

  message("Droppping NA values reduced the dataset from ", x.orig, ' to ', x.comp, 'entries')


  VAF.df = reshape::melt(x[, VAF.columns])
  DP.df = reshape::melt(x[, DP.columns])
  NV.df = reshape::melt(x[, NV.columns])

  require(ggplot2)

  min.VAF.above0 = min(VAF.df$value[VAF.df$value > 0], na.rm = T)

  VAF.df = VAF.df[VAF.df$value >= min.VAF.above0, ]

  max.VAF = max(VAF.df$value)
  if(max.VAF < 1) max.VAF = 1

  VAF = ggplot(VAF.df, aes(value, fill = variable)) +
    geom_histogram(binwidth = 0.01) +
    facet_wrap(~variable, nrow = 1) +
    guides(fill = FALSE, color = FALSE) +
    xlim(0, max.VAF)

  VAF.low.range.left =  0.13

  VAF.df = VAF.df[VAF.df$value <= VAF.low.range.left, ]
  VAF.low = ggplot(VAF.df, aes(value, fill = variable)) +
    geom_histogram(binwidth = 0.001) +
    facet_wrap(~variable, nrow = 1) +
    guides(fill = FALSE, color = FALSE) +
    geom_vline(xintercept = seq(0.01, 0.05, 0.01), linetype = 'longdash', color = 'gray')

  DP = ggplot(DP.df, aes(value, fill = variable)) +
    geom_histogram(binwidth = 1) +
    facet_wrap(~variable, nrow = 1) +
    guides(fill = FALSE, color = FALSE)

  NV.df = NV.df[NV.df$value > 0, ]
  NV = ggplot(NV.df, aes(value, fill = variable)) +
    geom_histogram(binwidth = 1) +
    facet_wrap(~variable, nrow = 1) +
    guides(fill = FALSE, color = FALSE)


  figure = ggpubr::ggarrange(
    VAF,
    VAF.low,
    DP,
    NV,
    ncol = 1,
    nrow = 4)

  figure = ggpubr::annotate_figure(figure, top = title)

  figure
}



plot_diagonal_MOBSTER = function(MB, samples, cex = 1)
{
  lapply(seq(MB),
         function(w)
           plot(
             MB[[w]],
             silent = TRUE,
             palette = 'Set1',
             cex = cex,
             histogram.main = paste("MOBSTER ", samples[w])
           )$mainHist)
}


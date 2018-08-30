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
#'
#' @return data for input in sciClone
#' @export
#'
#' @examples
sciClone_input = function(data, CN, samples)
{
  cn = colnames(data)

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

  cat("Original dataframe\n")
  print(tibble::as.tibble(data))

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


  inputs = list(data = data, SNVs = input.SNVs, CNs = input.CNs)

  inputs
}

# input = sciClone_input(data, samples = c("B1", "B2"), CN = NULL)
#
# lapply(input$SNVs, head)
# lapply(input$CNs, head)


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
sciClone_fit = function(data, samples, minimumDepth = 30, ...) {

  library(sciClone)
  sciClone.fit = sciClone(
    vafs = data$SNVs,
    copyNumberCalls = data$CNs,
    sampleNames = samples,
    clusterMethod = 'bmm',
    verbose = TRUE,
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
MOBSTER_sciClone_fit = function(data, samples, minimumDepth = 30, projection.tails = "Global", ...) {

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
      # Run MOBSTER on each sample
      inputs = data$SNVs[[s]]
      inputs = inputs$VAF/100
      inputs = inputs[inputs > 0]

      if(length(inputs) == 0) stop("0 entries with VAF > 0 in one MOBSTER input??")

      library(dbpmm)

      dbpmm::dbpmm.fit(
        X = inputs,
        K = 1:3,
        samples = 3,
        init = 'peaks',
        tail = c(TRUE, FALSE),
        fit.type = 'MM',
        maxIter = 6000,
        epsilon = 1e-8,
        parallel = FALSE,
        top = 10
      )

    })
  names(MOBSTER.fits) = samples
  data$MOBSTER.fits = MOBSTER.fits

  best.MOBSTER = lapply(MOBSTER.fits, function(w) w[[1]])
  names(best.MOBSTER) = samples
  data$best.MOBSTER = best.MOBSTER

  pio::pioTit("MOBSTER fits")
  print(best.MOBSTER)

  # MOBSTER clustering assignments
  for(s in seq(samples))
  {
    data$data[
      data$SNVs[[s]]$VAF > 0,
      cnames[s]] = best.MOBSTER[[s]]$labels
  }
  head(data$data)

  cat("\nTail detection -- projection type: ", projection.tails)
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
    data$data$any.tail = any.tail

    for(s in seq(samples))
    {
      data$projected.SNVs[[s]]$VAF[any.tail] = 0 # projection
      data$data = cbind(data$data, data$projected.SNVs[[s]]$VAF)
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
plot_2DVAF = function(data, x, y, cluster, marginal = FALSE) {

  data = data[!is.na(data[, cluster]), ]
  data[, cluster] = paste(data[, cluster])

  require(ggplot2)

  p = ggplot(data, aes(x = eval(parse(text = x)), colour = eval(parse(text = cluster)), y = eval(parse(text = y)))) +
    # theme_minimal() +
    theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              # panel.grid.minor = element_blank(),
              axis.line = element_line(size = 0.5, linetype = "solid",
                                       colour = "black")) +
    # geom_density_2d(aes(fill = ..level..), geom = "polygon", alpha = .5) +
    geom_point(alpha = 0.6) +
    labs(
      title = paste(x, "vs", y),
      subtitle = "",
      x = x, y = y) +
    xlim(0, 1) +
    ylim(0, 1) +
    scale_color_brewer(palette = "Set1") +
    guides(colour = guide_legend(title = cluster)) +
    theme(legend.position="bottom")

  if(!marginal) return(p)

  require(ggExtra)
  ggMarginal(p, type="histogram", fill = "gainsboro", binwidth = 0.01, alpha = 1,
             aes = aes(
               x = eval(parse(text = x)),
               colour = eval(parse(text = cluster)),
               y = eval(parse(text = y))
               ),
             data = data,
             xparams = list(colour = "black", size = 0.1),
             yparams = list(colour = "black", size = 0.1))

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
plot_grid = function(data, samples) {

  require(ggplot2)

  #  Diagonal
  MB = data$best.MOBSTER
  plots = lapply(
    seq(MB),
    function(w) plot(MB[[w]], silent = TRUE, main = paste("MOBSTER ", samples[w]) )[[1]])

  for(s in seq(samples)) {
    for(w in s:length(samples)) {
      if(s != w) {

        pl = plot_2DVAF(
          data$data,
          x = paste0('VAF.', samples[s]),
          y = paste0('VAF.', samples[w]),
          cluster = 'sciClone.cluster')

        plots = append(plots, list(pl))
      }
    }
  }

  for(s in seq(samples)) {
    for(w in s:length(samples)) {
      if(s != w) {

        pl = plot_2DVAF(
          data$data,
          x = paste0('VAF.projected.', samples[s]),
          y = paste0('VAF.projected.', samples[w]),
          cluster = 'MOBSTER.sciClone.cluster')

        plots = append(plots, list(pl))
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

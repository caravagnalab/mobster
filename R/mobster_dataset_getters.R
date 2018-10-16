####################### Getters for the plain tibble with the stored data

#' Extract VAF values.
#' 
#' @description 
#' 
#' Extract from the internal representation of an object all the entries
#' that refer to the VAF values. These can be subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A \code{mbst_data} object
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries
#' @export
#'
#' @examples
VAF = function(x,
               ids = keys(x),
               samples = x$samples)
{
  x$data %>%
        filter(variable == 'VAF' &
                 id %in% ids &
                 sample %in% samples)
}

#' Extract DP values.
#' 
#' @description 
#' 
#' Extract from the internal representation of an object all the entries
#' that refer to the DP values. These can be subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A \code{mbst_data} object
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries
#' @export
#'
#' @examples
DP = function(x,
              ids = keys(x),
              samples = x$samples)
{
  x$data %>%
    filter(variable == 'DP' &
             id %in% ids &
             sample %in% samples)
}

#' Extract NV values.
#' 
#' @description 
#' 
#' Extract from the internal representation of an object all the entries
#' that refer to the NV values. These can be subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A \code{mbst_data} object
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries
#' @export
#'
#' @examples
NV = function(x,
              ids = keys(x),
              samples = x$samples)
{
  x$data %>%
    filter(variable == 'NV' &
             id %in% ids &
             sample %in% samples)
}

#' Extract annotation values.
#' 
#' @description 
#' 
#' Extract from the internal representation of an object all the entries
#' that refer to the annotation values. These are all values which are not
#' any of VAF, DP or NV, or any of the location coordinates for the input
#' mutations. As for the other getters, you can subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A \code{mbst_data} object
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A tibble of the required entries
#' @export
#'
#' @examples

Annotations = function(x,
                       ids = keys(x),
                       variables = unique(x$annotations$variable))
{
  x$annotations %>%
    filter(variable %in% variables &
             id %in% ids)
}

####################### Getters for the tibble, which we trasnform to a table with the stored data

#' Tabular VAF values.
#' 
#' @description 
#' 
#' Similarly to \code{VAF}, this function however returns a spread-like table 
#' with one column per VAF value. As \code{VAF} the entries can be subset as
#' required, and the columns named with a particular format.
#'
#' @param x A \code{mbst_data} object
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param suffix Every column will be called as \code{sample_name.VAF} when \code{suffix = '.VAF'}
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
VAF_table = function(x,
               ids = keys(x),
               samples = x$samples,
               suffix = '.VAF')
{
  output = tibble(`id` = ids)
  
  for(s in samples) {
    entry = VAF(x, ids = ids, samples = s) %>% spread(variable, value) %>% select(id, VAF) 
    # %>% rename(!!s := VAF)

    output = full_join(entry, output, by = 'id')  
  }
  
  colnames(output) = c('id', paste0(samples, suffix))
  
  output
}

#' Tabular NV values.
#' 
#' @description 
#' 
#' Similarly to \code{NV}, this function however returns a spread-like table 
#' with one column per VAF value. As \code{NV} the entries can be subset as
#' required, and the columns named with a particular format.
#'
#' @param x A \code{mbst_data} object
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param suffix Every column will be called as \code{sample_name.NV} when \code{suffix = '.NV'}
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
NV_table = function(x,
                     ids = keys(x),
                     samples = x$samples,
                     suffix = '.NV')
{
  output = tibble(`id` = ids)
  
  for(s in samples) {
    entry = NV(x, ids = ids, samples = s) %>% spread(variable, value) %>% select(id, NV) 
    # %>% rename(!!s := NV)
    
    output = full_join(entry, output, by = 'id')  
  }
  
  colnames(output) = c('id', paste0(samples, suffix))
  
  output
}


#' Tabular DP values.
#' 
#' @description 
#' 
#' Similarly to \code{DP}, this function however returns a spread-like table 
#' with one column per VAF value. As \code{DP} the entries can be subset as
#' required, and the columns named with a particular format.
#'
#' @param x A \code{mbst_data} object
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param suffix Every column will be called as \code{sample_name.DP} when \code{suffix = '.DP'}
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
DP_table = function(x,
                    ids = keys(x),
                    samples = x$samples,
                    suffix = '.DP')
{
  output = tibble(`id` = ids)
  
  for(s in samples) {
    entry = DP(x, ids = ids, samples = s) %>% spread(variable, value) %>% select(id, DP) 
    # %>% rename(!!s := NV)
    
    output = full_join(entry, output, by = 'id')  
  }
  
  colnames(output) = c('id', paste0(samples, suffix))
  
  output
}

#' Tabular VAF, DP and NV values.
#' 
#' @description 
#' 
#' This is just a wrapper to a combined call of function \code{VAF_table}, \code{DP_table} and 
#' \code{NV_table}. The resulting outputs are bound by column; ususal subset options are available.
#'
#' @param x A \code{mbst_data} object
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
Data_table = function(x,
                    ids = keys(x),
                    samples = x$samples)
{
  full_join(
    full_join(
      VAF_table(x, ids, samples),
      DP_table(x, ids, samples),
      by = 'id'),
    NV_table(x, ids, samples),
    by = 'id'
  )
}


####################### Getters for the cohort size, etc.

#' Return the number of mutations in the dataset.
#'
#' @param x A \code{mbst_data} object
#'
#' @return The number of mutations in the dataset.
#' @export
#'
#' @examples
N = function(x) {
  nrow(VAF(x, samples = x$samples[1]))
}

####################### Getters for the clustering computed


#' Getter for MOBSTER clustering
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
MClusters = function(x, annotations = FALSE)
{
  if (is.null(x$fit.MOBSTER))
    stop("MOBSTER clusters are not available!")
  
  list.best = lapply(x$fit.MOBSTER,
                     function(w)
                       return(w$best$data %>% select(-sample,-VAF)))
  
  MOBSTER_clusters = list.best[[1]]
  
  for (i in 2:length(list.best)) {
    MOBSTER_clusters = full_join(
      MOBSTER_clusters,
      list.best[[i]],
      by = 'id',
      suffix =
        c(paste0('.', names(x$fit.MOBSTER)[i - 1]),
          paste0('.', names(x$fit.MOBSTER)[i]))
    )
  }
  
  if (annotations)
  {
    annotations = Annotations(x, ids = MOBSTER_clusters$id) %>%
      spread(variable, value)
    
    MOBSTER_clusters = full_join(MOBSTER_clusters, annotations, by = 'id')
    MOBSTER_clusters = full_join(MOBSTER_clusters, x$map_mut_seg, by = 'id')
  }
  
  MOBSTER_clusters
}

#' Getter for sciClone clustering
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
SClusters = function(x, annotations = FALSE)
{
  if (is.null(x$fit.sciClone))
    stop("sciClone clusters are not available!")
  
  # list.best = lapply(
  #   x$fit.MOBSTER,
  #   function(w) return(w$best$data %>% select(-sample, -VAF))
  # )
  #
  # MOBSTER_clusters = list.best[[1]]
  #
  # for(i in 2:length(list.best)) {
  #   MOBSTER_clusters = full_join(MOBSTER_clusters,
  #                                list.best[[i]],
  #                                by = 'id',
  #                                suffix =
  #                                  c(
  #                                    paste0('.', names(x$fit.MOBSTER)[i - 1]),
  #                                    paste0('.', names(x$fit.MOBSTER)[i])
  #                                  ))
  # }
  #
  # if(annotations)
  # {
  #   annotations = Annotations(x, ids = MOBSTER_clusters$id) %>%
  #     spread(variable, value)
  #
  #   MOBSTER_clusters = full_join(MOBSTER_clusters, annotations, by = 'id')
  #   MOBSTER_clusters = full_join(MOBSTER_clusters, x$map_mut_seg, by = 'id')
  # }
  #
  # MOBSTER_clusters
}

#
# MOBSTER_clusters$anyTail =
#   apply(MOBSTER_clusters, 1, function(w) any(w == 'Tail', na.rm = TRUE) )
#





####################### Private getters
keys = function(x) {
  unique(x$data$id)
}

byLoc = function(x, loc.id) {
  muts_CN = x$locations %>%
    spread(variable, value)
  
  thisSeg = x$segments %>%
    filter(id == loc.id)
  
  which_muts = muts_CN %>%
    filter(chr == thisSeg$chr[1] &
             from >= thisSeg$from[1] &
             to <= thisSeg$to[1])
  
  which_muts
}

minor = function(x, seg_id, samples = x$samples)
{
  x$segments %>% filter(id == seg_id &
                          sample == samples &
                          variable == 'minor') %>% pull(value)
}

Major = function(x, seg_id, samples = x$samples)
{
  x$segments %>% filter(id == seg_id &
                          sample == samples &
                          variable == 'Major') %>% pull(value)
}

# Converter for sciCLone inputs
convert_sciClone_input = function(x)
{
  seg_ids = unique(x$segments$id)
  
  # CNA
  CN = NULL
  for (s in x$samples)
  {
    CN.copies = mobster:::minor(x, seg_id = seg_ids, samples = s) +
      mobster:::Major(x, seg_id = seg_ids, samples = s)
    
    CN.segment = x$segments %>%
      select(chr, from, to) %>%
      distinct()
    colnames(CN.segment) = c("chr", 'start', 'stop')
    
    CN.segment$chr = as.numeric(sapply(CN.segment$chr, function(w)
      substr(w, 4, nchar(w))))
    
    CN.segment$segment_mean = CN.copies
    
    CN = append(CN, list(as.data.frame(CN.segment)))
  }
  
  names(CN) = x$samples
  
  # MUTAITONS
  locs = x$locations %>% spread(variable, value)
  
  MUTS = NULL
  
  for (s in x$samples)
  {
    sample.data = bind_rows(VAF(x, samples = s),
                            DP(x, samples = s),
                            NV(x, samples = s))
    
    sample.data = sample.data %>% spread(variable, value)
    sample.data = sample.data[complete.cases(sample.data),]
    
    sample.data = sample.data %>% select(id, DP, NV, VAF)
    sample.data = sample.data %>%
      mutate(refCount = DP - NV,
             varCount = NV,
             vaf = VAF * 100) %>%
      select(id, refCount, varCount, vaf)
    
    sample.locs = locs %>%
      filter(id %in% sample.data$id) %>%
      mutate(start = from) %>%
      select(id, chr, start)
    sample.locs$start = as.numeric(sample.locs$start)
    
    df = full_join(sample.data, sample.locs, by = 'id')
    df = df %>% select(chr, start, refCount, varCount, vaf)
    
    df$chr = as.numeric(sapply(df$chr, function(w)
      substr(w, 4, nchar(w))))
    
    MUTS = append(MUTS,
                  list(as.data.frame(df)))
  }
  
  list(CN = CN, MUTS = MUTS)
}

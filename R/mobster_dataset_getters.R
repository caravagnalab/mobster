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

    output = full_join(output, entry, by = 'id')  
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
    
    output = full_join(output, entry, by = 'id')  
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
    
    output = full_join(output, entry, by = 'id')  
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

# K = 

####################### Getters for the clustering computed

#' Getter for clustering
#'
#' @param x MOBSTER dataset.
#' @param cluster Use one of MOBSTER (M), sciClone (S), pyClone (P), Binomial (B).
#' @param annotations TRUE or FALSE will return also the latent variables if possible.
#'
#' @return
#' @export
#'
#' @examples
Clusters = function(x, cluster, annotations = FALSE){
  
  switch(cluster,
         MOBSTER = return(MClusters(x, annotations)),
         M = return(MClusters(x, annotations)),
         sciClone = return(SClusters(x, annotations)),
         S = return(MClusters(x, annotations)),
         pyClone = return(PClusters(x, annotations)),
         P = return(PClusters(x, annotations)),
         Binomial = return(BClusters(x, annotations)),
         B = return(BClusters(x, annotations)))
  
  
  stop("Cluster code not recognized; use one of MOBSTER (M), sciClone (S), pyClone (P), Binomial (B)")
}

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
  
  MOBSTER_clusters = list.best %>% purrr::reduce(full_join, by = "id")
  colnames(MOBSTER_clusters)[2:ncol(MOBSTER_clusters)] = paste0('cluster.', names(x$fit.MOBSTER)) 
  
  # Reduce(
  #   function(w) {full_join},
  #   list.best
  # )
  
  # MOBSTER_clusters = list.best[[1]]
  # 
  # for (i in 2:length(list.best)) {
  #   MOBSTER_clusters = full_join(
  #     MOBSTER_clusters,
  #     list.best[[i]],
  #     by = 'id',
  #     suffix =
  #       c(paste0('.', names(x$fit.MOBSTER)[i - 1]),
  #         paste0('.', names(x$fit.MOBSTER)[i]))
  #   )
  # }
  
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
  
  sciClone.fit = x$fit.sciClone
  
  clusters = as_tibble(sciClone.fit@vafs.merged)
  
  clusters = clusters %>% select(NA1, cluster, starts_with('cluster'))
  colnames(clusters)[c(1,2)] = c('id', 'cluster.sciClone')
  
  if (annotations)
  {
    annotations = Annotations(x, ids = clusters$id) %>%
      spread(variable, value)
    
    clusters = full_join(clusters, annotations, by = 'id')
    clusters = full_join(clusters, x$map_mut_seg, by = 'id')
  }
  
  clusters
}


#' Getter for pyClone clustering
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
PClusters <- function(x, annotations = FALSE) {
  
  if (is.null(x$fit.pyClone)) 
    stop("pyClone clusters are not available!\n")
  
  fit = x$fit.pyClone
  
  clusters = unique(fit$loci[,c("mutation_id","cluster_id")])
  clusters = as_tibble(clusters)
  colnames(clusters)[c(1,2)] = c("id", "cluster.pyClone")
  
  clusters = clusters[match(mobster:::keys(x), clusters$id),]
  
  if (annotations)
  {
    annotations = Annotations(x, ids = clusters$id) %>%
      spread(variable, value)
    
    clusters = full_join(clusters, annotations, by = 'id')
    clusters = full_join(clusters, x$map_mut_seg, by = 'id')
  }
  
  return(clusters)
}

#' Getter for Binomial clustering
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
BClusters <- function(x, annotations = FALSE) {
  
  if (is.null(x$fit.Binomial)) 
    stop("Binomial clusters are not available!\n")
  
  clusters = x$fit.Binomial$X
  
  if (annotations)
  {
    annotations = Annotations(x, ids = clusters$id) %>%
      spread(variable, value)
    
    clusters = full_join(clusters, annotations, by = 'id')
    clusters = full_join(clusters, x$map_mut_seg, by = 'id')
  }
  
  clusters
}

#
# MOBSTER_clusters$anyTail =
#   apply(MOBSTER_clusters, 1, function(w) any(w == 'Tail', na.rm = TRUE) )
#





####################### Private getters
keys = function(x) {
  unique(x$data$id)
}


split_segments_around_centromers = function()
{
  data("cytoband_hg19")
  
  new.segments = NULL
  
  for(i in 1:nrow(x$segments)) {
    
    seg = x$segments[i, ]
    
    SEG.cytoband = cytoband[cytoband$chr == seg$chr, ]
    SEG.cytoband = SEG.cytoband[SEG.cytoband$region == 'acen', ]
    
    ends_before = seg$to < SEG.cytoband$from[1]
    starts_after = seg$from > SEG.cytoband$to[2]
    
    # 
    if(ends_before | starts_after) new.segments = rbind(new.segments, seg)
    else {
      
      
    }
    
    
  }
  
  
}

byLoc = function(x, loc.id, offset_around_centromers) {
  muts_CN = x$locations %>%
    spread(variable, value)
  
  thisSeg = x$segments %>%
    filter(id == loc.id)
  
  # Enforce the correct types!
  muts_CN$chr = as.character(muts_CN$chr)
  muts_CN$from = as.numeric(muts_CN$from)
  muts_CN$to = as.numeric(muts_CN$to)
  
  thisSeg$chr = as.character(thisSeg$chr)
  thisSeg$from = as.numeric(thisSeg$from)
  thisSeg$to = as.numeric(thisSeg$to)
  
  # Within segment
  which_muts = muts_CN %>%
    filter(chr == thisSeg$chr[1] &
             from >= thisSeg$from[1] &
             to <= thisSeg$to[1])
  
  ws = nrow(which_muts)

  # thisSeg
  pio::pioStr("Chromosome", thisSeg$chr[1], suffix = '')
  pio::pioStr(" from", thisSeg$from[1], suffix = '')
  pio::pioStr(" to", thisSeg$to[1], suffix = '')
  pio::pioStr(" : n =", ws, suffix = '')
  
  # load centromers data -- these are in absolute location format, so we update the from/ to locations
  if(offset_around_centromers > 0) {
    
    # data("cytoband_hg19")
    # 
    # SEG.cytoband = cytoband[cytoband$chr == thisSeg$chr[1], ]
    # SEG.cytoband = SEG.cytoband[SEG.cytoband$region == 'acen', ]
    # 
    # which_muts = which_muts %>%
    #   filter(from < SEG.cytoband$from[1] |
    #            to > SEG.cytoband$to[2])
   
    # Get coordinates and offset the centromers by "offset"
    data('chr_coordinate_hg19', package = 'mobster')
    
    chr_coordinate_hg19 = chr_coordinate_hg19 %>% 
      mutate(
        centromerStart = centromerStart - offset_around_centromers,
        centromerEnd = centromerEnd + offset_around_centromers
        ) %>% 
      filter(chr == thisSeg$chr[1])
    
    which_muts = which_muts %>%
      filter(to < chr_coordinate_hg19$centromerStart[1] |
               from > chr_coordinate_hg19$centromerEnd[1])
    
    ws_nc = nrow(which_muts)
    pio::pioStr("off centromer", ws_nc, suffix = '\n')
  }
  
  which_muts
}

minor = function(x, seg_id, samples = x$samples)
{
  as.numeric(
    x$segments %>% filter(id == seg_id &
                          sample %in% samples &
                          variable == 'minor') %>% pull(value)
    )
}

Major = function(x, seg_id, samples = x$samples)
{
  as.numeric(
    x$segments %>% filter(id == seg_id &
                          sample %in% samples &
                          variable == 'Major') %>% pull(value)
  )
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
    df = df %>% select(chr, start, refCount, varCount, vaf, id)
    
    df$chr = as.numeric(sapply(df$chr, function(w)
      substr(w, 4, nchar(w))))
    
    MUTS = append(MUTS,
                  list(as.data.frame(df)))
  }
  
  list(CN = CN, MUTS = MUTS)
}


# # Phylogenetic tree REVOLVER
# as_revolver_dataset = function(x, 
#                                clonal,
#                                drivers,
#                                patientID = "MyPatient",
#                                variantID = 'id',
#                                remove_clusters = NULL
# ) {
#   
#   # # VAF  + Annotations 
#   # vaf = VAF_table(x)
#   # 
#   # matched = full_join(vaf, annotations, by = 'id')  
#   # 
#   
#   # Clusters
#   matched = BClusters(x) %>% select(id, cluster.Binomial)
#   matched = matched %>% rename(cluster = cluster.Binomial)
#   if(!is.null(remove_clusters)) matched = matched %>% filter(!(cluster %in% remove_clusters))
#   
#   
#   matched$Misc = ""
#   
#   # patientID
#   matched$patientID = patientID
#   
#   # drivers and clonal
#   annotations = Annotations(x) %>% spread(variable, value) %>% 
#     filter(gene %in% drivers) %>% filter(type == 'exonic') %>% pull(id)
#   
#   matched = matched %>% mutate(is.clonal = cluster == clonal)
#   matched = matched %>% mutate(is.driver = id %in% annotations)
#   
#   # CCF
#   vaf = VAF(x) %>% mutate(entry = paste0(sample, ':', value))
#   vaf = vaf %>% group_by(id) %>% summarize(CCF = paste(entry, collapse = ';'))
#   
#   matched = matched %>% left_join(vaf, by = 'id')  
#   
#   # variant key
#   matched = matched %>% rename(variantID = !!variantID)
#   
#   matched
# }
# 
# fitgp_dataset = as_revolver_dataset(fitgp, clonal = 'C1', drivers, patientID = "AfterMOB_all")
# fitgp_prv_dataset = as_revolver_dataset(fitgp, clonal = 'C1', drivers, patientID = "PRV_AfterMOB_all", remove_clusters = prioritize_Clusters(fitgp))
# 
# fit_dataset = as_revolver_dataset(fit, clonal = 'C2', drivers, patientID = "NoMOB_all")
# fit_prv_dataset = as_revolver_dataset(fit, clonal = 'C2', drivers, patientID = "PRV_NoMOB_all", remove_clusters = prioritize_Clusters(fit))
# 
# all_datasets = bind_rows(fitgp_dataset, fitgp_prv_dataset, fit_dataset, fit_prv_dataset)
# 
# library(revolver)
# cohort = revolver_cohort(as.data.frame(all_datasets))
# 
# cohort = revolver_compute_phylogenies(cohort, 'AfterMOB_all')
# cohort = revolver_compute_phylogenies(cohort, 'NoMOB_all')
# cohort = revolver_compute_phylogenies(cohort, 'PRV_AfterMOB_all')
# cohort = revolver_compute_phylogenies(cohort, 'PRV_NoMOB_all')
# 
# 
# revolver_report_patient(cohort, 'MyPatient')


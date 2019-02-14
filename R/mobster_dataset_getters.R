#' Extract keys to id mutations loaded inside a MOBSTER dataset.
#'
#' @param x A MOBSTER \code{mbst_data} object
#'
#' @return The keys that id mutations.
#' @export
#'
#' @examples
keys = function(x) {
  unique(x$data$id)
}

#' Extract VAF values from a MOBSTER dataset.
#' 
#' @description 
#' 
#' Extract from the internal representation of an object all the entries
#' that refer to the VAF values. These can be subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A MOBSTER \code{mbst_data} object
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

#' Extract depth values from a MOBSTER dataset.
#' 
#' @description 
#' 
#' Extract from the internal representation of an object all the entries
#' that refer to the DP values. These can be subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A MOBSTER \code{mbst_data} object.
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

#' Extract the number of reads with the alternative allele from a MOBSTER dataset.
#' 
#' @description 
#' 
#' Extract from the internal representation of an object all the entries
#' that refer to the NV values. These can be subset by sample and mutation
#' ID; by default all entries are returbed. The output is a tibble; no other
#' transformations are executed.
#'
#' @param x A MOBSTER \code{mbst_data} object.
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
#' @param x A MOBSTER \code{mbst_data} object.
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

#' Tabular VAF values.
#' 
#' @description 
#' 
#' Similarly to \code{VAF}, this function however returns a spread-like table 
#' with one column per VAF value. As \code{VAF} the entries can be subset as
#' required, and the columns named appending a custom suffix to sample names.
#'
#' @param x A MOBSTER \code{mbst_data} object.
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
#' required, and the columns named appending a custom suffix to sample names.
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
#' required, and the columns named appending a custom suffix to sample names.
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

    output = full_join(output, entry, by = 'id')  
  }
  
  colnames(output) = c('id', paste0(samples, suffix))
  
  output
}


#' Tabular VAF, DP and NV values, with annotations.
#' 
#' @description 
#' 
#' This is just a wrapper to a combined call of function \code{VAF_table}, \code{DP_table} and 
#' \code{NV_table}. The resulting outputs are bound by column; ususal subset options are available.
#' The output can be augmented with annotations from each available mutation.
#'
#' @param x A \code{mbst_data} object
#' @param ids The IDs of the mutations to extract, all by default.
#' @param samples The samples for which we extract the data, all by default.
#' @param annotations False by default; if true annotations are also returned.
#'
#' @return A spread tibble of the required entries.
#' @export
#'
#' @examples
Data_table = function(x,
                    ids = keys(x),
                    samples = x$samples,
                    annotations = FALSE)
{
  # Joined data tables
  dt = full_join(
    full_join(
      VAF_table(x, ids, samples),
      DP_table(x, ids, samples),
      by = 'id'),
    NV_table(x, ids, samples),
    by = 'id'
  )
  
  # All variables
  an = Annotations(x, ids)
  
  full_join(dt, an, by = 'id')
}


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








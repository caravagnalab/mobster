#' Title
#'
#' @param x 
#' @param ids 
#' @param samples 
#'
#' @return
#' @export
#'
#' @examples
VAF = function(x, ids = keys(x), samples = x$samples)
{
  x$data %>%
    filter(variable == 'VAF' & 
             id %in% ids & 
             sample %in% samples)
}

#' Title
#'
#' @param x 
#' @param ids 
#' @param samples 
#'
#' @return
#' @export
#'
#' @examples
DP = function(x, ids = keys(x), samples = x$samples)
{
  x$data %>%
    filter(variable == 'DP' & 
             id %in% ids &  
             sample %in% samples)
}

#' Title
#'
#' @param x 
#' @param ids 
#' @param samples 
#'
#' @return
#' @export
#'
#' @examples
NV = function(x, ids = keys(x), samples = x$samples)
{
  x$data %>%
    filter(variable == 'NV' & 
             id %in% ids & 
             sample %in% samples)
}

#' Title
#'
#' @param x 
#' @param ids 
#' @param variable 
#' @param samples 
#'
#' @return
#' @export
#'
#' @examples
Annotations = function(x, ids = keys(x), 
                       variables = unique(x$annotations$variable))
{
  
  x$annotations %>%
    filter(variable %in% variables &  
             id %in% ids)
}


#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
N = function(x){
  nrow(VAF(x, samples = x$samples[1]))
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
  if(is.null(x$fit.MOBSTER)) stop("MOBSTER clusters are not available!")
  
  list.best = lapply(
    x$fit.MOBSTER, 
    function(w) return(w$best$data %>% select(-sample, -VAF))
    )
  
  MOBSTER_clusters = list.best[[1]]
  
  for(i in 2:length(list.best)) {
    MOBSTER_clusters = full_join(MOBSTER_clusters, 
                                 list.best[[i]], 
                                 by = 'id', 
                                 suffix = 
                                   c(
                                     paste0('.', names(x$fit.MOBSTER)[i - 1]),
                                     paste0('.', names(x$fit.MOBSTER)[i])
                                   ))
  }
  
  if(annotations) 
  {
    annotations = Annotations(x, ids = MOBSTER_clusters$id) %>%
      spread(variable, value)
    
    MOBSTER_clusters = full_join(MOBSTER_clusters, annotations, by = 'id')
    MOBSTER_clusters = full_join(MOBSTER_clusters, x$map_mut_seg, by = 'id')
  }
  
  MOBSTER_clusters
}

# 
# MOBSTER_clusters$anyTail = 
#   apply(MOBSTER_clusters, 1, function(w) any(w == 'Tail', na.rm = TRUE) )
# 





####################### Private getters
keys = function(x) { unique(x$data$id) }

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
  for(s in x$samples)
  {

    CN.copies = mobster:::minor(x, seg_id = seg_ids, samples = s) +
      mobster:::Major(x, seg_id = seg_ids, samples = s)
      
    CN.segment = x$segments %>% 
      select(chr, from, to) %>%
      distinct() 
    colnames(CN.segment) = c("chr", 'start', 'stop')
    
    CN.segment$segment_mean = CN.copies
      
    CN = append(CN, list(as.data.frame(CN.segment)))
  }
  
  names(CN) = x$samples

  # MUTAITONS
  locs = x$locations %>% spread(variable, value)
  
  MUTS = NULL
  
  for(s in x$samples)
  {
    sample.data = bind_rows(
      VAF(x, samples = s),
      DP(x, samples = s),
      NV(x, samples = s))
    
    sample.data = sample.data %>% spread(variable, value)
    sample.data = sample.data[complete.cases(sample.data), ]
    
    sample.data = sample.data %>% select(id, DP, NV, VAF)
    sample.data = sample.data %>%
      mutate(refCount = DP - NV,
             varCount = NV,
             VAF = VAF * 100) %>%
      select(id, refCount, varCount, VAF)
  
    sample.locs = locs %>%  
      filter(id %in% sample.data$id) %>%
      mutate(start = from) %>%
      select(id, chr, start)
    
    
    MUTS = append(MUTS, 
                  list(
                    as.data.frame(
                      full_join(sample.data, sample.locs, by = 'id')
                    )))
  }
  
  list(CN = CN, MUTS = MUTS)
}


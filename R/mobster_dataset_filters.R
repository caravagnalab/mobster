#' Title
#'
#' @param x 
#' @param type 
#'
#' @return
#' @export
#'
#' @examples
mobster_flt_projection = function(x, type = 'global')
{
  clusters = MClusters(x) 
  
  colnames(clusters) = gsub('cluster.', '', colnames(clusters))
  
  cols.clusters = 2:ncol(clusters)

  # Return object
  y = x
  
  if(type == 'global')
  {
    clusters$projected = apply(
      clusters[, cols.clusters, drop = FALSE], 
      1, 
      function(w) any(w == 'Tail', na.rm = TRUE))
    
    pio::pioTit("The following will be projected with a global strategy")
    to_cancel = clusters %>% filter(projected)
    
    print(to_cancel)
    
    y = mobster:::delete_entries(x, to_cancel %>% pull(id))
    
    y = logOp(y, paste0("Projected read counts; type = ", type, ""))
    
  }
  
  if(type == 'local') stop("TODO")
    
  
  return(y)
}


#'Subsample data
#'
#' @param obj 
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
mobster_flt_downsample = function(obj, n)
{
  if(n < N(obj))
  {
    pio::pioTit(paste0("Subsampling ", n, "entries out of ", N(obj)))
    
    k = mobster:::keys(obj)
    ids = sample(k, length(k) - n)
    
    obj = mobster:::delete_entries(obj, ids)
    obj = logOp(obj, paste0("Subsampled to N = ", n, ""))
  }
  
  obj
}


# filter for VAF below `cutoff` parameters. This sets to 0 everything both the VAF
# and the NV entry, bute retains the coverage at the locus (DP). This function accepts
# also a purity vector in input to adjust accordingly the cutoff
#' Title
#'
#' @param x 
#' @param cutoff 
#'
#' @return
#' @export
#'
#' @examples
mobster_flt_minfreq = function(x, cutoff)
{
  
  # Adjusted cutoff
  ad.cutoff = cutoff/x$purity
  
  pio::pioTit("Cutoff(s) adjusted for purity")
  pio::pioDisp(ad.cutoff)
  
  for(s in x$samples)
  {
    ids = VAF(x, samples = s) %>% 
      filter(value < ad.cutoff[s] & value > 0) %>% pull(id)
    
    pio::pioStr(
      paste0("Sample ", s, ' - Number of entries below cutoff is N ='),
      length(ids), prefix = '\n')
    
    x$data = x$data %>% mutate_cond(
      id %in% ids & sample == s & variable %in% c("VAF", "NV"), value = 0)
  }
  
  # Log update
  x = logOp(x, paste0("VAF entries below (adjusted) cutoff ", cutoff, ' have been removed.'))
  
  all_zeroes(x)
}

# filter all mutations with VAF outside "VAF.range". Notice that the lower cutoff is 
# useless beacuse already we scaled for purity and removed some SNVs with the previous
# function. However, we can eliminate weird mutations which are in LOH regions wrongly
# called as diploid
#' Title
#'
#' @param x 
#' @param min.vaf 
#' @param max.vaf 
#'
#' @return
#' @export
#'
#' @examples
mobster_flt_vafrange = function(x, min.vaf, max.vaf)
{
  pio::pioTit(
    paste0("Keeping only mutations with VAF in (", min.vaf, ' - ', max.vaf, ') in all samples')
  )
  
  ids = VAF(x) %>% 
    filter(value > 0 & 
             (value < min.vaf | value > max.vaf)) %>% 
    select(id) %>%
    distinct() %>% pull(id)
  
  pio::pioStr(
    paste0("Numner of entries outside range is N ="),
    length(ids), prefix = '\n', suffix = '\n')
  
  x = delete_entries(x, ids)
  
  # Log update
  x = logOp(x, paste0("Entries with VAF outside range ", min.vaf, '-', max.vaf, ' removed.'))
  
  all_zeroes(x)
}

# Ad above but for DP
#' Title
#'
#' @param x 
#' @param min.DP 
#' @param max.DP 
#'
#' @return
#' @export
#'
#' @examples
mobster_flt_dprange = function(x, min.DP, max.DP)
{
  pio::pioTit(
    paste0("Cutting any mutation with DP outside range", min.DP, '-', max.DP)
  )
  
  ids = DP(x) %>% 
    filter(value > 0 & 
             (value < min.DP | value > max.DP)) %>% 
    select(id) %>%
    distinct() %>% pull(id)
  
  pio::pioStr(
    paste0("Number of entries outside range is N ="),
    length(ids), prefix = '\n')
  
  x = delete_entries(x, ids)
  
  x = logOp(x, paste0("Entries with DP outside range ", min.DP, '-', max.DP, ' removed.'))
  
  all_zeroes(x)
}




# Private functions to subset the data
delete_entries = function(x, ids)
{
  x$data = x$data %>% 
    filter(!id %in% ids)

  if(!is.null(x$annotations))
    x$map_mut_seg = x$map_mut_seg %>% 
      filter(!id %in% ids)
  
  if(!is.null(x$map_mut_seg))
    x$annotations = x$annotations %>% 
      filter(!id %in% ids)
  
  if(!is.null(x$VAF_cn_adjustment))
    x$VAF_cn_adjustment = x$VAF_cn_adjustment %>% 
      filter(!id %in% ids)
  
  if(!is.null(x$locations))
    x$locations = x$locations %>% 
      filter(!id %in% ids)
  
  x
}


all_zeroes = function(x)
{
  ids = VAF(x) %>%
    filter(value == 0) %>%
    group_by(id) %>%
    summarise(nsamples = length(unique(sample))) %>%
    filter(nsamples == length(x$samples)) %>%
    pull(id)
  
  if(length(ids) > 0) {
    
    message("\nN = ", length(ids), ' entries have VAF = 0 in all samples, will be removed.')
  
  }
  
  delete_entries(x, ids)
}

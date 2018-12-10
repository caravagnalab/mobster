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
  # get clusters
  clusters = MClusters(x) 
  
  colnames(clusters) = gsub('cluster.', '', colnames(clusters))
  
  cols.clusters = 2:ncol(clusters)

  if(type == 'global')
  {
    clusters$projected = apply(
      clusters[, cols.clusters, drop = FALSE], 
      1, 
      function(w) any(w == 'Tail', na.rm = TRUE))
    
    pio::pioTit("Global projection strategy")
    
    to_cancel = clusters %>% filter(projected)

    pio::pioStr("Entries to project", nrow(to_cancel), suffix = '\n')
    print(to_cancel)
    
    x = mobster:::delete_entries(x, to_cancel %>% pull(id))
    
    x = logOp(x, paste0("Projected read counts; type = ", type, ""))
    
  }
  
  if(type == 'local') {
    
    clusters = clusters %>% reshape2::melt(id = 'id') %>% as_tibble
    clusters$variable = gsub('cluster.', '', clusters$variable)
    
    clusters = clusters %>% filter(value == 'Tail')
    # clusters %>% group_by(variable) %>% summarize(n = n())
    
    pio::pioStr("Entries to project", nrow(clusters), suffix = '\n')
    print(clusters)
    
    pio::pioStr("Distributed as ...", '\n')
    print(clusters %>% group_by(variable) %>% summarize(n = n()))
    
    for(i in 1:nrow(clusters))
    {
      entry = clusters[i, ]
      
      x$data = x$data %>% 
        mobster:::mutate_cond(
          id == entry$id & sample == entry$variable & variable == "VAF", 
          value = 0)
      
      x$data = x$data %>% 
        mobster:::mutate_cond(
          id == entry$id & sample == entry$variable & variable == "NV", 
          value = 0)
    }
    
    x = all_zeroes(x)
    
    x = logOp(x, paste0("Projected read counts; type = ", type, ""))
  }
    
  return(x)
}


#'Subsample data
#'
#' @param x 
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
mobster_flt_downsample = function(x, n)
{
  if(n < N(x))
  {
    pio::pioTit(paste0("Subsampling ", n, "entries out of ", N(x)))
    
    k = mobster:::keys(x)
    ids = sample(k, length(k) - n)
    
    x = mobster:::delete_entries(x, ids)
    x = logOp(x, paste0("Subsampled to N = ", n, ""))
  }
  
  all_zeroes(x)
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
mobster_flt_minvaf = function(x, cutoff)
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


# Samve as above but with NV
#'
#' @param x 
#' @param cutoff 
#'
#' @return
#' @export
#'
#' @examples
mobster_flt_minnv = function(x, cutoff)
{
  for(s in x$samples)
  {
    ids = NV(x, samples = s) %>% 
      filter(value < cutoff & value > 0) %>% pull(id)
    
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


# UNUSED!
break_around_centromers = function(segments, offset = 1e6, relative = TRUE)
{
  # Get coordinates and offset the centromers by "offset"
  data('chr_coordinate_hg19', package = 'mobster')
  
  chr_coordinate_hg19 = chr_coordinate_hg19 %>% 
    mutate(
      centromerStart = centromerStart - offset,
      centromerEnd = centromerEnd + offset)
  
  # Work out where the segments map, first rescale them if relative
  if(relative)
  {
    segments = segments %>% left_join(chr_coordinate_hg19, by = 'chr')
    segments = segments %>% mutate(from = from.x + from.y, to = from.y + to.x)
    # segments = segments %>% select(-from.x, -to.x, -from.y, -to.y, -centromerStart, -centromerEnd)
  }
  
  
  # chopper
  segment_chop = function(a, b, x, y) {
    a = as.numeric(a)
    b = as.numeric(b)
    x = as.numeric(x)
    y = as.numeric(y)
    if(a >= y || b <= x) return(data.frame(from = a, to = b, stringsAsFactors = FALSE))
    if(b >= x & b <= y) return(data.frame(from = a, to = x, stringsAsFactors = FALSE))
    if(a >= x & a <= y) return(data.frame(from = y, to = b, stringsAsFactors = FALSE))
    if(a <= x & b >= y) return(data.frame(from = c(a, y), to = c(x, b), stringsAsFactors = FALSE))
    if(a >= x & b <= y) return(NULL)
    stop("error")
  }
  
  # segments$from = 
  
  newsegments = apply(segments, 1, function(w) {
    
    print(w)
    
    seg = segment_chop(a = w['from'], b = w['to'], x = w['centromerStart'], y = w['centromerEnd'])
    if(is.null(seg)) return(NULL)
    
    replicates = do.call("rbind", replicate(nrow(seg), w, simplify = FALSE))
    replicates = as_tibble(replicates)
    
    replicates = replicates %>% select(-from.x, -to.x, -from.y, -to.y, -centromerStart, -centromerEnd, -length)
    replicates = bind_cols(replicates, seg)
    
    replicates %>% mutate(length = to - from)
  })
  
  newsegments = Reduce(rbind, newsegments)
  rownames(newsegments) = NULL

  newsegments
}

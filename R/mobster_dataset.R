#' Title
#'
#' @param data 
#' @param samples 
#' @param segments 
#' @param purity 
#' @param N.min 
#' @param na.rm 
#' @param description 
#'
#' @return
#' @export
#'
#' @examples
mobster_dataset = function(
  data,
  samples,
  description = "My own dataset",
  segments,
  purity,
  N.min = 500,
  na.rm = TRUE
)
{
  pio::pioHdr("MOBSTER dataset")
  
  if(!is.data.frame(data)) stop("Data must be a dataframe")
  
  pio::pioStr('Number of mutations (N) =', nrow(data))
  
  # Columns with data
  DP.columns = paste0(samples, ".DP")
  NV.columns = paste0(samples, ".NV")
  VAF.columns = paste0(samples, ".VAF")
  
  all.columns = c(DP.columns, NV.columns, VAF.columns)
  
  if(any(!(DP.columns %in% colnames(data)))) stop("Missing DP columns in data.")
  if(any(!(NV.columns %in% colnames(data)))) stop("Missing NV columns in data.")
  if(any(!(VAF.columns %in% colnames(data)))) stop("Missing VAF columns in data.")
  
  # mutations must have valid locations
  if(any(!(c('chr', 'from', 'to') %in% colnames(data)))) stop("The position of the annotated mutations is missing.")

  # NAs in any of the columns for values
  if(any(is.na(data[, all.columns, drop = FALSE])))
  {
    data = data[
      complete.cases(data[, all.columns, drop = FALSE]), , drop = FALSE 
    ]
    pio::pioStr('Removed NAs, N =', nrow(data))
  }
  
  # ids
  if('id' %in% colnames(data)) {
    stop("Data contains an $id column, remove it.")
  }
  else {
    data$id = paste0('__mut', 1:nrow(data))
  }
  
  # Prepare a tibble for the input data 
  tib_data = as_tibble(data)
  
  pio::pioTit("Creating tibble for input data")
  
  tib_data = tib_data %>% 
    reshape2::melt(id = 'id') %>%
    as_tibble

  # make everything a chr
  tib_data$variable = paste(tib_data$variable)

  ##=============================##
  # Create a mbst_data object  #
  ##=============================##
  x = list()
  class(x) <- "mbst_data"
  
  x$samples = samples

  x$description = description
  
  # Annotations: anything but VAF/ CN values
  x$annotations = tib_data %>%
    filter(!(variable %in% 
               c(all.columns, 'chr', 'from', 'to')))
  
  # Data: DP, NV and VAF
  x$data = tib_data %>%
    filter(variable %in% all.columns)
  
  x$data$value = as.numeric(x$data$value)
  
  sp = strsplit(x$data$variable, '\\.')
  x$data$variable = sapply(sp, function(w) w[2])
  x$data$sample = sapply(sp, function(w) w[1])
  
  # Locations: CN position of mutations
  x$locations = tib_data %>%
    filter(variable %in%  c('chr', 'from', 'to'))
  
  # Log creation
  x = logOp(x, "Initialization")
  
  ########## Mapping mutations to CNAs
  Major.columns = paste0(x$samples, ".Major")
  minor.columns = paste0(x$samples, ".minor")
  
  if(any(!(names(purity) %in% x$samples))) stop("Missing purity for some samples?")
  
  # CN
  if(
    !is.data.frame(segments) |
    !all(c('chr', 'from', 'to') %in% colnames(segments)) |
    !all(Major.columns %in% colnames(segments)) |
    !all(minor.columns %in% colnames(segments))
  ) {
    stop("Your segments do not look in the right format.")
  }
  
  # get tibble for the segments
  tib_segmnets = as_tibble(segments)
  tib_segmnets$id = paste0("_CN_seg", 1:nrow(tib_segmnets))
  
  tib_segmnets = tib_segmnets %>%
    reshape2::melt(id = c('id', 'chr', 'from', 'to')) %>%
    as_tibble
  
  tib_segmnets$variable = paste(tib_segmnets$variable)
  
  sp = strsplit(tib_segmnets$variable, '\\.')
  tib_segmnets$variable = sapply(sp, function(w) w[2])
  tib_segmnets$sample = sapply(sp, function(w) w[1])
  
  # store segments in the obj
  x$segments = tib_segmnets  
  
  # store purity
  x$purity = purity
  
  # check that all samples have segments available
  if(!all(x$samples %in% unique(x$segments$sample))) stop("Some samples are missing from the list of segments")
  
  if(any(!(unique(x$segments$sample) %in% x$samples))) 
  {
    message("Your list of segments contains more sample IDs that those used -- will be dropped.")
    
    x$segments = x$segments %>%
      filter(sample %in% x$samples)
  }
  
  # we create a map for each mutation to each segment
  pio::pioTit("Mapping mutations to segments")
  pb = txtProgressBar(min = 1, max = nrow(segments), style = 3) 
  segments_ids = unique(x$segments$id)
  
  x$map_mut_seg = NULL
  
  for(s in seq(segments_ids)) {
    setTxtProgressBar(pb, s)
    
    mapped = mobster:::byLoc(x, segments_ids[s])
    if(nrow(mapped) == 0) next;
    
    mapped$seg_id = segments_ids[s]
    
    x$map_mut_seg = bind_rows(x$map_mut_seg, mapped)
  }
  
  if(nrow(x$map_mut_seg) == 0) {
    stop("None of the input mutations map to a segment, cannot do anything.")
  }

  # These have not been map to any segment 
  unmapped = VAF(x) %>%
    filter(!id %in% x$map_mut_seg$id) %>% pull(id)
  unmapped = unique(unmapped)
  
  if(length(unmapped) > 0)
  {
    num = length(unmapped)
    pnum = num/N(x) * 100
    
    pio::pioStr(
      'Mutations outside the input segments', 
      paste0('N = ', num),
      suffix = paste0('(', round(pnum, 0), '%)')
    )
    
    x = mobster:::delete_entries(x, unmapped)
  }
    
  # A summary table with the number of mutations per segment
  size_table = x$map_mut_seg %>%
    group_by(seg_id) %>%
    summarise(N = length(unique(id)))
  
  # we reject those with N < N.min
  rejected = size_table %>%
    filter(N < N.min) %>%
    arrange(desc(N)) %>%
    inner_join(y = x$segments, by=c("seg_id" = "id")) %>%
    select(seg_id, N, chr, from, to) %>%
    distinct
  
  accepted = size_table %>%
    filter(N >= N.min) %>%
    arrange(desc(N)) %>%
    inner_join(y = x$segments, by=c("seg_id" = "id")) %>%
    select(seg_id, N, chr, from, to) %>%
    distinct
  
  pio::pioTit(paste0("Discarded ", nrow(rejected)," segments with N < ", N.min))
  print(rejected)
  
  pio::pioTit(paste0("Will use ", nrow(accepted)," segments with N >= ", N.min))
  print(accepted)
  
  if(nrow(accepted) == 0)
    stop("No segments can be used, aborting.")
  
  # Subset to match accepted ... Non-accepted mutations can be immediately removed
  rejected = x$map_mut_seg %>%
    filter(seg_id %in% rejected$seg_id) %>% 
    pull(id)
  
  # N(x)
  
  if(length(rejected) > 0)
    x = mobster:::delete_entries(x, unique(rejected))
  
  pio::pioTit(paste0("Adjusting the VAF for ",
                     nrow(x$map_mut_seg)," entries across ", length(x$samples), " samples"))
  
  new.VAF = NULL
  
  # for every segment to map
  seg_to_match = unique(x$map_mut_seg$seg_id)
  pb = txtProgressBar(min = 1, max = length(seg_to_match), style = 3) 
  
  for(seg in seq(seg_to_match))  
  {
    setTxtProgressBar(pb, seg)
    
    # for every sample 
    for(s in x$samples) 
    {
      # get VAF values for these entries
      values = VAF(x,
                   ids = x$map_mut_seg %>% filter(seg_id == seg_to_match[seg]) %>% pull(id), 
                   samples = s
      ) %>% # Adjust VAF according to minor/ Major (CN=m+M)
        mutate(
          Major = mobster:::Major(x, seg_to_match[seg], s),
          minor = mobster:::minor(x, seg_to_match[seg], s),
          CN = Major + minor,
          purity = purity[s],
          adj_VAF = mobster:::vaf_adjustCN(value, minor, Major, purity)
        )
      
      new.VAF = bind_rows(new.VAF, values)
    }
  }
  
  # Save a copy of the mapping
  x$VAF_cn_adjustment = new.VAF
  
  new.VAF$value = new.VAF$adj_VAF  
  new.VAF = new.VAF %>% select(id, variable, value, sample)
  
  x$data = bind_rows(
    new.VAF,
    DP(x),
    NV(x)
  )
  
  # As CN, we can keep only those accepted
  x$segments = x$segments %>%
    filter(id %in% accepted$seg_id)
  
  # Log update
  x = logOp(x, "Added CN and adjusted VAF")
  
  all_zeroes(x)
  
  return(x)
}

#' Title
#'
#' @param x 
#' @param show.log 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
print.mbst_data = function(x, 
                           show.log = TRUE, 
                           cuts.DP = c(50, 100),
                           cuts.NV = c(5, 15),
                           cuts.VAF = c(0.05, 0.15),
                           ...)
{
  
  pio::pioHdr(paste("MOBSTER dataset - ", x$description))
  
  pio::pioStr("Samples", paste(x$samples, collapse = ', '), prefix = '\t- ', suffix = '\n')
  pio::pioStr("   N  =", N(x), prefix = '\t- ', suffix = '\n')
   
  prt = x$data %>%
    group_by(sample, variable) %>%
    filter(value > 0) %>%
    summarize(mean = mean(value), median = median(value), min = min(value), max = max(value)) %>%
    arrange(variable)
  
  # prt = prt %>%
  #   nest() %>%
  #   unlist(recursive = FALSE)
  # 
  # prt
  # 
  
  # cat(
  #   sprintf(
  #     ''
  #     )
  # )
  # 
  # for(i in 1:nrow(prt))
  # {
  #   cat(prt$sample[i], "")
  #       
  # }
  
  # prt = function(x, 
  #                cols,
  #                v)
  # {
  #   d = x$data[, cols]
  #   f = apply(d, 2, summary)
  #   f = apply(f, 2, round, digit = 2)
  #   
  #   colnames(f) = x$samples
  #   
  #   f = f[c(1,3,4,6), , drop = FALSE]
  # 
  #   frmt = function(f, v)
  #   {
  #     if(f < v[1]) cat(crayon::red(sprintf('%8s', f)))
  #     if(f >= v[1] & f <= v[2]) cat(crayon::yellow(sprintf('%8s', f)))
  #     if(f > v[2]) cat(crayon::green(sprintf('%8s', f)))
  # 
  #   }
  #   
  #   for(i in 1:nrow(f)) {
  #     cat(sprintf('  %10s', rownames(f)[i]))
  #     for(j in 1:ncol(f)) {
  #       frmt(f[i, j], v)
  #     }  
  #     
  #     cat('\n')
  #   }
  #   
  # }
  # 

  pio::pioTit('Coverage (DP), reads with mutant (NV) and allele frequency (VAF)')
  print(prt)
    
  if(!is.null(x$annotation))
  {
    lbl = head(unique(x$annotation$variable))
    
    # pio::pioTit('ANNOTATION')
    # print((x$annotation))
    pio::pioStr('Annotations', prefix = '\n\n',
                paste0( c(lbl, '...'), collapse = ', '))
  }
  
  if(show.log)
  {
    cat('\n')
    pio::pioTit('LOGGED OPERATIONS')
    
    print(x$operationLog)
  }
}

#' Title
#'
#' @param x 
#' @param samples 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
mobster_fit_multivariate = function(x, samples = x$samples, ...) 
{
  x$fit.MOBSTER = lapply(samples,
                         function(w)
                         {
                           pio::pioTit(paste("Fitting with MOBSTER sample", w))
                           
                           data = mobster:::VAF(x, samples = w) %>% 
                             filter(value > 0) %>% 
                             spread(variable, value)
                           
                           print(data)
                           
                           mobster_fit(
                             x = data,
                             annotation = w,
                             ...)
                         })
  
  names(x$fit.MOBSTER) = x$samples
  
  cat("\n")
  pio::pioTit("MOBSTER fits")
  
  for(s in x$samples) print(x$fit.MOBSTER[[s]]$best)
  
  # Log update
  x = logOp(x, paste0("Fit MOBSTER to", x$samples, collapse = ', '))
  
  x
}


#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' 
#' @export
#' 
#' @import sciClone
#'
#' @examples
mobster_fit_sciClone = function(x, ...) 
{
  inputs = convert_sciClone_input(x)
  
  library(sciClone)
  
  pio::pioTit("Fitting data with sciClone")
  
  sciClone.fit = sciClone(
    vafs = inputs$MUTS,
    copyNumberCalls = inputs$CN,
    sampleNames = x$samples,
    clusterMethod = 'bmm',
    ...
    )
  
  x$fit.sciClone = sciClone.fit
  
  # factor.values = sort(unique(sciClone.fit@vafs.merged$cluster))
  # data$data$sciClone.cluster = factor(sciClone.fit@vafs.merged$cluster, levels = factor.values)
  # 
  # data$sciClone.fit = sciClone.fit
  
  
  # Log update
  x = logOp(x, paste0("Fit sciClone to", x$samples, collapse = ', '))
  
  x
}


#' Title
#'
#' @param x 
#' @param peak.range 
#' @param m.peak 
#'
#' @return
#' @export
#'
#' @examples
mobster_purity = function(x, peak.range = c(0.2, 0.5), m.peak = 10)
{
  pio::pioTit(paste("Guessing purity via peak detection in", peak.range[1], '--', peak.range[2], 'for diploid SNVs'))
  
  message("This thing makes sense only if your using all diploid SNVs \n\nWill improve it to take as input one segment id!")
  
  results = NULL
  
  for(s in x$samples)
  {
    v = VAF(x, samples = s) %>% pull(value)
    
    # Detect peaks
    v = v[v > peak.range[1]]
    v = v[v < peak.range[2]]
    
    h = hist(v, breaks = seq(0, 1, 0.01), plot = FALSE)
    
    peaks = mobster:::.find_peaks(h$density, m.peak)
    x.peaks = (peaks * 0.01)
    
    # diploid tumour and normal
    guess = data.frame(samples = s, purity = 2 * x.peaks)
    
    results = rbind(results, guess)
  }
  
  results = as_tibble(results)
  
  results  
}


###########################################################################   
####################### Private auxiliary functions #######################
###########################################################################  

logOp = function(obj, op)
{
  new.entry = tribble(~time, ~operation, Sys.time(),  op)
  obj$operationLog = bind_rows(obj$operationLog, new.entry)
                               
  obj
}


# VAF Adjustment for CN an purity -- switches to CCF, gets back to VAF
vaf_adjustCN = function(v, m, M, p, mut.allele =1)
{
  CN = m+M

  0.5 * v * ((CN-2) * p + 2) / (mut.allele * p) 
}

# Conditional mutate
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

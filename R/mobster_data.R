#' Title
#'
#' @param data 
#' @param samples 
#' @param description 
#' @param DP.columns 
#' @param NV.columns 
#' @param VAF.columns 
#'
#' @return
#' @export
#'
#' @examples
mobster_dataset = function(
  data,
  samples,
  description = "My own dataset",
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
  else data$id = paste0('__mut', 1:nrow(data))
  
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
  obj = list()
  class(obj) <- "mbst_data"
  
  obj$samples = samples
  obj$N = nrow(data)
  
  obj$description = description
  
  # Annotations: anything but VAF/ CN values
  obj$annotations = tib_data %>%
    filter(!(variable %in% 
               c(all.columns, 'chr', 'from', 'to')))
  
  # Data: DP, NV and VAF
  obj$data = tib_data %>%
    filter(variable %in% all.columns)
  
  obj$data$value = as.numeric(obj$data$value)
  
  sp = strsplit(obj$data$variable, '\\.')
  obj$data$variable = sapply(sp, function(w) w[2])
  obj$data$sample = sapply(sp, function(w) w[1])
  
  # Locations: CN position of mutations
  obj$locations = tib_data %>%
    filter(variable %in%  c('chr', 'from', 'to'))
  
  # Log creation
  obj = logOp(obj, "Initialization")
  
  return(obj)
}


######### Subsample data
#' Title
#'
#' @param obj 
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
mobster_downsample = function(obj, N)
{
  if(N < obj$N)
  {
    pio::pioTit(paste0("Subsampling ", N, "entries out of ", obj$N))
    
    ids = sample(unique(obj$data$id), N)

    obj$data = obj$data %>%
      filter(id %in% ids)
    
    obj$N = N
    obj$annotations = obj$annotations %>%
      filter(id %in% ids)
    
    obj = logOp(obj, paste0("Subsampled to N = ", N, ""))
  }
  
  obj
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
                           ...) {
  
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
#' @param segments 
#' @param N.min 
#'
#' @return
#' @export
#'
#' @examples
mobster_add_CN = function(x, 
                          segments,
                          purity,
                          N.min = 500
                          ) 
{
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
    
    mapped = byLoc(x, segments_ids[s])
    if(nrow(mapped) == 0) next;
    
    mapped$seg_id = segments_ids[s]
    
    x$map_mut_seg = bind_rows(x$map_mut_seg, mapped)
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
  
  # Subset to match accepted
  x$map_mut_seg = x$map_mut_seg %>%
    filter(seg_id %in% accepted$seg_id)
  
  # Non-accepted mutations can be immediately removed
  x = delete_entries(x, 
                     setdiff(
                       keys(x), 
                       unique(x$map_mut_seg$id)))
  
  
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
          Major = Major(x, seg_to_match[seg], s),
          minor = minor(x, seg_to_match[seg], s),
          CN = Major + minor,
          purity = purity[s],
          adj_VAF = vaf_adjustCN(value, minor, Major, purity)
        )

      new.VAF = bind_rows(new.VAF, values)
    }
  }
  
  # cat("   10 Random entries\n\n")
  # print(new.VAF[sample(1:nrow(new.VAF), 10), ])
  # 
  # Save a copy of the mapping
  x$VAF_cn_adjustment = new.VAF
     
  new.VAF$value = new.VAF$adj_VAF  
  new.VAF = new.VAF %>% select(id, variable, value, sample)
  
  x$data = bind_rows(
    new.VAF,
    DP(x),
    NV(x)
  )
  
  # Log update
  x = logOp(x, "Added CN and adjusted VAF")
  
  all_zeroes(x)
}






####################### Public getters
keys = function(x) { unique(x$data$id) }

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

Annotations = function(x, ids = keys(x), 
                       variable = unique(x$annotations$variable), samples = x$samples)
{
  x$annotations %>%
    filter(variable %in% variable &  
           id %in% ids & 
           sample %in% samples)
}

####################### Private getters

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

####################### Aux functions

logOp = function(obj, op)
{
  new.entry = tribble(~time, ~operation, Sys.time(),  op)
  obj$operationLog = bind_rows(obj$operationLog, new.entry)
                               
  obj
}


# as in GEL
vaf_adjustCN = function(v, m, M, p, mut.allele =1)
{
  CN = m+M
  
  
  0.5 * v * ((CN-2) * p + 2) / (mut.allele * p) 
}

####################### Plot functions


plot_raw_data = function(x, 
                         # by_genotype = FALSE, 
                         title = x$description,
                         scales = 'fixed')
{
  by_genotype = FALSE
    
  if(!by_genotype)
  { 
    # Full VAF plot
    vaf = VAF(x) %>% filter(value > 0)
    
    max.VAF = max(vaf$value)
    if(max.VAF < 1) max.VAF = 1
    
    pl_vaf = ggplot(vaf, aes(value, fill = sample)) +
      geom_histogram(binwidth = 0.01) +
      # geom_density() +
      facet_wrap(~sample, nrow = 1, scales = scales) +
      theme_light(base_size = 8) +
      guides(fill = FALSE) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")
      ) +
      geom_vline(aes(xintercept = 0.5), colour = 'red', linetype = "longdash", size = .3) +
      geom_vline(aes(xintercept = 0.25), colour = 'red', linetype = "longdash", size = .3) +
      geom_vline(aes(xintercept = 0.33), colour = 'blue', linetype = "longdash", size = .3) +
      geom_vline(aes(xintercept = 0.66), colour = 'blue', linetype = "longdash", size = .3) +
      geom_vline(aes(xintercept = 1), colour = 'black', linetype = "longdash", size = .2) +
      # geom_density() +
      xlim(0, max.VAF) + 
      labs(
        title = "Variant Allele Frequency (VAF)",
        x = 'VAF')
    
    # Low-frequency VAF
    vaf = vaf %>% filter(value < 0.2)
    
    pl_vaf_lfreq = ggplot(vaf, aes(value, fill = sample)) +
      geom_histogram(binwidth = 0.01) +
      # geom_density() +
      facet_wrap(~sample, nrow = 1, scales = scales) +
      theme_light(base_size = 8) +
      guides(fill = FALSE) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")
      ) +
      labs(
        title = "Low-freq. VAF < 0.20",
        x = 'VAF')
    
    DP_val = DP(x) %>% filter(value > 0)
    
    stats_DP_val = quantile(DP_val %>% pull(value))
    stats_DP_val = data.frame(value = stats_DP_val, quantile = names(stats_DP_val))
    
    pl_DP = ggplot(DP_val, aes(value, fill = sample)) +
      geom_histogram(binwidth = 1) +
      geom_vline(xintercept = stats_DP_val$value, size = .3, linetype = 'dashed', show.legend = TRUE) +
      geom_vline(xintercept = median(DP_val$value), size = .5, colour = 'steelblue') +
      facet_wrap(~sample, nrow = 1,  scales = scales) +
      guides(fill = FALSE, color = guide_legend(title="Quantile")) +
      labs(
        title = "DP histogram (raw data)",
        subtitle = paste0('Quantiles (dashed):', 
                          paste0(stats_DP_val$quantile, collapse = ', '),
                          '. Median (blue)  = ', round(median(DP_val$value), 0), 
                          'x (var. ', round(var(DP_val$value), 0), ')'),
        x = 'Coverage (DP)') +
      theme_light(base_size = 8) +
      theme(legend.position="bottom",
            legend.text = element_text(size = 8),
            legend.key.size = unit(.3, "cm")
      ) 
    
    # 
    # v = DP_val %>%
    #   group_by(sample) %>%
    #   summarize(
    #     # n = n(),
    #     min = min(value),
    #     max = max(value),
    #     mean = mean(value),
    #     sd = sd(value),
    #     median = median(value)
    #   ) 
    # 
    # # v.size = v %>% select(-sample)
    # # C.v.size = apply(v.size, 2, max)
    # # v.size = data.frame(t(t(v.size)/C.v.size))
    # # v.size$sample = v$sample
    # # 
    # v = reshape2::melt(v)
    # v$value = round(v$value)
    # 
    # ggballoonplot(v, 
    #               x = "sample", y = "variable",
    #               size = "value", fill = "value", size.range = c(1, 5), ) 
    # scale_fill_gradientn(colors = my_cols) +
      # guides(size = FALSE) +
      # facet_wrap(~variable, scales = 'free')
    
    
    figure = ggpubr::ggarrange(
      ggpubr::ggarrange(pl_vaf, pl_vaf_lfreq, ncol = 1, heights = c(2, 1.5), nrow = 2),
      pl_DP,
      ncol = 1,
      nrow = 2)
    
    figure = ggpubr::annotate_figure(figure, top = title)
    
    return(figure)
  }
  
  
  values = x$VAF_cn_adjustment %>% 
    mutate(Genotype = paste0(Major, ':', minor))
  
  pio::pioStr('Available genotypes',
              paste0(unique(values$Genotype), collapse = ', '))
  
  ggplot(values, aes(value, fill = sample)) +
    geom_histogram(binwidth = 0.01) +
    facet_wrap(Genotype ~ sample, scales = 'free') +
    theme_light(base_size = 8) +
    guides(fill = FALSE) +
    theme(legend.position="bottom",
          legend.text = element_text(size = 8),
          legend.key.size = unit(.3, "cm")
    ) +
    geom_vline(aes(xintercept = 0.5), colour = 'red', linetype = "longdash", size = .3) +
    geom_vline(aes(xintercept = 0.25), colour = 'red', linetype = "longdash", size = .3) +
    geom_vline(aes(xintercept = 0.33), colour = 'blue', linetype = "longdash", size = .3) +
    geom_vline(aes(xintercept = 0.66), colour = 'blue', linetype = "longdash", size = .3) +
    geom_vline(aes(xintercept = 1), colour = 'black', linetype = "longdash", size = .2) +
    # xlim(0, max.VAF) + 
    labs(
      title = "VAF histogram (raw data)",
      x = 'VAF')
  
  stop("Not yet implemented -- plot histograms by CN genotype")
  # 
  # segments = x$segments %>% group_by(sample) %>% spread(variable, value) %>% mutate(CN = paste0(Major, ':', minor))
  # 
  # for(s in x$samples)
  # {
  #   sample.segments = segments %>% filter(sample == s)
  #   
  #   
  #   sample.segments %>% group_by(id) %>% mutate(CN = minor + Major)
  #   
  #   
  # }
  

  
}

MOBSTER_guessPurity = function(x, peak.range = c(0.2, 0.5), m.peak = 10)
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




# Fitting 
mobster_fit_multivariate = function(x, samples = x$samples, ...) {
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
  
  x
}

# Filters

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


# filter for VAF below `cutoff` parameters. This sets to 0 everything both the VAF
# and the NV entry, bute retains the coverage at the locus (DP). This function accepts
# also a purity vector in input to adjust accordingly the cutoff
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
 

 delete_entries = function(x, ids)
 {
   x$data = x$data %>% 
     filter(!id %in% ids)
   
   if(!is.null(x$map_mut_seg))
     x$map_mut_seg = x$map_mut_seg %>% 
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
     pio::pioStr("Number of entries that have VAF = 0 in all samples; N =", length(ids))
   }
     
   delete_entries(x, ids)
 }

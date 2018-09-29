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
  if(!is.data.frame(data)) stop("Data must be a dataframe")
  
  DP.columns = paste0(samples, ".DP")
  NV.columns = paste0(samples, ".NV")
  VAF.columns = paste0(samples, ".VAF")
  
  if(any(!(DP.columns %in% colnames(data)))) stop("Missing DP columns in data.")
  if(any(!(NV.columns %in% colnames(data)))) stop("Missing NV columns in data.")
  if(any(!(VAF.columns %in% colnames(data)))) stop("Missing VAF columns in data.")
  
  all.columns = c(DP.columns, NV.columns, VAF.columns)
  
  if(any(is.na(data[, all.columns, drop = FALSE])))
  {
    pio::pioTit('Removing NAs')
    
    data = data[
      complete.cases(data[, all.columns, drop = FALSE]), , drop = FALSE 
    ]
  }
  
  if('id' %in% colnames(data)) {
    warning("Data contains an $id column, will try to use it as mutation id.
            If that is not a valid unique id, random things might happen.")
  }
  else data$id = paste0('__mut', 1:nrow(data))
  
  tibdata = as_tibble(data)
  
  pio::pioTit("Creating tibble for input data")
  tibdata = tibdata %>% 
    reshape2::melt(id = 'id')

  tibdata$variable = paste(tibdata$variable)
  tibdata = as_tibble(tibdata)
  
  ##=============================##
  # Create a mbst_data object  #
  ##=============================##
  obj           = list()
  class(obj) <- "mbst_data"
  
  obj$samples = samples
  obj$N = nrow(data)
  
  obj$description = description
  
  obj$DP.columns = DP.columns
  obj$NV.columns = NV.columns
  obj$VAF.columns = VAF.columns
  
  obj$all.columns = all.columns
  
  # Get a copy of the data, where there are no VAF values etc.
  # obj$annotations = data[, !(colnames(data) %in% obj$all.columns), drop = FALSE]
  # obj$annotations$id = data$id
  
  obj$annotations = tibdata %>%
    filter(!(variable %in% all.columns))
  
  # Store data
  obj$data = tibdata %>%
    filter(variable %in% all.columns)
  
  sp = strsplit(obj$data$variable, '\\.')
  obj$data$variable = sapply(sp, function(w) w[2])
  obj$data$sample = sapply(sp, function(w) w[1])
  

  # Handy format
  # obj$data = obj$data %>% 
  #   mutate(sample = strsplit(variable, '\\.')[[1]][2]) %>% 
  #   mutate(variable = strsplit(variable, '\\.')[[1]][1]) 
  
  obj$data$value = as.numeric(obj$data$value)
  
  obj$operationLog = ""
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

logOp = function(obj, op)
{
  Sys.time()
  
  obj$operationLog = paste(
    obj$operationLog,
    paste0(Sys.time(), " | ", op, '\n')
  )
  
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
  pio::pioStr("      N", x$N, prefix = '\t- ', suffix = '\n')
  # 
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

  prt = obj$data %>%
    group_by(sample, variable) %>%
    filter(value > 0) %>%
    summarize(mean = mean(value), median = median(value), min = min(value))
  
  pio::pioTit('Coverage (DP), reads with mutant (NV) and allele frequency (VAF)')
  print(prt)
    
  if(!is.null(x$annotation))
  {
    pio::pioTit('ANNOTATION')
    print((x$annotation))
  }
  
  if(show.log)
  {
    pio::pioTit('LOGGED OPERATIONS')
    
    cat(obj$operationLog)
  }
}

mobster_add_CN = function(x, 
                          segments,
                          N.min = 500
                          ) 
{
  Major.columns = paste0(x$samples, ".Major")
  minor.columns = paste0(x$samples, ".minor")
  
  # CN
  if(
    !is.data.frame(segments) |
    !all(c('chr', 'from', 'to') %in% colnames(segments)) |
    !all(Major.columns %in% colnames(segments)) |
    !all(minor.columns %in% colnames(segments))
  ) {
    stop("Your segments do not look in the right format.")
  }
  
  # mutations 
  if(
    all(!is.null(x$annotations)) |
    any(!(c('chr', 'from', 'to') %in% x$annotations$variable))
    )
  {
    stop("To add CN you sould have included annotations in your data
         with the position of each of the annotated mutations.")
  }
  
  ######## Actual processing
  
  # get tibble
  tibsegmnets = as_tibble(segments)
  tibsegmnets$id = paste0("_CN_seg", 1:nrow(tibsegmnets))
  
  tibsegmnets = tibsegmnets %>%
    reshape2::melt(id = c('id', 'chr', 'from', 'to')) 
  
  tibsegmnets$variable = paste(tibsegmnets$variable)
  
  sp = strsplit(tibsegmnets$variable, '\\.')
  tibsegmnets$variable = sapply(sp, function(w) w[2])
  tibsegmnets$sample = sapply(sp, function(w) w[1])
  
  tibsegmnets = as_tibble(tibsegmnets)
  
  obj$data$value = as.numeric(obj$data$value)
  
  # get mutation mapping stored in the object
  muts_CN = x$annotations %>%
    filter(variable %in% c('chr', 'from', 'to')) %>%
    spread(variable, value)
  
  muts_CN$from = as.numeric(muts_CN$from)
  muts_CN$to = as.numeric(muts_CN$to)
  
  muts_CN = muts_CN %>%
    mutate(minor = NA, Major = NA)
  
  pb = txtProgressBar(min = 1, max = nrow(segments)) 
  
  pio::pioTit("Mapping mutations to CN segments")
  pio::pioStr("Number of segments", nrow(segments), suffix = '\n')
  pio::pioStr("         Minimum N", N.min, suffix = '\n')
  
  tib_new.data = NULL
  
  for(s in unique(tibsegmnets$id))
  {
    setTxtProgressBar(pb, s)
    
    thisSeg = tibsegmnets %>%
      filter(id == s)
    
    which_muts = muts_CN %>%
      filter(chr == thisSeg$chr[1], 
             from >= thisSeg$from[1], 
             to <= thisSeg$to[1])
    
    pio::pioStr(sprintf('\n%15s', s), nrow(which_muts), '\n\n')
    
    # cat("Segment", s, "N = ", nrow(which_muts), '\n')
    
    if(nrow(which_muts) == 0 |
       nrow(which_muts) < N.min
       ) next;
    
    new.data = inner_join(which_muts, x$data, by = 'id') %>%
      filter(variable == 'VAF')
    
    toPrint = NULL
    
    for(samp in x$samples) {
      
      new.data$minor = thisSeg %>% filter(sample == samp, variable == 'minor') %>% pull(value)
      new.data$Major = thisSeg %>% filter(sample == samp, variable == 'Major') %>% pull(value)
      
      new.values = new.data %>%
          filter(sample == samp) %>%
          mutate(adj_VAF = 
                   vaf_adjustCN(value,
                                m = 1, # minor,
                                M = Major,
                                p = 1)
          )

      tib_new.data = bind_rows(tib_new.data, new.values)
      toPrint = bind_rows(toPrint, new.values[1, ])
    }
    
    print(toPrint)
  }
  
  tib_new.data$value = tib_new.data$adj_VAF  
  
  # merge to the old data where we change the VAF value
  tib_old = x$data
  tib_old = tib_old %>%
    filter(variable != 'VAF')
  
  tib_old = bind_rows(tib_old,
    tib_new.data %>%
      select(id, variable, value)
  )
      
  
  tib_old
}

# m= 1!
vaf_adjustCN = function(v, m, M, p)
{
  v * ((M-2) * p + 2) / (m*p)
}

map_SNV_2_CNA_segment = function(segments,
                                 segment.id,
                                 muts
                                 )
{
  stopifnot(all(unlist(labels$Chromosome) %in% colnames(segments)))
  stopifnot(all(unlist(labels$Mutations) %in% colnames(muts)))
  
  map = NULL
  map = muts[muts[, labels$Mutations['chromosome']] == segments[segment.id, labels$Chromosome['chromosome']], ]
  map = map[map[, labels$Mutations['position']] >= segments[segment.id, labels$Chromosome['from']], ]
  map = map[map[, labels$Mutations['position']] <= segments[segment.id, labels$Chromosome['to']], ]
  
  map
}
  

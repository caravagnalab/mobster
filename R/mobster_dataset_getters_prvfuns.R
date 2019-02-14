####################### Private getters

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

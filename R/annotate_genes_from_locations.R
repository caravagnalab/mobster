# This function uses
# GenomicRanges
# GenomicFeatures
# TxDb.Hsapiens.UCSC.hg19.knownGene
# org.Hs.egSYMBOL
annotate_genes_from_locations = function(x)
{
  crash_ifnotinstalled(c('Homo.sapiens', 'GenomicRanges', 'GenomicFeatures', 'TxDb.Hsapiens.UCSC.hg19.knownGene', 'org.Hs.egSYMBOL'))
    
  # This loads most of the required stuff
  require(Homo.sapiens)
  
  if(!all(c('chr', 'from', 'to', 'ref', 'alt') %in% colnames(x$data)))
    stop("Missing genomic coordinates (chr, from, ref, alt) from the input data, cannot annotate genes.")
  
    if(any(!grepl('chr', x$data$chr)))
    {
      message("Adding `chr` to the chromosome names (1 -> chr1)")
      x$data = x$data %>%
        dplyr::mutate(
          chr = ifelse(grepl('chr', chr), chrom, paste0('chr', chr))
        )
    }
    
  mutations = mobster::Clusters(x) %>%
    dplyr::select(chr, from, ref, alt) %>%
    dplyr::rename(chrom = chr, start = from) %>% 
    dplyr::mutate(
      start = as.numeric(start),
      end = start + nchar(alt)
    ) %>%
    dplyr::select(chrom, start, end) 
  
  GR_df_mutations = GenomicRanges::makeGRangesFromDataFrame(
    mutations %>% 
    data.frame(stringsAsFactors = FALSE)
    )
  
  GR_df_overlaps = GenomicRanges::subsetByOverlaps(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene), GR_df_mutations) %>%
    data.frame(stringsAsFactors = FALSE) %>% 
    tibble::as_tibble() %>%
    dplyr::mutate(
      seqnames = paste(seqnames),
      strand = paste(strand)
    ) %>%
    dplyr::left_join(as.data.frame(org.Hs.egSYMBOL), by = 'gene_id') %>%
    dplyr::select(-strand, -gene_id, -width)
  
  colnames(GR_df_overlaps) = c('chr', 'from', 'to', 'gene')
  colnames(mutations) = c('chr', 'from', 'to')
  
  final_mapping = lapply(1:nrow(GR_df_overlaps),
                         function(x){
                           mutations %>%
                             dplyr::filter(
                               chr == GR_df_overlaps$chr[x],
                               from >= GR_df_overlaps$from[x],
                               to <= GR_df_overlaps$to[x]
                             ) %>%
                             dplyr::mutate(gene = GR_df_overlaps$gene[x])
                         }) %>%
    Reduce(f = bind_rows)
  
  mutations_enriched = mutations %>%
    dplyr::left_join(final_mapping, by = c('chr', 'from', 'to')) %>%
    dplyr::distinct(chr, from, to, .keep_all = T)
  
  # table(mutations_enriched$gene) %>% sort(decreasing = TRUE) %>% pioDisp()
  mobster::Clusters(x)  %>%
    dplyr::mutate(
      from = as.numeric(from),
      to = from + nchar(ref)
      ) %>%
    dplyr::left_join(mutations_enriched, by = c('chr', 'from', 'to')) %>%
    dplyr::select(chr, from, to, ref, alt, gene, everything())
}
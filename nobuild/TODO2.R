# # biocLite("Homo.sapiens")
# library(Homo.sapiens)
# genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
# 
# M = apply(mutations, 1, paste0, collapse = ':', sep ='') 
# M = Reduce(append, M)
# 
# mycoords.list =list('2:2542384:2542385')
# mycoords.list = list("1:4864876:5864876", "1:14283067:15283067", "1:21786817:22786817", 
#      "1:33465769:34465769", "1:45300539:46300539", "1:54333815:55333815", 
#      "1:65114236:66114236", "1:75194833:76194833", "1:86037468:87037468", 
#      "1:96462256:97462256", "1:105259436:106259436", "1:116234756:117234756", 
#      "1:120842170:121842170", "1:145808064:146808064", "1:155459582:156459582", 
#      "1:166112356:167112356", "1:174453227:175453227", "1:185347260:186347260", 
#      "1:194299241:195299241", "1:205731116:206731116")
# 
# mycoords.list = mycoords.list[1:10]
# mycoords.gr = lapply(mycoords.list, function (x) {res=strsplit(x, ':')}) %>%
#   unlist %>%
#   as.numeric %>%
#   matrix(ncol=3, byrow=T) %>%
#   as.data.frame %>%
#   dplyr::select(chrom=V1, start=V2, end=V3) %>%
#   mutate(chrom=paste0('chr', chrom)) %>%
#   makeGRangesFromDataFrame
# 
# overlaps = subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), mycoords.gr)
# 
# genes = data.frame(gene_id = overlaps$gene_id, stringsAsFactors = FALSE) %>%
#   left_join(as.data.frame(org.Hs.egSYMBOL), by = 'gene_id') %>%
#   as_tibble()
# 
# library(tidyverse)
# mutations = mobster::LUFF76_lung_sample$best$data %>%
#   dplyr::select(chr, from, ref, alt) %>%
#   dplyr::rename(chrom = chr, start = from) %>% 
#   dplyr::mutate(
#     start = as.numeric(start),
#     end = start + nchar(alt)
#   ) %>%
#   dplyr::select(chrom, start, end) 
# 
# GR_df_mutations = mutations %>% 
#   data.frame(stringsAsFactors = FALSE) %>%
#   makeGRangesFromDataFrame 
# 
# GR_df_overlaps = subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), GR_df_mutations) %>%
#   data.frame(stringsAsFactors = FALSE) %>% 
#   as_tibble() %>%
#   mutate(
#     seqnames = paste(seqnames),
#     strand = paste(strand)
#     ) %>%
#   left_join(as.data.frame(org.Hs.egSYMBOL), by = 'gene_id') %>%
#   select(-strand, -gene_id, -width)
# 
# colnames(GR_df_overlaps) = c('chr', 'from', 'to', 'gene')
# colnames(mutations) =c('chr', 'from', 'to')
# 
# final_mapping = lapply(1:nrow(GR_df_overlaps),
#        function(x){
#          mutations %>%
#            filter(
#              chr == GR_df_overlaps$chr[x],
#              from >= GR_df_overlaps$from[x],
#              to <= GR_df_overlaps$to[x]
#            ) %>%
#            mutate(gene = GR_df_overlaps$gene[x])
#        }) %>%
#   Reduce(f = bind_rows)
#          
# mutations %>%
#   left_join(final_mapping, by = c('chr', 'from', 'to'))
#   

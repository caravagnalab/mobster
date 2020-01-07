# 
# # annotate(
# library(evoverse.datasets)
# data('TRACERx_NEJM_2017', package = 'evoverse.datasets')
# 
# data('LUFF76_lung_sample', package = 'mobster')
# 
# 
# library(tidyverse)
# mutations = LUFF76_lung_sample$best$data %>%
#   # filter(row_number() < 10) %>%
#   mutate(sample = 'XXX') %>%
#   select(sample, chr, from, ref, alt) %>%
#   data.frame(stringsAsFactors = FALSE)
# 
# x = mutations %>% as_tibble()
# 
# if(any(grepl('chr', x$chr)))
#   x$chr = gsub(x = x$chr, 'chr', '')
# 
# x = x %>%
#   mutate(
#     sample = paste(sample),
#     from = as.numeric(from)
#   ) %>%
#   filter(ref != alt) %>%
#   rename(
#     sampleID = sample,
#     pos = from,
#     mut = alt
#   )
# 
# gene_list = c("TP53", "CDKN2A", "KRAS", "EGFR")
# gene_list = mobster::cancer_genes_dnds$Tarabichi_drivers
# 
# # mutations = mutations[,1:5] # Restricting input matrix to first 5 columns
# # mutations[,c(1,2,3,4,5)] = lapply(mutations[,c(1,2,3,4,5)], as.character) # Factors to character
# # mutations[[3]] = as.numeric(mutations[[3]]) # Chromosome position as numeric
# # mutations = mutations[mutations[,4]!=mutations[,5],] # Removing mutations with identical reference and mutant base
# # colnames(mutations) = c("sampleID","chr","pos","ref","mut")
# 
# 
# # [Input] Reference database
# data("refcds_hg19", package="dndscv")
# if (any(gene_list=="CDKN2A")) { # Replace CDKN2A in the input gene list with two isoforms
#   gene_list = unique(c(setdiff(gene_list,"CDKN2A"),"CDKN2A.p14arf","CDKN2A.p16INK4a"))
# }
# 
# # [Input] Gene list (The user can input a gene list as a character vector)
# if (is.null(gene_list)) {
#   gene_list = sapply(RefCDS, function(x) x$gene_name) # All genes [default]
# } else { # Using only genes in the input gene list
#   allg = sapply(RefCDS,function(x) x$gene_name)
#   nonex = gene_list[!(gene_list %in% allg)]
#   if (length(nonex)>0) { stop(sprintf("The following input gene names are not in the RefCDS database: %s", paste(nonex,collapse=", "))) }
#   RefCDS = RefCDS[allg %in% gene_list] # Only input genes
#   gr_genes = gr_genes[gr_genes$names %in% gene_list] # Only input genes
# }
# 
# # Expanding the reference sequences [for faster access]
# for (j in 1:length(RefCDS)) {
#   RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
#   RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
#   RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
#   if (!is.null(RefCDS[[j]]$seq_splice)) {
#     RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
#     RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
#     RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
#   }
# }
# 
# 
# ## 2. Mutation annotation
# message("[2] Annotating the mutations...")
# 
# ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
# gr_genes_ind = ind[gr_genes$names]
# 
# # Warning about possible unannotated dinucleotide substitutions
# if (any(diff(mutations$pos)==1)) {
#   warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
# }
# 
# # Warning about multiple instances of the same mutation in different sampleIDs
# if (nrow(unique(mutations[,2:5])) < nrow(mutations)) {
#   warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
# }
# 
# # Start and end position of each mutation
# x$end = x$start = x$pos
# l = nchar(x$ref)-1 # Deletions of multiple bases
# x$end = x$end + l
# ind = substr(x$ref,1,1)==substr(x$mut,1,1) & nchar(x$ref)>nchar(x$mut) # Position correction for deletions annotated in the previous base (e.g. CA>C)
# x$start = x$start + ind
# 
# 
# # Mapping mutations to genes
# gr_muts = GenomicRanges::GRanges(x$chr,
#                                  IRanges::IRanges(x$start, x$end))
# 
# ol = as.data.frame(GenomicRanges::findOverlaps(gr_muts, gr_genes, type =
#                                                  "any", select = "all"))
# x = x[ol[, 1], ] # Duplicating subs if they hit more than one gene
# x$geneind = gr_genes_ind[ol[, 2]]
# x$gene = sapply(RefCDS, function(y)
#   y$gene_name)[x$geneind]
# x = unique(x)
# 
# 
# 
# 

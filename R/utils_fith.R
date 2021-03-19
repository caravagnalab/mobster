tensorize <- function(x){

  torch <- reticulate::import("torch")
  x <- lapply(x, function(y) torch$tensor(y))
  return(x)
}

get_purity <- function(x){

  return(x$cnaqc$purity)

}

format_data_mobsterh_QC <-  function(x, kar = c("1:0", "1:1", "2:1", "2:0", "2:2"), vaf_t = 0.05){

  valid_karyo <- x$QC$QC_table %>% dplyr::filter(QC == "PASS", type == "Peaks") %>% pull(karyotype)

  valid_karyo <-  intersect(kar, valid_karyo)

  res <- x$cnaqc$snvs %>% filter(karyotype %in% valid_karyo,Variant_Type == "SNP") %>%  filter(VAF >= vaf_t) %>% mutate(id = paste(chr, from, to, sep = ":")) %>% select(VAF, karyotype, id)

  return(split_and_tolist(res))

}



format_data_mobsterh_DF <-  function(x, kar = c("1:0", "1:1", "2:1", "2:0", "2:2"), vaf_t = 0.05){

  res <- x$cnaqc$snvs %>% filter(karyotype %in% kar,VAF >= 0.05, VAF < 1, VAF > 0) %>% mutate(id = paste(chr, from, to, sep = ":")) %>% select(VAF, karyotype, id)

  return(split_and_tolist(res))

}


split_and_tolist <- function(x){
  res <- split(x, x$karyotype)

  nm <-  lapply(res, function(y) y$id)
  res <-  lapply(res, function(y) y$VAF)

  for(i in 1:length(res))
    names(res[[i]]) <- nm[[i]]

  filt <- sapply(res, function(y) length(y) > 100)

  nm <- names(res)[which(filt)]
  res <- as.list(res[which(filt)])
  names(res) <- nm

  return(res)
}

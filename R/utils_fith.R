tensorize <- function(x) {
  torch <- reticulate::import("torch")
  x <- lapply(x, function(y)
    torch$tensor(y))
  return(x)
}

get_purity <- function(x) {
  return(x$cnaqc$purity)
}

format_data_mobsterh_QC <-
  function(x,
           kar = c("1:0", "1:1", "2:1", "2:0", "2:2"),
           vaf_t = 0.05,
           n_t = 100,
           NV_filter = 5,
           enforce_QC_PASS = TRUE
           ) {
    if (enforce_QC_PASS)
      valid_karyo <-
        x$QC$QC_table %>% dplyr::filter(QC == "PASS", type == "Peaks") %>% pull(karyotype)
    else
      valid_karyo <-  kar

    valid_karyo <-  intersect(kar, valid_karyo)

    if (length(valid_karyo) < 1) {
      reason = case_when(
        is.null(x$QC$QC_table) ~ "There are no QC tables for input 'x', rerun the data QC pipeline",
        (nrow(Peaks_entries) == 0) ~ "There are no 'Peaks' in the QC tables for input 'x', rerun the data QC pipeline",
        all(Peaks_entries$QC != "PASS") ~ "All peaks in the input data are failed, check your segmentation and CN calls.",
        TRUE ~ paste0(
          "Unknown error - the following karyotypes are PASS: ",
          paste(QC_peaks, collapse = ', '),
          '.'
        )
      )

      cat("\n")
      cat(
        cli::boxx(
          paste0("There is nothing to perform deconvolution here! ", reason),
          padding = 1,
          col = 'white',
          float = 'center',
          background_col = "brown"
        )
      )
      cat("\n")
      return(NULL)
    }

    res <-
      x$cnaqc$snvs %>% filter(karyotype %in% valid_karyo, type == "SNV") %>%
      mutate(VAF = NV / DP) %>%  filter(NV > NV_filter) %>%
      filter(VAF >= vaf_t, VAF < 1) %>% mutate(id = paste(chr, from, to, sep = ":")) %>%
      select(NV, DP, karyotype, id)

    valid_k_n <-
      res %>%  dplyr::group_by(karyotype) %>% dplyr::summarize(n = dplyr::n()) %>%  dplyr::filter(n > n_t) %>% dplyr::pull(karyotype)

    nremoved <- length(res$karyotype %>%  unique()) - length(valid_k_n %>% unique())
    if(nremoved > 0 )
      cli::cli_alert_warning("Removing {length(nremoved)}  karyotypes, containing less than {n_t} mutations")

    return(split_and_tolist(res %>% filter(karyotype %in% valid_k_n)))

  }



format_data_mobsterh_DF <-
  function(x,
           kar = c("1:0", "1:1", "2:1", "2:0", "2:2"),
           vaf_t = 0.05,
           NV_filter = 5,
           n_t = 100) {

    if("cluster" %in% colnames(x)){
      cli::cli_alert_warning("A coloumn names cluster already exists, overwriting it!")
      x$cluster <-  NULL
    }

    res <- x %>%
      mutate(VAF = NV / DP) %>%
      filter(karyotype %in% kar, VAF > vaf_t, VAF < 1, VAF > 0, NV > NV_filter) %>%
      mutate(id = paste(chr, from, to, sep = ":")) %>%
    select(NV, DP, karyotype, id)

    valid_k_n <-
      res %>%  dplyr::group_by(karyotype) %>% dplyr::summarize(n = dplyr::n()) %>%  dplyr::filter(n > n_t) %>% dplyr::pull(karyotype)
    nremoved <- length(res$karyotype %>%  unique()) - length(valid_k_n %>% unique())
    if(nremoved > 0 )
      cli::cli_alert_warning("Removing {length(nremoved)} karyotypes, containing less than {n_t} mutations")

    return(split_and_tolist(res %>% filter(karyotype %in% valid_k_n)))

  }


split_and_tolist <- function(x) {
  res <- split(x, x$karyotype)

  nm <-  lapply(res, function(y)
    y$id)
  res <-  lapply(res, function(y)
    cbind(y$NV, y$DP))

  for (i in 1:length(res))
    rownames(res[[i]]) <- nm[[i]]

  nm <- names(res)
  res <- as.list(res)
  names(res) <- nm

  return(res)
}


has_tail <-  function(x) {
  if (is_mobsterhL(x))
    return(!is.null(x$model_parameters[[1]]$tail_scale))
  else
    return(x$fit.tail)

}

has_subclones <- function(x) {
  return(x$run_parameters$K > 0)
}

is_moyal <- function(x) {
  return(x$run_parameters$subclonal_prior == "Moyal")
}

is_truncated <-  function(x) {
  return(has_tail(x) & x$run_parameters$truncated_pareto)

}

get_simple_karyotypes <- function() {
  
  return(c("1:0", "1:1", "2:0", "2:1", "2:2") )
  
}

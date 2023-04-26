tensorize <- function(x) {
  torch <- reticulate::import("torch")
  x <- lapply(x, function(y)
    torch$tensor(y)$float())
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
           enforce_QC_PASS = TRUE,
           filter_indels = TRUE
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
      x$cnaqc$mutations %>% filter(karyotype %in% valid_karyo) %>%
      mutate(VAF = NV / DP) %>%  filter(NV > NV_filter) %>%
      filter(VAF >= vaf_t, VAF < 1) %>% mutate(id = paste(chr, from, to, sep = ":")) %>%
      select(NV, DP, karyotype, id)
    
    if(filter_indels) res <- res %>%  filter(type == "SNV")

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
           n_t = 100, 
           filter_indels = TRUE) {

    if("cluster" %in% colnames(x)){
      cli::cli_alert_warning("A coloumn names cluster already exists, overwriting it!")
      x$cluster <-  NULL
    }

    if(filter_indels) x <- x %>% filter(to - from == 1)
    
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


filter_vaf <- function(x, min_VAF = 0, max_VAF = 1){
  stl <- split_and_tolist(x$data %>% dplyr::filter(id %in% x$used_mutations) )
  stl_nms <- lapply(stl,function(x) rownames(x))
  x$data <- x$data %>% filter(VAF > min_VAF, VAF < max_VAF)
  x$model_parameters <- mapply(x$model_parameters,stl_nms,FUN =  function(y,z) {
    colnames(y$cluster_probs) <- z
    y$cluster_probs <- y$cluster_probs[,colnames(y$cluster_probs) %in% x$data$id]
    names(y$cluster_assignments) <- z
    y$cluster_assignments <- y$cluster_assignments[names(y$cluster_assignments) %in% x$data$id]
    return(y)
  }, SIMPLIFY = FALSE)
  return(x)
}

mutation_rate_by_VAF_cut <- function(x, VAF_interval = seq(0.01, 0.1, by = 0.01), ...) {
  
  mu_table <- lapply(VAF_interval, function(i) mu_posterior(filter_vaf(x, i), ...)) %>% do.call(rbind,.)
  mu_table$min_VAF <- VAF_interval
  elbow_slope <- c(NA, diff(mu_table$mean)/diff(mu_table$min_VAF)) %>% abs()
  elbow_slope[1] <- 0 
  mu_table$slope <- elbow_slope
  elbow_point <- mu_table$min_VAF[(elbow_slope %>% which.max()) - 1]
  plot <- ggplot(mu_table, aes(y = mean, x = min_VAF )) + geom_line() + 
    geom_point(aes(color = min_VAF == min_VAF[(elbow_slope %>% which.max()) - 1] )) +
    theme_bw() + scale_color_manual("Elbow point?", values = c("grey60", "indianred"))

  return(list(mu_table, plt = plot, elbow_point = elbow_point))

} 




#' Run a dN/dS analysis on Mobster clusters
#'
#' @description This function takes a MOBSTER fit and runs dndscv (https://github.com/im3sanger/dndscv) to calculate dN/dS 
#' values per cluster. It computes global dN/dS and per gene dN/dS values and makes a plot. dN/dS values are computed with
#' the best fitting MOBSTER model.
#' 
#' @param fit MOBSTER fit object
#' @param gene_list An optional vector of gene names to infer dN/dS values, 
#' default is to use the whole exome
#' @param refdb The genome referene to use, default is to use hg19. Other references are available from 
#' https://github.com/im3sanger/dndscv_data
#' 
#' @return the fit object is returned with additional dnds list to the best fit model.
#'
#'
#'
#' @export
dnds <- function(fit, gene_list = NULL, refdb = "hg19"){
  cl <- Clusters(fit$fit$best) %>%
    dplyr::mutate(dummysample = "sample") %>%
    dplyr::select(dummysample, CHROM, POS,REF,ALT, everything())
  
  if (refdb == "hg19" & stringr::str_detect(cl$CHROM[1], "chr")){
    message("Removing chr from chromosome names for hg19 reference compatability")
    cl$CHROM <- stringr::str_sub(cl$CHROM, 4)
  }
  
  clusters <- unique(cl$cluster)
  globaldndstable <- tibble::tibble()
  dndscvtable <- tibble::tibble()
  
  for (i in clusters){
    message(paste0("Calculating dN/ds values for ", i))
    dndsout <- cl %>%
      dplyr::filter(cluster == i) %>%
      dndscv::dndscv(., gene_list = gene_list)
    globaldndstable <- dplyr::bind_rows(globaldndstable, dndsout$globaldnds %>%
                                   mutate(cluster = i))
    dndscvtable <- dplyr::bind_rows(dndscvtable, dndsout$sel_cv %>%
                                   mutate(cluster = i))
  }
  
  globaldndsplot <- globaldndstable %>%
    filter(name == "wall") %>%
    ggplot2::ggplot(ggplot2::aes(x = cluster, y = mle, ymin = cilow, ymax = cihigh)) +
    #ggplot2::geom_point() +
    ggplot2::geom_pointrange() +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylab("dN/dS") +
    ggplot2::geom_hline(yintercept = 1.0, lty = 2) +
    ggplot2::ggtitle("Global dN/dS values by cluster")
  
  fit$fit$best$dnds <- list(globaldndstable = globaldndstable, dndscvtable = dndscvtable, dndsplot = globaldndsplot)
  return(fit)
}
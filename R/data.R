#' Example MOBSTER fit.
#'
#' @description Data from an example MOBSTER fit; this object is the
#' result of running `mobster_fit` function on the input dat.a
#' The ouput is a list which contains the best fit, the top runs and 
#' a table summarising the fit scores. The input data has been simulated
#' with a stochastic branching process model.
#' 
#' @docType data
#'
#' @usage data(fit_example)
#'
#' @format Output from `mobster_fit`.
#'
#' @keywords datasets
#'
#' @examples
#' data(fit_example)
#' print(fit_example$best)
#' plot(fit_example$best)
"fit_example"


#' List of cancer genes to compute dnds values.
#'
#' @description List of gene ids (HuGO format) to compute dnds values.
#' There are 2 lists available of putative driver genes, and 2 lists of
#' essential genes; these lists have been subset to include only
#' genes in the RefCDS database.
#' 
#' Drivers have been compiled in 
#' `Martincorena, et al. Cell 171.5 (2017): 1029-1041`, and in 
#' `Tarabichi, et al. Nature Genetics 50.12 (2018): 1630`.
#' Essential genes have been compiled using two different cell lines in
#' `Wang et al. Science 350.6264 (2015): 1096-1101.` and
#' `Bloomen et al. Science 350.6264 (2015): 1092-1096`. All lists are 
#' available and named accordingly; use `names(cancer_genes_dnds)` to
#' see the available names.
#'
#' @docType data
#'
#' @usage data(cancer_genes_dnds)
#'
#' @format List of cancer genes to compute dnds values.
#'
#' @keywords datasets
#'
#' @examples
#' data(cancer_genes_dnds)
#' print(lapply(cancer_genes_dnds, head))
"cancer_genes_dnds"


#' MOBSTER fit for the LUFF76 lung sample.
#'
#' @description MOBSTER fit for the LUFF76 lung sample; this object is the
#' result of running `mobster_fit` function on the input data described in the 
#' main MOBSTER paper. The data consists of diploid mutations for the 
#' sample available at the Comprehensive Omics Archive of Lung Adenocarcinoma
#' (http://genome.kaist.ac.kr/).
#'
#' @docType data
#'
#' @usage data(LUFF76_lung_sample)
#'
#' @format Output from `mobster_fit`.
#'
#' @keywords datasets
#'
#' @examples
#' data(LUFF76_lung_sample, package = 'mobster')
#' print(LUFF76_lung_sample$best)
#' plot(LUFF76_lung_sample$best)
"LUFF76_lung_sample"


#' MOBSTER fit for the LU4 lung sample.
#'
#' @description MOBSTER fit for the LU4 lung sample; this object is the
#' result of running `mobster_fit` function on the input data described in the 
#' main MOBSTER paper. The data consists of diploid mutations for the 
#' sample available at the Comprehensive Omics Archive of Lung Adenocarcinoma
#' (http://genome.kaist.ac.kr/).
#'
#' @docType data
#'
#' @usage data(LU4_lung_sample)
#'
#' @format Output from `mobster_fit`.
#'
#' @keywords datasets
#'
#' @examples
#' data(LU4_lung_sample, package = 'mobster')
#' print(LU4_lung_sample$best)
#' plot(LU4_lung_sample$best)
"LU4_lung_sample"


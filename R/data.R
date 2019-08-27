#' Example MOBSTER fit.
#'
#' @description Data from an example MOBSTER fit. This is a list returned by the fitting
#' algorithm, which contains the best fit, the top runs and a table summarising the fit scores.
#'
#' @docType data
#'
#' @usage data(fit_example)
#'
#' @format Data from an example MOBSTER fit.
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
#' @description List of gene ids (HuGO format) developed from Martincorena et al. to compute dnds values.
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
#' print(cancer_genes_dnds$best)
"cancer_genes_dnds"

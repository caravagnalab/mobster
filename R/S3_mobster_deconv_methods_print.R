#' Summary for an object of class \code{'mobster_deconv'} is a print.
#'
#' @param object An obj of class \code{'mobster_deconv'}.
#' @param ...
#'
#' @return See \code{\link{print}}.
#' @exportS3Method summary mobster_deconv
#' @export summary.mobster_deconv
#'
#' @examples
#' data(fit_example)
#' summary(fit_example)
summary.mobster_deconv = function(object, ...) {
  print.mobster_deconv(object, ...)
}

#' Print a MOBSTER[h] model-selection object.
#'
#' @param x An obj of class \code{'mobster_deconv'}.
#' @param ...
#'
#' @return nothing.
#' @exportS3Method print mobster_deconv
#' @export print.mobster_deconv
#' @importFrom clisymbols symbol
#'
#' @examples
#' data("fit_example_mobsterh", package = 'mobster')
print.mobster_deconv = function(x, ...)
{
  stopifnot(inherits(x, "mobster_deconv"))

  ####  INFORMATION ABOUT THE BEST MODEL

  cli::cli_h1("MOBSTER[h] subclonal deconvolution")
  cat('\n')

  print(x$best)

  cli::cli_h3("Model selection with other {.value {x$fits.table %>% nrow()}} fits.")
}

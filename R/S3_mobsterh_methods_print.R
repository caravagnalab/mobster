#' Summary for an object of class \code{'dbpmmh'} is a print.
#'
#' @param object An obj of class \code{'dbpmmh'}.
#' @param ...
#'
#' @return See \code{\link{print}}.
#' @exportS3Method summary dbpmmh
#' @export summary.dbpmmh
#'
#' @examples
#' data(fit_example)
#' summary(fit_example$best)
summary.dbpmmh = function(object, ...) {
  print.dbpmmh(object, ...)
}

#' Print a MOBSTERh object.
#'
#' @param x An obj of class \code{'dbpmmh'}.
#' @param ...
#'
#' @return nothing.
#' @exportS3Method print dbpmmh
#' @export print.dbpmmh
#' @importFrom clisymbols symbol
#'
#' @examples
print.dbpmmh = function(x, ...)
{
  stopifnot(inherits(x, "dbpmmh"))

  ####  INFORMATION ABOUT THE BEST MODEL

  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black("[ MOBSTERh ] {.value {x$description}}")),
      '{.value {length(x$model_parameters)}} with {.field S = {x$run_parameters$K}} subclonal Beta(s) {.value {ifelse(x$run_parameters$tail == 1, crayon::green("and a tail"), crayon::red("without tail"))}}'
    )
  )



  #clonality_interpreter(x) %>% print()
}

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
#' data("fit_example_mobsterh", package = 'mobster')
print.dbpmmh = function(x, ...)
{
  stopifnot(inherits(x, "dbpmmh"))

  ####  INFORMATION ABOUT THE BEST MODEL

  cli::cli_rule(
    paste(
      crayon::bgMagenta(crayon::white("[ MOBSTERh ] {.value {x$description}}")),
      '{.field {length(x$model_parameters)}} karyotypes, {.field {x$run_parameters$K}} subclonal Beta(s) {.value {ifelse(x$run_parameters$tail == 1, crayon::green("with tail"), crayon::red("without tail"))}}'
    )
  )

  # x$model_parameters %>% names
  # x$run_parameters
  # get_beta(x)
  # mobster:::get_pareto(x)

  cln = clonality_interpreter(x)
  cln = ifelse(
    "Subclone" %in% cln$what,
    " Polyclonal " %>% crayon::magenta() %>% crayon::bgWhite(),
    ' Monoclonal ' %>% crayon::black() %>% crayon::bgGreen()
  )

  cli::cli_h3("Best model is {.value {crayon::bold(cln)}}")
  cat('\n')

  sapply(x$model_parameters %>% names,
         function(y) {
           w = clonality_interpreter(x) %>%
             filter(karyotype == y) %>%
             pull(cluster) %>%
             sort()

           w_add = function(cl)
           {
             x$data %>%
               filter(is_driver, karyotype == y, cluster == cl) %>%
               pull(driver_label) %>%
               paste(sep = ', ')
           }

           w_map = sapply(w, w_add)

           # Add label fun
           w_addl = function(cl){
             if(w_map[[cl]] %>% length() == 0)
               return(paste("", cl, ""))
             paste0(" ", cl, ' (', w_map[[cl]],') ')
           }

           # Clonal etc.
           w_cl = w[grep("C", w)]
           w_cl = sapply(w_cl, w_addl)  %>% crayon::bgBlue()


           w_sl = w[grep("S", w)]
           w_sl = sapply(w_sl, w_addl)  %>% crayon::bgMagenta()

           w_t = w[grep("T", w)]
           w_t = sapply(w_t, w_addl)  %>% crayon::bgYellow() %>% crayon::black()

           if (length(w_sl) == 0)
             w_sl = 'x' %>% crayon::red()

           cli::cli_alert('{.field {y}}: {.value {w_cl}} + {.value {w_sl}} + {.value {w_t}}')

         })

  # cli::cli_h3("Drivers")

  if(x$data %>% filter(is_driver) %>% nrow == 0)
  {
    cli::cli_alert_danger(
    "Missing drivers annotations. \\
    Consider using CNAqc [{crayon::italic('caravagnalab.github.io')}] to annotate them.")
    return()
  }

  n_d = x$data %>% filter(is_driver, !is.na(cluster)) %>% nrow
  if(n_d > 0)
  {
    cat('\n')
    cli::cli_alert_success("Succesfully clusterd {crayon::green(n_d)} driver(s).")
    cat('\n')
  }

  n_dna = x$data %>% filter(is_driver, is.na(cluster)) %>% nrow
  if(n_dna > 0)
  {
    cli::cli_alert_warning("Cannot cluster {crayon::red(n_dna)} driver(s).")
    cat('\n')

    # writeLines(paste0("  ", capture.output(print(df, row.names = F))))

    writeLines(paste0("      ",
                      capture.output(
                        x$data %>% filter(is_driver, is.na(cluster)) %>%
                          select(chr, from, to, ref, alt, driver_label, karyotype, cluster) %>%
                          print(row.names = F)
                      )))

  }

  #clonality_interpreter(x) %>% print()
}

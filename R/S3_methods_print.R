#' Summary for an object of class \code{'dbpmm'} is a print.
#'
#' @param x An obj of class \code{'dbpmm'}.
#' @param ...
#'
#' @return See \code{\link{print}}.
#' @export
#'
#' @examples
#' data(fit_example)
#' summary(fit_example$best)
summary.dbpmm = function(x, ...) {
  print.dbpmm(x, ...)
}

#' Summaries for an object of class \code{'dbpmm'} is like a print.
#'
#' @param x An obj of class \code{'dbpmm'}.
#' @param ...
#'
#' @return nothing.
#' @export
#' @import crayon
#' @import cli
#'
#' @examples
#' data(fit_example)
#' print(fit_example$best)
print.dbpmm = function(x, ...)
{
  stopifnot(inherits(x, "dbpmm"))
  
  # cli::cli_rule(
  #   paste(
  #     crayon::bgYellow(crayon::black("[ MOBSTER ]")),
  #     'n = {.value {x$N}}, fit by {.field {x$fit.type}} in {.value {length(x$all.NLL)}} steps',
  #     ifelse(
  #       x$status,
  #       crayon::green('(converged).'),
  #       crayon::red('(interrupted).')
  #     )
  #   )
  # )
  cli::cli_rule(
    paste(
      crayon::bgYellow(crayon::black("[ MOBSTER ] {.value {x$description}}")),
      'n = {.value {x$N}}, fit by {.field {x$fit.type}} in {.value {length(x$all.NLL)}} steps',
      ifelse(
        x$status,
        crayon::green('(converged).'),
        crayon::red('(interrupted).')
      )
    )
  )
  
  ####################### Pi
  clus.size = table(x$data$cluster)
  clus.size = clus.size[order(clus.size)]
  
  clus.size.pi = x$Clusters %>%
    dplyr::select(-init.value) %>%
    dplyr::filter(type == 'Mixing proportion') %>%
    dplyr::mutate(fit.value = formatC(fit.value, digits = 2))
  
  pi.tail = (clus.size.pi %>% filter(cluster == 'Tail') %>% pull(fit.value)  %>% as.numeric()) * 100
  n.tail = ifelse('Tail' %in% names(clus.size), clus.size['Tail'], crayon::red('0'))
  shape.tail = x$Clusters %>%
    dplyr::filter(cluster == 'Tail', type == 'Shape') %>%
    dplyr::mutate(fit.value = formatC(fit.value, digits = 2)) %>% pull(fit.value)
  
  # Print components
  # cat(crayon::black(crayon::bgYellow("\n  Components (mixture)  \n")))
  # cat("\n")
  
  
  if (!x$fit.tail)
    cli::cli_alert_danger("No tail fit.")
  else
    cli::cli_li(
      paste(
        sprintf('%7s', 'Tail'),
        "[n = {.value {n.tail}}, {.value {pi.tail}}%] with alpha = {.value {shape.tail}}."
      )
    )
  
  # if (!x$fit.tail)
  #   cat(sprintf('%9s', 'Tail'), crayon::red('OFF\n'))
  # else
  #   cat(paste0(
  #     sprintf('%9s', 'Tail'),
  #     '\tn = ', n.tail, ' (', pi.tail  , ') \t Shape = ', shape.tail), '\n')
  
  B.comp = x$Clusters %>%
    dplyr::filter(cluster != 'Tail', type == 'Mean' |
                    type == 'Mixing proportion') %>%
    dplyr::select(-init.value) %>%
    spread(type, fit.value)
  
  B.comp$`Mixing proportion` = round(B.comp$`Mixing proportion`, 2) * 100
  B.comp$Mean = round(B.comp$Mean, 2)
  
  for (i in 1:nrow(B.comp))
    cli::cli_li(
      paste0(
        "Beta {.field {B.comp$cluster[i]}} [n = {.value {clus.size[B.comp$cluster[i]]}}, {.value {B.comp$`Mixing proportion`[i]}}%] with mean = {.value {B.comp$Mean[i]}}."
      )
    )
  
  # for (i in 1:nrow(B.comp))
  #   cat(
  #     paste0(
  #       sprintf('%9s', paste('Beta', B.comp$cluster[i])),
  #       ' \tn = ',
  #       ifelse(
  #         B.comp$cluster[i] %in% names(clus.size),
  #         clus.size[B.comp$cluster[i]],
  #         crayon::red('0')
  #       ),
  #       ' (',
  #       B.comp$`Mixing proportion`[i],
  #       ') \t Mean = ',
  #       B.comp$Mean[i],
  #       '\n'
  #     )
  #   )
  
  
  
  ####################### Scores
  # cat(crayon::black(crayon::bgYellow("\n  Scores\n")))
  # # cat('\n')
  # print(x$scores, row.names = '')
  
  cli::cli_end()
  
  # cat('\n')
  xs = round(x$scores, 2)
  cli::cli_alert_info(
    'Score(s): NLL = {.value {xs["NLL"]}}; ICL = {.value {xs["ICL"]}} ({.value {xs["reICL"]}}), H = {.value {xs["entropy"]}} ({.value {xs["reduced.entropy"]}}), '
  )
}

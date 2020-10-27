#' Summary for an object of class \code{'dbpmm'} is a print.
#'
#' @param object An obj of class \code{'dbpmm'}.
#' @param ...
#'
#' @return See \code{\link{print}}.
#' @exportS3Method summary dbpmm
#' @export summary.dbpmm
#'
#' @examples
#' data(fit_example)
#' summary(fit_example$best)
summary.dbpmm = function(object, ...) {
  print.dbpmm(object, ...)
}

#' Summaries for an object of class \code{'dbpmm'} is like a print.
#'
#' @param x An obj of class \code{'dbpmm'}.
#' @param ...
#'
#' @return nothing.
#' @exportS3Method print dbpmm
#' @export print.dbpmm
#' @importFrom crayon white red green yellow black bgYellow blue bold
#' @importFrom cli cli_rule cli_text
#' @importFrom clisymbols symbol
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
      'n = {.value {x$N}} with {.field k = {x$Kbeta}} Beta(s) {.value {ifelse(x$fit.tail, crayon::green("and a tail"), crayon::red("without tail"))}}' 
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
  
  # Pi proportions
  sor_p = sort(x$pi, decreasing = TRUE)
  sor_p = names(sor_p)
  
  pi = round(x$pi[sor_p], digits = 2) * 100
  pi = pi[pi > 0]
  pi_label = paste0(pi, '% [', yellow(names(pi)), ']')
  
  mobster:::m_txt(
    ' Clusters: \u03C0 = {.value {pi_label}}, with \u03C0 > 0.',
    symbol = 'clisymbols::symbol$bullet'
  ) %>% cli::cli_text()
  
  
  if (!x$fit.tail)
    cli::cli_text("{crayon::red(clisymbols::symbol$cross)} No tail fit.") %>% cli::cli_text()
  else
    mobster:::m_txt(
      paste(
        sprintf('%7s', 'Tail'),
        "[n = {.value {n.tail}}, {.value {pi.tail}}%] with alpha = {.value {shape.tail}}."
      ),
      symbol = 'clisymbols::symbol$bullet'
    ) %>% cli::cli_text()
  
  B.comp = x$Clusters %>%
    dplyr::filter(cluster != 'Tail', type == 'Mean' |
                    type == 'Mixing proportion') %>%
    dplyr::select(-init.value) %>%
    tidyr::spread(type, fit.value)
  
  B.comp$`Mixing proportion` = round(B.comp$`Mixing proportion`, 2) * 100
  B.comp$Mean = round(B.comp$Mean, 2)
  
  for (i in 1:nrow(B.comp))
    mobster:::m_txt(
      paste0(
        "Beta {.field {B.comp$cluster[i]}} [n = {.value {clus.size[B.comp$cluster[i]]}}, {.value {B.comp$`Mixing proportion`[i]}}%] with mean = {.value {B.comp$Mean[i]}}."
      ),
      symbol = 'clisymbols::symbol$bullet'
    ) %>% cli::cli_text()
  
  # cat('\n')
  xs = round(x$scores, 2)
  mobster:::m_inf('Score(s): NLL = {.value {xs["NLL"]}}; ICL = {.value {xs["ICL"]}} ({.value {xs["reICL"]}}), H = {.value {xs["entropy"]}} ({.value {xs["reduced.entropy"]}}). Fit {.value {ifelse(x$status, crayon::green("converged"), crayon::red("interrupted"))}} by {.field {x$fit.type}} in {.value {length(x$all.NLL)}} steps.') %>% cli::cli_text()
  
  if(has_drivers_annotated(x))
  {
    drv_list = x$data %>%
      dplyr::filter(is_driver) 
    
    if(nrow(drv_list) > 0)
    {
      mobster:::m_inf('The fit object model contains also drivers annotated.') %>% cli::cli_text()
      pio::pioDisp(drv_list)
    }
  }
  
}

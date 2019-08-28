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
#'
#' @examples
#' data(fit_example)
#' print(fit_example$best)
print.dbpmm = function(x, ...)
{
  stopifnot(inherits(x, "dbpmm"))
  
  cat(
    crayon::bgYellow(crayon::black("[ MOBSTER ]")),
    'N =',
    x$N,
    crayon::cyan("points with"),
    'K =',
    x$Kbeta,
    crayon::cyan("Beta clusters, fit by"),
    crayon::yellow(x$fit.type),
    crayon::cyan('in'),
    length(x$all.NLL),
    crayon::cyan("steps"),
    ifelse(
      x$status,
      crayon::green('[CONVERGED]'),
      crayon::red('[NON CONVERGED]')
    ),
    "\n"
  )
  
  ####################### Pi
  clus.size = table(x$data$cluster)
  clus.size = clus.size[order(clus.size)]
  
  clus.size.pi = x$Clusters %>%
    dplyr::select(-init.value) %>%
    dplyr::filter(type == 'Mixing proportion') %>%
    dplyr::mutate(fit.value = formatC(fit.value, digits = 2))
  
  
  pi.tail = clus.size.pi %>% filter(cluster == 'Tail') %>% pull(fit.value)
  n.tail = ifelse('Tail' %in% names(clus.size), clus.size['Tail'], crayon::red('0'))
  shape.tail = x$Clusters %>%
    dplyr::filter(cluster == 'Tail', type == 'Shape') %>%
    dplyr::mutate(fit.value = formatC(fit.value, digits = 2)) %>% pull(fit.value)
  
  
  # Print components
  cat(crayon::black(crayon::bgYellow("\n  Components (fit)  \n")))
  
  if (!x$fit.tail)
    cat(sprintf('%9s', 'Tail'), crayon::red('OFF\n'))
  else
    cat(paste0(
      sprintf('%9s', 'Tail'),
      '\tn = ', n.tail, ' (', pi.tail  , ') \t Shape = ', shape.tail), '\n')
  
  B.comp = x$Clusters %>%
    filter(cluster != 'Tail', type == 'Mean' |
             type == 'Mixing proportion') %>%
    select(-init.value) %>%
    spread(type, fit.value)

  B.comp$`Mixing proportion` = round(B.comp$`Mixing proportion`, 2)
  B.comp$Mean = round(B.comp$Mean, 2)
  
  for (i in 1:nrow(B.comp))
    cat(
      paste0(
        sprintf('%9s', paste('Beta', B.comp$cluster[i])),
        ' \tn = ',
        ifelse(
          B.comp$cluster[i] %in% names(clus.size),
          clus.size[B.comp$cluster[i]],
          crayon::red('0')
        ),
        ' (',
        B.comp$`Mixing proportion`[i],
        ') \t Mean = ',
        B.comp$Mean[i],
        '\n'
      )
    )
  
  
  
  ####################### Scores
  cat(crayon::black(crayon::bgYellow("\n  Scores (model selection)  \n")))
  print(x$scores, row.names = FALSE)
}

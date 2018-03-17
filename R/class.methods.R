#' Summaries for an object of class dbpmm is like a print.
#'
#' @param x the obj of class dbpmm
#' @param ...
#'
#' @return see print
#' @export
#'
#' @examples something..
summary.dbpmm = function(x, ...) {print.dbpmm(x, ...)}

#' Print an object of class dbpmm.
#'
#' @param x an object of class dbpmm..
#' @param ...
#'
#' @return nothing.
#' @export
#' @import crayon
#'
#' @examples something..
print.dbpmm = function(x, ...)
{
  cat(crayon::bgYellow(crayon::black("[ dbpmm ]")),
      'N =', x$N, crayon::cyan("samples with"),
      'K =', x$Kbeta, crayon::cyan("Beta clusters, fit by"),
      crayon::yellow(x$fit.type), crayon::cyan('in'), length(x$all.NLL), crayon::cyan("steps"),
      ifelse(x$status, crayon::green('[CONVERGED]'), crayon::red('[NON CONVERGED]')), "\n")

  cat(crayon::cyan("              Beta  -->  "))

  pars = .MeanVarBeta(x$beta['a', ], x$beta['b', ])
  pars = lapply(pars, round, digits = 2)

  greeks=c(alpha='\u03b1', tau='\u03c4', sigma='\u03c3',
           beta='\u03b2', gamma='\u03b3', mu ='\u00b5',
           pi = '\u03C0', propto = '\U221D')

  cat(crayon::yellow(greeks['mu'], '='), paste(pars$mean, collapse = ', '), ' ',
      crayon::yellow(greeks['sigma'], '='), paste(pars$var, collapse = ', '), '\n')

  cat(crayon::cyan("            Pareto  -->  "))
  if(all(is.na(x$tail))) cat(crayon::red('Off\n'))
  else cat(
    crayon::yellow(greeks['tau'], '='), x$tail$shape,
    crayon::yellow(greeks['gamma'], '='), x$tail$scale, '\n')

  cat(crayon::cyan("Mixing proportions  -->  "))
  pi = round(x$pi, 2)
  cat(crayon::yellow(greeks['pi'], '='), paste(pi, collapse = ', '), "\n")

  cat(crayon::cyan("Clusters dimension  -->  "))
  if(all(is.null(x$N.k))) cat(crayon::red('Unavaiable'), "\n")
  else cat(paste(x$N.k, collapse = ', '), "\n")

  cat(crayon::cyan("               NLL  --> "), crayon::yellow('D | \U03B8', greeks['propto']),  round(x$NLL, 2), "\n")

  if(!is.null(x$scores))
  {
    sc = round(x$scores, 2)
    sc = sc[, setdiff(colnames(sc), 'NLL')]

    cat(crayon::cyan("    Model selection -->  "))
    for(s in colnames(sc))
      cat(crayon::yellow(s, '='), sc[, s], ' ')

    cat('\n')
  }
}

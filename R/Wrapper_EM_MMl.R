#' Fit a model with MOBSTER.
#' 
#' @description This function fits the mixure of Beta distributions with a power-law Pareto
#' Type-I tail (optional). The function performs model selection for different mixtures, which
#' the user specify with the input parmeters. The function return a list of all fits computed
#' (objects of class \code{dbpmm}), the best fit, a table with the results of the fits and a
#' variable that specify which score has been used for model selection. 
#'
#' @param x Input tibble (or data.frame) which is required to have a VAF column which reports the 
#' frequency of the mutant allele (this should be computed adjusting the raw VAF for tumour purity 
#' and copy number status). See also package \code{mmobster} which contains functions to prepare
#' common inputs for this computation.
#' @param K A vector with the number of Beta components to use. All values of \code{K} must be positive 
#' and strictly greater than 0; they are combined with the value of \code{tail} to define all model
#' configurations tested for model selection
#' @param samples Number of fits that should be attempted for each configuration of the model tested.
#' @param init Initial values for the paremeters of the model. With \code{"ranodm"} the mean and variance 
#' of each Beta component are randomply sampled in the interval (0,1). With \code{"peaks"} a peak detection
#' heuristic is used to place the Beta means to match the peaks; in that case the variance is still randomised.
#' In both cases the power-law shape is randomised.
#' @param tail If \code{FALSE} the fit will not use a tail, if \code{TRUE} it will.
#' @param epsilon Tolerance for convergency estimation. For MLE fit this is compared to the differential of the 
#' negative log-likelihood (NLL); for MM fit the largest differential among the mixing proportions (pi) is used.
#' @param maxIter Maximum number of steps for a fit. If convergency is not achieved before these steps, the fit is interrupted.
#' @param is_verbose Verbose output. This also shows a sort of progress bar for the fitting.
#' @param fit.type A string that determines the type of fit. \code{"MLE"} is for the Maximum Likelihood Estimate 
#' of the Beta paraneters, \code{"MM"} is for the Momemnt Matching; MLE is numerical, and thus slower. In both cases
#' the estimator of the tail shape - if required - is its MLE and its analytical.
#' @param seed Seed for the random numbers generator
#' @param model.selection Score to minimize to select the best model; this has to be one of \code{'reICL'}, \code{'ICL'}, 
#' \code{'BIC'}, \code{'AIC'} or \code{'NLL'}. We advise to use only reICL and ICL
#' @param trace If \code{TRUE}, a trace of the model fit is returned that can be used to animate the model run.
#' @param parallel Optional parameter to run the fits in parallel (default), or not.
#'
#' @return A list of all fits computed (objects of class \code{dbpmm}), the best fit, a table with the results of the fits and a
#' variable that specify which score has been used for model selection.
#' @export
#'
#' @import crayon
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import magrittr
#'
#' @examples
#' TODO
mobster_fit = function(x,
                       K = 1:3,
                       samples = 5,
                       init = 'peaks',
                       tail = c(TRUE, FALSE),
                       epsilon = 1e-10,
                       maxIter = 2000,
                       fit.type = 'MM',
                       seed = 12345,
                       model.selection = 'reICL',
                       trace = FALSE,
                       parallel = TRUE)
{
  # Check for basic input requirements
  check_input(x, K, samples, init, tail, epsilon, maxIter, fit.type, seed, model.selection, trace)

  X = tibble::as.tibble(x)

  ###################### Initializations
  set.seed(seed)
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  # Storage variables
  best = obj = runs = NULL

  # Configurations that will be used for model selection
  tests = expand.grid(K = K, Run = 1:samples, tail = tail, stringsAsFactors = FALSE)
  tests = tests[order(tests$tail, tests$K), ]

  ntests = nrow(tests)

  ###################### Print headers
  pio::pioHdr("MOBSTER fit")

  cat(
    cyan('\n\t-'),
    'N =',
    nrow(x),
    cyan('samples with'),
    'K =',
    K,
    cyan('Beta; Pareto Type-I power-law tail:'),
    ifelse(tail, green('YES'), red('NO'))
  )
  cat(
    cyan('\n\t- Fit   : '),
    yellow(fit.type),
    cyan('for'),
    maxIter,
    cyan('max. steps with'),
    ifelse(all(is.character(init)), init, 'custom'),
    cyan('initialization;'),
    yellow('\u03B5 ='),
    epsilon,
    cyan('. Model selection: '),
    yellow(model.selection)
  )
  cat(
    cyan('\n\t- Runs  : '),
    samples,
    ' x ',
    length(K),
    ' x ',
    length(tail),
    '=',
    yellow(ntests),
    '\n'
  )

  # Input Checks
  # if (any(x$VAF == 0)) {
  #   cat(crayon::red('\n[VAFs = 0] setting them to 1e-9 to avoid numerical errors'))
  # 
  #   X$VAF[X$VAF == 0] = 1e-9
  # }
  # 
  # if (any(X$VAF == 1)) {
  #   cat(crayon::red('\n[VAFs = 1] setting them to 1-1e-9 to avoid numerical errors'))
  #   X$VAF[X$VAF == 1] = 1 - 1e-9
  # }
  # END: Input Checks

  
  # Inputs in the easypar format - list of lists
  inputs = lapply(1:nrow(tests), 
                  function(r)
                  list(
                    X = x,
                    K = tests[r, 'K'],
                    init = init,
                    tail = tests[r, 'tail'],
                    epsilon = epsilon,
                    maxIter = maxIter,
                    fit.type = fit.type,
                    trace = trace
                    ))
  
  # Fits are obtained using the easypar package
  # which allows easy parallelization of R functions
  #
  # https://github.com/caravagn/easypar
  #
  runs = easypar::run(
    FUN = .dbpmm.EM,
    PARAMS = inputs,
    packages = c("crayon", "pio", "tidyverse"),
    export = ls(globalenv(), all.names = TRUE),
    cores.ratio = .8,
    parallel = parallel,
    cache = NULL
  )
  
  # Polish errors if any
  nerrs = easypar::numErrors(runs)
  if(nerrs == samples) {
    
    lapply(runs, function(w) print(w$message))
    
    stop("All task returned errors, no fit available, raising error.")
  }
  
  runs = easypar::filterErrors(runs)
  
  # timing
  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")
  cat(bold("\n\nCOMPLETED."), round(TIME, 2), cyan('mins'), '\n')

  # Get all scores
  tests = cbind(tests, 
                Reduce(rbind, lapply(runs, function(w)
    w$scores)))

  # clean up some repeated results -- show unique fits
  tests$Run = NULL
  scores.columns = colnames(runs[[1]]$scores)
  tests[, scores.columns] = apply(tests[, scores.columns], 2, round, digits = 2)


  runs = runs[!duplicated(tests)] # remove duplicated entries..
  tests = tests[!duplicated(tests), ] # remove duplicated entries..
  rownames(tests) = NULL

  # Model selection -- this will be returned later ..
  model = model_selection(runs, scores.suitable = model.selection)
  model = model$model.selection[[model.selection]]
  model$model.selection = model.selection

  ###### SHOW BEST FIT
  cat("\n", bold("BEST FIT:", model.selection, "\n\n"))
  print.dbpmm(model$best)

  return(model)
}


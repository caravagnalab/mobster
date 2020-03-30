#' Fit a model with MOBSTER.
#'
#' @description This function fits the mixure of Beta distributions with a power-law Pareto
#' Type-I tail (optional). The function performs model selection for different mixtures, which
#' the user specify with the input parmeters. The function return a list of all fits computed
#' (objects of class \code{dbpmm}), the best fit, a table with the results of the fits and a
#' variable that specify which score has been used for model selection. The fitting of each model
#' also runs the function \code{choose_clusters} which implements a simple heuristic to filter
#' out small clusters from the fit output.
#'
#' Note: You can also use the "auto setup" functionality that, with one keyword, loads preset
#' parameter values in order to implement different analysis.
#'
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
#' @param fit.type A string that determines the type of fit. \code{"MLE"} is for the Maximum Likelihood Estimate
#' of the Beta paraneters, \code{"MM"} is for the Momemnt Matching; MLE is numerical, and thus slower. In both cases
#' the estimator of the tail shape - if required - is its MLE and its analytical.
#' @param seed Seed for the random numbers generator
#' @param model.selection Score to minimize to select the best model; this has to be one of \code{'reICL'}, \code{'ICL'},
#' \code{'BIC'}, \code{'AIC'} or \code{'NLL'}. We advise to use only reICL and ICL
#' @param trace If \code{TRUE}, a trace of the model fit is returned that can be used to animate the model run.
#' @param parallel Optional parameter to run the fits in parallel (default), or not.
#' @param pi_cutoff Parameter passed to function \code{choose_clusters}, which determines the minimum mixing proportion of a
#' cluster to be returned as output.
#' @param N_cutoff Parameter passed to function \code{choose_clusters}, which determines the minimum number of mutations
#' assigned to a cluster to be returned as output.
#' @param auto_setup Overrides all the parameters with an predined set of values, in order to implement different analyses.
#' Availables keys: `FAST`, uses 1) max 2 clones (1 subclone), 2) random initial conditions 3) 2 samples per parameter set
#' 4) mild `epsilon` and `maxIter`, sequential run. For reference, the default set of parameters represent a more exhaustive
#' analysis.
#' @param description A textual description of this dataset.
#'
#' @return A list of all fits computed (objects of class \code{dbpmm}), the best fit, a table with the results of the fits and a
#' variable that specify which score has been used for model selection.
#'
#' @export
#'
#' @import crayon
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import magrittr
#'
#' @examples
#' # Generate a random dataset
#' x = random_dataset(seed = 123, Beta_variance_scaling = 100, N = 200)
#' print(x) # Contains a ggplot object
#'
#' # Fit, default models, changed epsilon for convergence
#' x = mobster_fit(x$data, epsilon = 1e-5)
#'
#' plot(x$best)
#' print(x$best)
#'
#' lapply(x$runs[1:3], plot)
mobster_fit = function(x,
                       K = 1:3,
                       samples = 5,
                       init = 'peaks',
                       tail = c(TRUE, FALSE),
                       epsilon = 1e-10,
                       maxIter = 250,
                       fit.type = 'MM',
                       seed = 12345,
                       model.selection = 'reICL',
                       trace = FALSE,
                       parallel = TRUE,
                       pi_cutoff = 0.02,
                       N_cutoff = 10,
                       auto_setup = NULL,
                       silent = FALSE,
                       description = "My MOBSTER model")
{
  pio::pioHdr(paste0("MOBSTER fit"))
  cat('\n')
  
  # Check for basic input requirements
  check_input(x,
              K,
              samples,
              init,
              tail,
              epsilon,
              maxIter,
              fit.type,
              seed,
              model.selection,
              trace)
  
  X = tibble::as_tibble(x)
  
  mobster:::m_ok("Loaded input data, n = {.value {nrow(x)}}.") %>% cli::cli_text()
  
  ###################### Auto setup of parameters
  if (!is.null(auto_setup))
  {
    # Get the parameters, checks they are known, throws errors.
    template = auto_setup(auto_setup)
    
    K = template$K
    samples = template$samples
    init = template$init
    tail = template$tail
    epsilon = template$epsilon
    maxIter = template$maxIter
    fit.type = template$fit.type
    seed = template$seed
    model.selection = template$model.selection
    trace = template$trace
    parallel = template$parallel
    pi_cutoff = template$pi_cutoff
    N_cutoff = template$N_cutoff
  }
  
  ###################### Initializations
  set.seed(seed)
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  
  # Storage variables
  best = obj = runs = NULL
  
  # Configurations that will be used for model selection
  tests = expand.grid(
    K = K,
    tail = tail,
    Run = 1:samples,
    stringsAsFactors = FALSE
  )
  ntests = nrow(tests)
  
  ###################### Print message
  mobster:::m_txt(
    "n = {.value {nrow(x)}}. Mixture with k = {.field {paste(K, collapse = ',')}} Beta(s). Pareto tail: {.field {tail}}. Output clusters with \u03c0 > {.value {pi_cutoff}} and n > {.value {N_cutoff}}."
  ) %>% cli::cli_text()
  
  if (!is.null(auto_setup))
    mobster:::m_wrn("mobster automatic setup {.field {auto_setup}} for the analysis.") %>% cli::cli_text()
  else
    mobster:::m_txt(
      'Custom fit by {.field {ifelse(fit.type == \'MM\', "Moments-matching", "Maximum-Likelihood")}} in up to {.value {maxIter}} steps, with \u03B5 = {.value {epsilon}} and {.field {init}} initialisation.'
    ) %>% cli::cli_text()
  
  mobster:::m_txt(
    'Scoring ({.value {ifelse(parallel, green("with parallel"), red("without parallel"))}}) {.value {samples}} x {.value {length(K)}} x {.value {length(tail)}} = {.field {ntests}} models by {.field {model.selection}}.'
  ) %>% cli::cli_text()
  cat('\n')
  
  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Fits are obtained using the easypar package
  # which allows easy parallelization of R functions
  #
  # https://github.com/caravagn/easypar
  #
  # Inputs in the easypar format - list of lists
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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
                      trace = trace,
                      pi_cutoff = pi_cutoff,
                      N_cutoff = N_cutoff, 
                      description = description
                    ))
  
  runs = easypar::run(
    FUN = .dbpmm.EM,
    PARAMS = inputs,
    # packages = c("crayon", "pio", "tidyverse"),
    export = ls(globalenv(), all.names = TRUE),
    cores.ratio = .8,
    parallel = parallel,
    cache = NULL,
    filter_errors = TRUE # Error moanagment is inside easypar
  )
  
  # Should not happen
  if (length(runs) == 0)
    stop("All task returned errors, no fit available, raising this error to interrupt the computation....")
  
  # What is successfull (id of the task)
  succesfull_tasks = names(runs) %>% as.numeric()
  tests = tests[succesfull_tasks, , drop = FALSE]
  
  # Report timing to screen
  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")
  
  cat('\n\n')
  mobster:::m_inf("{crayon::bold('MOBSTER fits')} completed in {.value {prettyunits::pretty_dt(TIME)}}.") %>% cli::cli_text()
  cat('\n')
  
  # Get all scores
  scores_succesfull_tasks = lapply(runs, function(w)
    w$scores)
  tests = bind_cols(tests, Reduce(bind_rows, scores_succesfull_tasks))
  
  # clean up some repeated results -- show unique fits
  tests$Run = NULL
  scores.columns = colnames(runs[[1]]$scores)
  tests[, scores.columns] = apply(tests[, scores.columns], 2, round, digits = 2)
  
  runs = runs[!duplicated(tests)] # remove duplicated entries..
  tests = tests[!duplicated(tests),] # remove duplicated entries..
  rownames(tests) = NULL
  
  # Model selection -- this will be returned later ..
  model = model_selection(runs, scores.suitable = model.selection)
  model = model$model.selection[[model.selection]]
  model$model.selection = model.selection
  
  ###### SHOW BEST FIT
  # cli::cli_alert_info(paste(bold("BEST:"), model.selection))
  print.dbpmm(model$best)
  
  return(model)
}

#' Fit a model with MOBSTERH.
#'
#' @description This function fits a bayesian hierarchical version of the MOBSTER model implemented in \code{mobster_fit}. We still
#' have a mixture of beta distributions and an optional Pareto type-one distribution to model the neutral tail. From a modelling
#' point of view the main difference here is that we are expanding that model over different karyotypes and we treat the problem from a bayesian point of view.
#' In this way we can grant information about the mutationt rate and the tail pooling from the different karyotypes and at the same time include the
#' strong prior knowledge we have about
#'
#'
#' @param x Input tibble (or data.frame) or an evopipe_qc object (preferred),
#' @param subclonal_clusters A vector with the number of Beta components to use. All values of \code{K} must be positive
#' and strictly greater than 0; they are combined with the value of \code{tail} to define all model
#' configurations tested for model selection
#' @param init Initial values for the paremeters of the model. With \code{"ranodm"} the mean and variance
#' of each Beta component are randomply sampled in the interval (0,1). With \code{"peaks"} a peak detection
#' heuristic is used to place the Beta means to match the peaks; in that case the variance is still randomised.
#' In both cases the power-law shape is randomised.
#' @param tail If \code{0} the fit will not use a tail, if \code{1} it will.
#' @param epsilon Tolerance for convergency estimation. For MLE fit this is compared to the differential of the
#' negative log-likelihood (NLL); for MM fit the largest differential among the mixing proportions (pi) is used.
#' @param maxIter Maximum number of steps for a fit. If convergency is not achieved before these steps, the fit is interrupted.
#' @param seed Seed for the random numbers generator
#' @param model.selection Score to minimize to select the best model; this has to be one of \code{'ICL'},
#' \code{'BIC'}, \code{'AIC'} or \code{'likelihood'}. We advise to use only reICL and ICL
#' @param parallel Optional parameter to run the fits in parallel (default), or not.
#' @param pi_cutoff Parameter passed to function \code{choose_clusters}, which determines the minimum mixing proportion of a
#' cluster to be returned as output.
#' @param N_cutoff Parameter passed to function \code{choose_clusters}, which determines the minimum number of mutations
#' assigned to a cluster to be returned as output.
#' @param description A textual description of this dataset.
#' @param lrd_gamma learning rate decay fator, final learning rate is gonna be lrd_gamma * lr
#' @param vaf_filter Discard mutations under a specific VAF threshold for the fitting procedure
#' @param n_t Discard karyotypes with less then a given number of mutations.
#' @return A list of all fits computed (objects of class \code{dbpmm}), the best fit, a table with the results of the fits and a
#' variable that specify which score has been used for model selection.
#'
#' @importFrom dplyr filter mutate select arrange desc pull row_number group_by
#' @importFrom dplyr summarise bind_cols rename bind_rows left_join distinct
#' @importFrom dplyr ungroup full_join
#' @importFrom tidyr spread gather tibble tribble as_tibble unnest
#' @importFrom magrittr %>%
#' @importFrom tibble enframe
#' @export
#
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
#'
mobsterh_fit = function(x,
                       subclonal_clusters = 1:2,
                       tail = c(TRUE, FALSE),
                       truncate_pareto= c(TRUE, FALSE),
                       purity = 1.,
                       samples = 1,
                       epsilon = 1e-5,
                       maxIter = 2000,
                       model.selection = 'ICL',
                       parallel = FALSE,
                       pi_cutoff = 0.05,
                       N_cutoff = 80,
                       silent = FALSE,
                       alpha_precision_concentration = 5,
                       alpha_precision_rate = 0.01,
                       number_of_trials_clonal_mean=300,
                       number_of_trials_k = 150,
                       prior_lims_clonal=c(0.1, 100000.),
                       prior_lims_k=c(0.1, 100000.),
                       lr = 0.05,
                       compile = FALSE,
                       CUDA = FALSE,
                       description = "My MOBSTERH model",
                       karyotypes = c("1:0", "1:1","2:1", "2:0", "2:2"),
                       lrd_gamma = 0.1,
                       vaf_filter = 0.05,
                       n_t = 100)
{
  pio::pioHdr(paste0("MOBSTERh fit"))
  cat('\n')


  data_raw <- x

  if(inherits(x, "evopipe_qc")){
    purity <- get_purity(x)
    x <- format_data_mobsterh_QC(x, vaf_t = vaf_filter, n_t = n_t)
  } else if (is.matrix(x) | is.data.frame(x)){
    if(!all(colnames(x) %in% c("VAF","karyotypes", "chr", "from", "to")))
      stop("Please provide a data.frame with the following columns: VAF, karyotypes, chr, from, to")
    x <- format_data_mobsterh_DF(x)
  } else {
    stop("Please provide a data.frame or a CNAqc object as input")
  }


  # Check for basic input requirements
  mobster:::check_inputh(
                        K,
                        subclonal_clusters,
                        tail,
                        epsilon,
                        maxIter,
                        fit.type,
                        lr,
                        model.selection,
                        karyotypes,
                        prior_lims_k,
                        prior_lims_clonal,
                        alpha_precision_concentration,
                        alpha_precision_rate
                        )


  mobster:::m_ok("Loaded input data, {.value {length(x)}} karyotypes.") %>% cli::cli_text()
  for(k in 1:length(x)){
    cli::cli_alert_info("Karyotype {names(x)[k]} with n = {.value {length(x[[k]])}} mutations")
  }

  ###################### Initializations
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  # Storage variables
  best = obj = runs = model = NULL

  # Configurations that will be used for model selection
  tests = expand.grid(
    subclonal_clusters = subclonal_clusters,
    tail = tail,
    truncate_pareto = truncate_pareto,
    stringsAsFactors = FALSE
  )
  ntests = nrow(tests)

  ###################### Print message

  mobster:::m_txt(
    "n = {.value {nrow(x)}}. Mixture with S = {.field {paste(subclonal_clusters, collapse = ',')}} subclonal Beta(s). Pareto tail: {.field {tail}}. Output clusters with \u03c0 > {.value {pi_cutoff}} and n > {.value {N_cutoff}}."
  ) %>% cli::cli_text()


  mobster:::m_txt(
    'Custom fit by Gradient-Based Stocastic Variational Inference in up to {.value {maxIter}} steps, with \u03B5 = {.value {epsilon}}..'
  ) %>% cli::cli_text()

  mobster:::m_txt(
    'Scoring ({.value {ifelse(parallel, green("with parallel"), red("without parallel"))}}) {.value {length(subclonal_clusters)}} x {.value {length(tail)}} = {.field {ntests}} models by {.field {model.selection}}.'
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
                      data = x,
                      data_raw = data_raw,
                      subclonal_clusters =  as.integer(tests[r, 'subclonal_clusters']),
                      tail = as.integer(tests[r, 'tail']),
                      truncate_pareto = tests[r, 'truncate_pareto'],
                      samples = samples,
                      purity = purity,
                      max_it = as.integer(maxIter),
                      alpha_precision_concentration = alpha_precision_concentration,
                      alpha_precision_rate = alpha_precision_rate,
                      number_of_trials_clonal_mean=number_of_trials_clonal_mean,
                      number_of_trials_k = number_of_trials_k,
                      prior_lims_clonal=prior_lims_clonal,
                      prior_lims_k=prior_lims_k,
                      lr = lr,
                      e = epsilon,
                      compile = compile,
                      CUDA = CUDA,
                      description = description,
                      lrd_gamma = lrd_gamma
                    ))


  runs = easypar::run(
    FUN = mobsterh_fit_aux,
    PARAMS = inputs,
    export = ls(globalenv(), all.names = TRUE),
    cores.ratio = .5,
    parallel = parallel,
    cache = NULL,
    filter_errors = TRUE # Error managment is inside easypar
    ,progress_bar = FALSE
  )



  # Report timing to screen
  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  cat('\n\n')
  mobster:::m_inf("{crayon::bold('MOBSTERH fits')} completed in {.value {prettyunits::pretty_dt(TIME)}}.") %>% cli::cli_text()
  cat('\n')

  # Get all scores
  scores_succesfull_tasks = lapply(runs, function(w)
    w$information_criteria  %>% as.data.frame())
  tests = bind_cols(tests, Reduce(bind_rows, scores_succesfull_tasks))


  # Model selection -- this will be returned later ..

  tests$likelihood <- -tests$likelihood

  sort_runs <- order(tests[, model.selection])
  runs <- runs[sort_runs]
  model$best = runs[[1]]
  model$model.selection = model.selection

  model$runs <-  runs
  model$fits.table <- tests

  ###### SHOW BEST FIT
  print.dbpmmh(model$best)

  return(model)
}




mobsterh_fit_aux <-  function( data,
                               data_raw,
                              subclonal_clusters,
                              tail,
                              truncate_pareto,
                              samples,
                              purity,
                              max_it,
                              seed,
                              alpha_precision_concentration,
                              alpha_precision_rate,
                              number_of_trials_clonal_mean,
                              number_of_trials_k,
                              prior_lims_clonal,
                              prior_lims_k,
                              lr,
                              e,
                              compile,
                              CUDA,
                              description,
                              lrd_gamma){


  data_u <- data
  data <- tensorize(data_u)

  mob <- reticulate::import("mobster")

  inf_res <- lapply(1:samples,function(s) mob$fit_mobster(data = data, K = as.integer(subclonal_clusters),
                                               tail = as.integer(tail),
                                               truncated_pareto = truncate_pareto,
                                               purity = purity,
                                               max_it = max_it,
                                               lr = lr,
                                               e = e,
                                               alpha_precision_concentration = alpha_precision_concentration,
                                               alpha_precision_rate = alpha_precision_rate,
                                               number_of_trials_clonal_mean = number_of_trials_clonal_mean,
                                               number_of_trials_k = number_of_trials_k,
                                               prior_lims_clonal = prior_lims_clonal
                                               ,prior_lims_k = prior_lims_k
                                               ,compile = compile, CUDA = CUDA, lrd_gamma = lrd_gamma))

  best_model <- which.max(sapply(inf_res, function(h) h$information_criteria$likelihood))
  inf_res <-  inf_res[[best_model]]

  assig_temp = lapply(1:length(data_u), function(i) return(data.frame(id = names(data_u[[i]]),
                                                                      cluster = inf_res$model_parameters[[i]]$cluster_assignments
                                                                      )))
  assig_temp = Reduce(assig_temp, f = rbind)
  if(inherits(data_raw, what = "evopipe_qc"))
    table = data_raw$cnaqc$snvs %>% select(chr, from, to, ref, alt, VAF,karyotype, is_driver, driver_label) %>%
      mutate(id = paste(chr,from, to, sep = ":"))
  else
    table = data_raw  %>%
    mutate(id = paste(chr,from, to, sep = ":"))

  inf_res$data = dplyr::left_join(table, assig_temp, by = "id", copy = T) %>% as.data.frame()

  inf_res$description <- description

  class(inf_res) <- "dbpmmh"

  return(inf_res)

}

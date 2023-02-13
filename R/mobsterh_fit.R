#' Fit a model with MOBSTERh.
#'
#' @description This function fits a bayesian hierarchical version of the MOBSTER model implemented in \code{mobster_fit}. We still
#' have a mixture of Beta distributions and an optional Pareto type-one distribution to model the neutral tail. From a modelling
#' point of view there are two main differenceshere:
#'   1. We are expanding that model over different karyotypes and we treat the problem from a bayesian point of view.
#'   2. We model do not consider the VAF deconvolution but the read count deconvolution (with a beta-binomial noise term)
#' In this way we can grant information about the mutation rate and the tail pooling from the different karyotypes and at the same time include the
#' strong prior knowledge we have about how clonal and subclonal clusters are supposed to be ddistributed along the VAF spectrum.
#' By modelling counts then, we also explicitly account for the binomial (with dispersion) sampling process that happens
#' during sequencing.
#'
#' All the Beta distributions in the prior model are not modelled using concentration parameters but using this parametrization:
#' \deqn{concentration1 = mean  * number_of_trials}
#' \deqn{concentration2 = (1 - mean)  * number_of_trials}
#'
#'
#' @param x Input tibble (or data.frame), cnaqc object ot evopipe_qc object (laste two preferred),
#' the input data.frame should have at least 6 coloumns named as chr, from, to, NV (number of variant), DP (depth) and karyotype.
#' @param subclonal_clusters A vector with the number of Beta components to use. All values of \code{K} must be positive
#' and strictly greater than 0; they are combined with the value of \code{tail} and to define all model
#' configurations tested for model selection
#' @param truncate_pareto When \code{TRUE} the fit will be done with a Truncated pareto distribution with probability density equal to 0
#' for values x greater than the mean of the smallest clonal cluster
#' @param tail If \code{TRUE} the fit will not use a Pareto to model the tail, if \code{FALSE} it will.
#' @param purity User provided tumor purity, used only when the input is a data.frame
#' @param samples If the number of samples is greater than 1, then for each tail-truncation-subclone configuration \code{samples} fit
#' are produced and the one with the highest \code{model.selection} values is taken.
#' @param enforce_QC_PASS if \code{TRUE} when using a cnaqc of evopipe_qc object fit just the karyotype that passe QC.
#' @param epsilon Tolerance for convergency estimation. As ELBO oscillations are common in gradient based VI, we will monitor the convergence of all the parameters in the model,
#' the inference stops when (abs(new-old) / abs(old)) < epsilon, for all the parameters
#' @param maxIter Maximum number of steps for a fit. If convergency is not achieved before these steps, the fit is interrupted.
#' @param model.selection Score to minimize to select the best model; this has to be one of \code{'ICL'},
#' \code{'BIC'}, \code{'AIC'} or \code{'likelihood'}. We advise to use only ICL
#' @param parallel Optional parameter to run the fits in parallel, or not (default).
#' @param alpha_precision_concentration Concentration value for the gamma modelling the prior shape of the Pareto
#' @param alpha_precision_rate Rate value for the gamma modelling the prior shape of the Pareto
#' @param number_of_trials_clonal_mean Number of trials for the Beta prior over the clonal clusters mean
#' @param number_of_trials_k Number of trials for the Beta prior over the subclonal clusters mean
#' @param prior_lims_clonal Bounds on the uniform prior over the number of trials for the clonal clusters
#' @param prior_lims_k Bounds on the uniform prior over the number of trials for the subclonal clusters
#' @param lr Learning rate used by the Adam oprimizer
#' @param compile Use the just-in-time (JIT) compiler
#' @param CUDA Use the default GPU to train the model (you need to setup PyTorch for this)
#' @param description A textual description of this dataset.
#' @param lrd_gamma learning rate decay fator, final learning rate is gonna be lrd_gamma * lr
#' @param vaf_filter Discard mutations under a specific VAF threshold for the fitting procedure
#' @param n_t Discard karyotypes with less then a given number of mutations.
#' @param quantile_filt Filter the mutations with VAF higher than those quantile
#' @param N_MAX subsample an N_MAX number of mutations, it keeps the drivers. Works only with a CNAqc input
#' @param assign_drivers assign drivers in non used karyotype using betabinomials and pareto-binomials mixtures
#' @param assign_mutation_posteriori assign mutation in used karyotypes not included in the main fit
#'
#' @return An object of class \code{mobster_deconv}, i.e. list of all fits computed (objects of class \code{dbpmm}), the best fit, a table with the results of the fits and a
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
#' data("fit_example_mobsterh", package = "mobster")
#'
#'x = fit_example_mobsterh$best$data
#'
#' # Fit, default model
#' x = mobster_fit(x, K = 0:1, truncate_pareto = FALSE)
#'
#' plot(x$best)
#' print(x$best)
#'
mobsterh_fit = function(x,
                        subclonal_clusters = 0:2,
                        tail = c(TRUE, FALSE),
                        truncate_pareto = TRUE,
                        subclonal_prior = "Beta",
                        multi_tail = FALSE,
                        purity = 1.,
                        samples = 1,
                        k_means_init = TRUE,
                        enforce_QC_PASS = TRUE,
                        epsilon = 0.005,
                        maxIter = 2000,
                        model.selection = 'ICL',
                        parallel = FALSE,
                        alpha_precision_concentration = 5,
                        alpha_precision_rate = 0.1,
                        number_of_trials_clonal_mean = 1000,
                        number_of_trials_k = 500,
                        number_of_trials_subclonal = 900,
                        prior_lims_clonal = c(0.1, 100000.),
                        prior_lims_k = c(0.1, 100000.),
                        lr = 0.01,
                        compile = FALSE,
                        CUDA = FALSE,
                        description = "My MOBSTERh model",
                        karyotypes = c("1:0", "1:1", "2:1", "2:0", "2:2"),
                        lrd_gamma = 0.1,
                        vaf_filter = 0.05,
                        NV_filter = 5,
                        n_t = 100,
                        quantile_filt = 1,
                        N_MAX = 50000,
                        assign_drivers = TRUE,
                        assign_mutation_posteriori = FALSE)
{
  pio::pioHdr(paste0("MOBSTERh fit"))
  cat('\n')

  data_raw <- NULL
  can_work = FALSE

  # Evoverse pipeline input
  if (inherits(x, "evopipe_qc"))
  {
    purity <- mobster:::get_purity(x)
    x$cnaqc <- CNAqc::subsample(x$cnaqc,N = N_MAX)
    data_raw <- x
    x <-
      format_data_mobsterh_QC(x,
                              vaf_t = vaf_filter,
                              n_t = n_t,
                              enforce_QC_PASS = enforce_QC_PASS,
                              NV_filter = NV_filter
                              )

    can_work = TRUE
  }

  # CNAqc object
  if (inherits(x, "cnaqc"))
  {
    x <- CNAqc::subsample(x, N = N_MAX)
    purity <- x$purity
    data_raw <- x$mutations
    x <- format_data_mobsterh_QC(x,
                              vaf_t = vaf_filter,
                              n_t = n_t,
                              enforce_QC_PASS = enforce_QC_PASS,
                              NV_filter = NV_filter
                              )

    can_work = TRUE
  }

  # Tibble
  if (is.matrix(x) | is.data.frame(x))
  {
    if (!all(c("NV", "DP", "karyotype", "chr", "from", "to") %in% colnames(x)))
      stop(
        "Please provide a data.frame with the following columns: NV, DP, karyotype, chr, from, to"
      )

    data_raw <- x
    x <- mobster:::format_data_mobsterh_DF(x,
                                 vaf_t = vaf_filter,
                                 NV_filter = NV_filter,
                                 n_t = n_t)
    cli::cli_alert_warning("Using input purity {.field {purity}}")

    can_work = TRUE
  }

  if (!can_work) {
    stop(
      "Input must be any of:\n\t- a data.frame, \n\t- an evoverse data QC pipeline output, \n\t- a CNAqc object."
    )
  }

  if (is.null(x))
    return(NULL)

  x <-  lapply(x, function(k) {
    qf <- quantile((k[,1] / k[,2]) , quantile_filt)
    k[(k[,1] / k[,2]) <= qf, ]
  })

  
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
  for (k in 1:length(x)) {
    cli::cli_alert_info("Karyotype {names(x)[k]} with n = {.value {length(x[[k]])}} mutations")
  }

  ###################### Initializations
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  # Storage variables
  best = obj = runs = model = NULL

  # Configurations that will be used for model selection
  tests = expand.grid(
    subclonal_clusters = subclonal_clusters,
    subclonal_prior = subclonal_prior,
    multi_tail = multi_tail,
    tail = tail,
    truncate_pareto = truncate_pareto,
    stringsAsFactors = FALSE
  )
  if(length(truncate_pareto) == 2)
    tests <- tests[-which(!tests$tail & tests$truncate_pareto), ]
  ntests = nrow(tests)

  ###################### Print message

  mobster:::m_txt(
    "n = {.value {nrow(tests)}}. Mixture with S = {.field {paste(subclonal_clusters, collapse = ',')}} subclonal Beta(s). Pareto tail: {.field {tail}}."
  ) %>% cli::cli_text()


  mobster:::m_txt(
    'Custom fit by Gradient-Based Stocastic Variational Inference in up to {.value {maxIter}} steps, with \u03B5 = {.value {epsilon}}..'
  ) %>% cli::cli_text()

  mobster:::m_txt(
    'Scoring ({.value {ifelse(parallel, green("with parallel"), red("without parallel"))}}) {.value {.field {ntests}}} models by {.field {model.selection}}.'
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
                      subclonal_prior = tests[r, 'subclonal_prior'],
                      multi_tail = tests[r, 'multi_tail'],
                      samples = samples,
                      purity = purity,
                      max_it = as.integer(maxIter),
                      alpha_precision_concentration = alpha_precision_concentration,
                      alpha_precision_rate = alpha_precision_rate,
                      number_of_trials_clonal_mean = number_of_trials_clonal_mean,
                      number_of_trials_k = number_of_trials_k,
                      prior_lims_clonal = prior_lims_clonal,
                      prior_lims_k = prior_lims_k,
                      lr = lr,
                      e = epsilon,
                      compile = compile,
                      CUDA = CUDA,
                      description = description,
                      lrd_gamma = lrd_gamma,
                      number_of_trials_subclonal = number_of_trials_subclonal,
                      k_means_init = k_means_init
                    ))


  runs = easypar::run(
    FUN = mobster:::mobsterh_fit_aux,
    PARAMS = inputs,
    export = ls(globalenv(), all.names = TRUE),
    cores.ratio = .5,
    parallel = parallel,
    cache = NULL,
    filter_errors = FALSE,
    progress_bar = FALSE,silent = F
  )

  filt <-  sapply(runs,function(w) !is.null(w$information_criteria))
  runs <-  runs[filt]



  # Report timing to screen
  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  cat('\n\n')
  mobster:::m_inf(
    "{crayon::bold('MOBSTERh fits')} completed in {.value {prettyunits::pretty_dt(TIME)}}."
  ) %>% cli::cli_text()
  cat('\n')

  # Get all scores
  scores_succesfull_tasks = lapply(runs, function(w)
    w$information_criteria  %>% as.data.frame())
  tests = bind_cols(tests[filt,] , Reduce(bind_rows, scores_succesfull_tasks))


  # Model selection -- this will be returned later ..

  tests$likelihood <- -tests$likelihood

  sort_runs <- order(tests[, model.selection])
  runs <- runs[sort_runs]
  model$best = runs[[1]]
  model$model.selection = model.selection

  model$runs <-  runs
  model$fits.table <- tests
  

  if(assign_drivers)
    model$best <- mobster:::assign_drivers(model$best)
  
  if(assign_mutation_posteriori)
    model$best <- mobster:::assign_mutations_posteriori(model$best)


  ###### SHOW BEST FIT
  print.dbpmmh(model$best)

  # Add S3 class
  class(model) <- "mobster_deconv"

  return(model)
}




mobsterh_fit_aux <-  function(data,
                              data_raw,
                              subclonal_clusters,
                              tail,
                              truncate_pareto,
                              subclonal_prior,
                              multi_tail,
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
                              lrd_gamma,
                              number_of_trials_subclonal,
                              k_means_init) {
  data_u <- data
  
  used_muts <- lapply(data_u, rownames) %>% do.call(c,.) %>% unname()
  
  data <- mobster:::tensorize(data_u)

  mob <- reticulate::import("mobster")

  inf_res <-
    lapply(1:samples, function(s)
      mob$fit_mobster(
        data = data,
        K = as.integer(subclonal_clusters),
        tail = as.integer(tail),
        truncated_pareto = truncate_pareto,
        subclonal_prior = subclonal_prior,
        purity = purity,
        max_it = max_it,
        lr = lr,
        e = e,
        alpha_precision_concentration = alpha_precision_concentration,
        alpha_precision_rate = alpha_precision_rate,
        number_of_trials_clonal_mean = number_of_trials_clonal_mean,
        number_of_trials_k = number_of_trials_k,
        prior_lims_clonal = prior_lims_clonal,
        prior_lims_k = prior_lims_k,
        compile = compile,
        CUDA = CUDA,
        lrd_gamma = lrd_gamma,
        k_means_init = k_means_init
      ))

  best_model <-
    which.max(sapply(inf_res, function(h)
      h$information_criteria$likelihood))
  inf_res <-  inf_res[[best_model]]

  for(i in 1:length(data_u)) {
    names(inf_res$model_parameters[[i]]$cluster_assignments) <- names(data_u[[i]])
    colnames(inf_res$model_parameters[[i]]$cluster_probs) <- names(data_u[[i]])
  }

  assig_temp = lapply(1:length(data_u), function(i){

    return(
      data.frame(
        id = rownames(data_u[[i]]),
        cluster = inf_res$model_parameters[[i]]$cluster_assignments
      )
    )
    })
  assig_temp = Reduce(assig_temp, f = rbind)

  table = NULL
  if (inherits(data_raw, what = "evopipe_qc"))
  {
    table = data_raw$cnaqc$mutations %>% 
      mutate(id = paste(chr, from, to, sep = ":"))
  } else
  {
    if (inherits(data_raw, what = "cnaqc"))
      table = data_raw$mutations %>%
        mutate(id = paste(chr, from, to, sep = ":"))
    else
      table = data_raw  %>%
        mutate(id = paste(chr, from, to, sep = ":"))
  }

  if (!("is_driver" %in% colnames(table)))
    table$is_driver <-  FALSE
  if (!("driver_label" %in% colnames(table)))
    table$driver_label <-  ""

  if("cluster" %in% colnames(table)){
    cli::cli_alert_warning("A coloumn named cluster already exists, overwriting it!")
    table$cluster <-  NULL
  }

  inf_res$data <-  dplyr::left_join(table, assig_temp, by = "id", copy = T) %>% as.data.frame()

  inf_res$data <- inf_res$data %>% dplyr::mutate(VAF = NV / DP)

  inf_res$description <- description

  inf_res$data$driver_posteriori_annot <-  FALSE
  
  inf_res$used_mutations <- used_muts

  class(inf_res) <- "dbpmmh"

  return(inf_res)

}

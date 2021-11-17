# Error checking inputs
check_inputh = function(
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
                       alpha_prior_concentration,
                       alpha_prior_rate)
{

  stopifnot(all(sapply(subclonal_clusters, function(k) k >= 0)))

  stopifnot(
    tail %in% c(1, 0) | is.null(tail)
  )

  if(lr > 0.05){
    cli::cli_alert_warning("You have selected a relatively high learning rate, consider that such a choice can cause instabilities especially in models with tail.")
  }

  stopifnot(maxIter > 0)
  stopifnot(!is.na(epsilon))
  stopifnot(epsilon > 0)

  stopifnot(
    alpha_prior_concentration > 0 & alpha_prior_rate > 0
  )

  stopifnot(all(sapply(prior_lims_clonal, function(k) k >= 0)))

  stopifnot(all(sapply(prior_lims_k, function(k) k >= 0)))


  stopifnot(model.selection %in% c('ICL', 'BIC', 'AIC', 'likelihood'))


}

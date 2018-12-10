
bootstrapped_statistics = function(fit, bootstrap.fits, alpha = 0.025)
{
  n = length(bootstrap.fits)
  
  pio::pioTit(paste("Bootstrap observations n =", n))
  
  res = NULL
  
  # Overall statistic: "model frequency"
  models.tab = lapply(seq(bootstrap.fits),
               function(w)
                 data.frame(
                   `resample` = w,
                   `tail` = bootstrap.fits[[w]]$fit.tail,
                   `K` =  bootstrap.fits[[w]]$Kbeta,
                   `Model` = paste0(
                     'K = ',
                     bootstrap.fits[[w]]$Kbeta,
                     ifelse(bootstrap.fits[[w]]$fit.tail,
                            " with tail",
                            " without tail")
                   ),
                   stringsAsFactors = FALSE
                 ))
  
  models.tab = Reduce(rbind, models.tab)
  
  model.frequency = tibble::as_tibble(
    sort(table(models.tab$Model)/n, decreasing = TRUE)
  )
  colnames(model.frequency) = c("Model", "Frequency")
  
  this.model = paste0(
    'K = ',
    fit$Kbeta,
    ifelse(fit$fit.tail,
           " with tail",
           " without tail")
  )
  
  model.frequency$fit.model = FALSE
  model.frequency$fit.model[model.frequency$Model == this.model] = TRUE
  
  # Parameter statistics : alpha-level empirical CIs for the bootstrap distribution
  bootstrap.values = lapply(seq(bootstrap.fits),
                            function(w)
                              cbind(bootstrap.fits[[w]]$Clusters, `resample` = w))
  
  bootstrap.values = tibble::as_tibble(Reduce(rbind, bootstrap.values)) %>%
    rename(statistics = type) 
  
  # Stats table
  fit$Clusters = fit$Clusters %>% rename(statistics = type)
  
  stats = bootstrap.values %>%
    group_by(cluster, statistics) %>%
    summarise(
      min = min(fit.value),
      lower_quantile = quantile(fit.value, alpha),
      higher_quantile = quantile(fit.value, 1-alpha),
      max = max(fit.value)
    ) %>%
    left_join(fit$Clusters, by = c('cluster', 'statistics'))
  
  
  pio::pioTit("Bootstrapped model frequency")
  pio::pioDisp(model.frequency)
  
  pio::pioTit("CI (empirical quantiles from bootstrap replicates)")
  
  pio::pioStr("\nMixing proportions", "\n")
  print(stats %>%
          filter(statistics == 'Mixing proportion'))
  
  pio::pioStr("\nTail shape/ scale", "\n")
  print(stats %>%
          filter(cluster == 'Tail' & statistics %in% c('Shape', 'Scale')))
  
  pio::pioStr("\nBeta peaks", "\n")
  print(stats %>%
          filter(statistics %in% c('Mean', 'Variance') & cluster != 'Tail'))
  
  ret = list(bootstrap.values = bootstrap.values,
             models.tab = models.tab,
             stats = stats
             )
  
  invisible(ret)
}

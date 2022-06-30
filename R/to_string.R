#' Return a tabular representation of a model.
#' 
#' @description 
#' 
#' This functionr returns a table with all the parameters fits (one per column), the scores
#' of the model, and the SSE of the fit versus data.
#'
#' @param x A MOBSTER fit.
#'
#' @return
#' @export
#'
#' @examples 
#' data('fit_example', package = 'mobster')
#' 
#' to_string(fit_example$best)
to_string = function(x)
{
  is_mobster_fit(x)
  
  vcz = function(w){
    lblt = w
    w = x$Clusters %>% dplyr::filter(type == w)
    if(nrow(w) == 0) return(NA)
    
    pio:::nmfy(
      paste0(lblt, '.', w  %>% dplyr::pull(cluster)),
      w %>% dplyr::pull(fit.value)
    )
  }
  
  values_f = data.frame(
    tail = x$fit.tail,
    K_beta = x$Kbeta,
    stringsAsFactors = FALSE
  )
  
  values_p = data.frame(
    N = x$N,
    N = x$N.k %>% data.frame %>% t,
    pi = x$pi %>% data.frame %>% t,
    stringsAsFactors = FALSE
  )
  
  values_m = vcz('Mean') %>% data.frame %>% t
  values_v = vcz('Variance') %>% data.frame %>% t
  
  values_sc = data.frame(Scale.Tail = NA)
  if(x$fit.tail) values_sc = vcz('Scale') %>% data.frame %>% t
  
  values_sh = data.frame(Shape.Tail = NA)
  if(x$fit.tail) values_sh = vcz('Shape') %>% data.frame %>% t
  
  # Reasonable clonal cluster
  rcc = sapply(x$z_nk %>% colnames,  mobster:::is_reasonable_clonal_cluster, x = x)
  rcc = pio:::nmfy(x$pi %>% names, rcc)
  rcc = data.frame(
    rcc = rcc %>% data.frame %>% t,
    stringsAsFactors = FALSE
  )
  
  values = cbind(
    values_f,
    values_p,
    values_m,
    values_v,
    values_sc,
    values_sh, 
    rcc
  ) %>%
    as_tibble()
  
  # Re-format and sort columns
  colnames(values) = gsub('\\.', '_',   colnames(values))
  values = values %>% dplyr::select(order(colnames(values)) %>% noquote %>% as.vector)
  
  # Add scores
  values = dplyr::bind_cols(values, x$scores)
  
  # SSE profile (data vs fit)
  sse = mobster:::.compute_fit_sqerr(x, binning = 1e-2)
  
  total_sse = sum(sse$y)
  ranged_sse = sapply(1:10, function(w){
    l = (w - 1) * 0.1
    r = w * 0.1
    v = sum(sse %>% dplyr::filter(x < r, x >= l) %>% dplyr::pull(y))
    names(v) = paste0('sse_', l, '_', r)
    v
  })
  
  ssedf = data.frame(
    sse_total = total_sse,
    t(data.frame(ranged_sse))
  )
  
  values = dplyr::bind_cols(values, ssedf)
  
  # Per cluster error measure
  # densities = per_cluster_err(x)
  #   
  # densities = densities %>%
  #   group_by(cluster) %>%
  #   summarise(median(err))

  return(values)
}

per_cluster_err = function(x)
{
  # Binning o the VAF spectrum
  binning = 1e-3
  x_domain = seq(1e-3, 1-binning, by = binning)
  
  # Hard assignments from the model
  domain = mobster::Clusters_denovo(x, data.frame(VAF = x_domain)) %>%
    dplyr::select(cluster)
  
  # Density of each mixture component in the domain
  densities = mobster:::template_density(
    x,
    x.axis = x_domain,
    binwidth = binning,
    reduce = TRUE) %>%
    tidyr::spread(cluster, y) %>%
    tibble::as_tibble()
  
  densities$cluster = domain$cluster
  
  # Empirical density
  empirical = hist(x$data$VAF, breaks = x_domain, plot = FALSE)$density
  empirical[length(empirical) + 1] = empirical[length(empirical)]
  empirical = empirical * binning # adjust for binwidth
  
  empirical = data.frame(x = x_domain, e = empirical) %>% as_tibble()
  
  densities = densities %>%
    dplyr::left_join(empirical, by = 'x') %>%
    dplyr::filter(x >= min(!!x$data$VAF))
  
  densities$M = NA
  for(j in 1:nrow(densities)) densities$M[j] = unlist(densities[j, densities$cluster[j]])
  
  densities$err = abs(densities$M - densities$e)
  
  densities 
}

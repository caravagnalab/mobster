#' Assign driver mutations in a general karyotype
#'
#' @param x A mobster hierarchical fit object
#' @param rho overdispersion parameter to use in the beta binomial
#'
#' @return
#' @export
#'
#' @examples
#' data()
#'
assign_drivers <- function(x, rho = 0.01){

  drivers <- x$data %>%  filter(is_driver, is.na(cluster), !is.na(karyotype), nchar(karyotype) > 1)  %>%
    tidyr::separate("karyotype", c("Major", "minor"), sep = ":") %>%
    mutate(Tot = Major %>%  as.numeric() + minor %>%  as.numeric())

  if(nrow(drivers) == 0){
    cli::cli_alert_info("No drivers to annotate, returning input object!")
    return(x)
  }

  purity <- x$run_parameters$purity

  p_clonal <- calculate_prob_clonal(drivers,purity, rho)

  if(x$run_parameters$tail){
    p_tail <- calculate_prob_pareto(drivers,lower =  x$data$VAF %>% min(), shape = min(get_pareto(x)$shape),
                                    scale = min(get_pareto(x)$scale),trunc = x$run_parameters$truncated_pareto,
                                    purity = purity)
  } else {
    p_tail = (p_clonal * 0)
  }

  if(x$run_parameters$K > 0){
    if(x$run_parameters$K == 1){
      p_subclonal <- calculate_prob_subclonal(x1 = drivers, ccfs = x$model_parameters[[1]]$ccf_subclones,
                                              purity = purity, rho =  rho) %>%  matrix(nrow = 1)
    } else {
      p_subclonal <- calculate_prob_subclonal(x1 = drivers, ccfs = x$model_parameters[[1]]$ccf_subclones,
                                              purity = purity, rho =  rho)
    }

  }

  names <-  c("C*", "Tail", paste0("S", seq_along(x$model_parameters[[1]]$ccf_subclones) ))
  probs <- list(p_clonal, p_tail)
  if(x$run_parameters$K > 0){
    for(i in 1:length(x$model_parameters[[1]]$ccf_subclones) ) {

      probs[[2 + i]] <- p_subclonal[i,]
    }
  } else {
    probs[[3]] <- (p_clonal * 0)
  }

  probs <-  data.frame(probs)

  colnames(probs) <-  names

  new_clust <- names[apply(probs,1,which.max)]

  x$data[x$data$is_driver & is.na(x$data$cluster) & !is.na(x$data$karyotype) & nchar(x$data$karyotype) > 1,] <-
    x$data %>%  filter(is_driver, is.na(cluster), !is.na(karyotype), nchar(karyotype) > 1) %>% mutate(cluster = new_clust, driver_posteriori_annot = T)

  return(x)

}


calculate_prob_pareto <- function(x1, lower, shape, scale, trunc, purity){


  p_tail <- apply(x1, 1, FUN = calculate_prob_pareto_aux, lower=lower, shape = shape,
                  scale = scale, trunc = trunc, purity = purity)


  return(p_tail)

}

calculate_prob_pareto_aux <- function(df, lower, shape, scale, trunc, purity){


  p_tail <- integrate(pareto_beta_bin_mix, lower = lower, upper = ifelse(trunc, get_expected_peak(( 1 - 1e-2) , purity,  df["Tot"] %>%  as.numeric()), 0.999),
                      NV = df["NV"]%>%  as.numeric(), DP = df["DP"]%>%  as.numeric(),
                      shape = shape * df["Tot"] %>%  as.numeric(), scale = scale)$value

}

pareto_beta_bin_mix <- function(p, NV, DP, shape, scale){

  dbinom(x = NV, size = DP, prob = p) * dpareto(p,shape = shape, scale = scale)

}

calculate_prob_subclonal <- function(x1, ccfs,purity, rho = 0.01){

  p_subclone <- apply(x1, 1, FUN = calculate_prob_subclonal_aux, rho = rho, ccfs = ccfs, purity = purity)


  return(p_subclone)

}

calculate_prob_subclonal_aux <-  function(df, ccfs,purity, rho) {

  res <-  VGAM::dbetabinom(x = df["NV"]%>%  as.numeric(),
                                                     size = df["DP"]%>%  as.numeric(),
                                                     prob = get_expected_peak((1 - 1e-2), purity, df["Tot"] %>%  as.numeric(), ccfs),
                                                     rho = rho)
  return(res)
}

calculate_prob_clonal <- function(x1, purity, rho = 0.01){

  p_clone <- apply(x1, 1, FUN = calculate_prob_clonal_aux, rho = rho, purity = purity)


  return(p_clone)

}


calculate_prob_clonal_aux <- function(df, rho, purity){

  res <-  sapply(1:(df["Major"]%>%  as.numeric()), function(i)
    VGAM::dbetabinom(x = df["NV"] %>%  as.numeric(), size = df["DP"]%>%  as.numeric(),
                     prob = get_expected_peak((i -1e-2) , purity,  df["Tot"] %>%  as.numeric()),
                     rho = rho))
  return(res %>% max)
}



get_expected_peak <- function(c,purity, tot,ccf = 1) {

  (c * purity * ccf) / (2 * (1-purity) + purity * tot)

}

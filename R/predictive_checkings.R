

predictive_checks <- function(obj, samples = 100, type = c("prior", "posterior")){
  
  if(type == "prior")
    return(prior_predictive_checks(obj, samples))
  if(type == "posterior")
    return(posterior_predictive_checks(obj, samples))
  
}



prior_predictive_checks <- function(obj, samples = 100){
  
  samples <- easypar::run(prior_predictive_checks_aux,lapply(1:samples, list),parallel = F,export = "obj")
  return(samples)
  
  
}


prior_predictive_checks_aux <- function(i){
  
  res <- vector(length = length(names(obj$model_parameters)), mode = "list")
  names(res) <- names(obj$model_parameters)
  
  if(!obj$run_parameters$multi_tail){
    alpha_prior = rgamma(shape = 3,rate = 3, n = 1)
  }else {
    alpha_prior = rgamma(shape = 10,rate = 10, n = 1)
  }
    
  
  if(obj$run_parameters$K > 0){
    subclonal_ccf = runif(0.05, 0.95, n = obj$run_parameters$K)
  }
  
  for(k in names(obj$model_parameters)){
    
    theo_peaks = (theo_clonal_num(k) * obj$run_parameters$purity - 1e-9) / (2 * (1 - obj$run_parameters$purity) + theo_clonal_tot(k) * obj$run_parameters$purity)
    
    
    prior_overdispersion = runif(obj$run_parameters$prior_lims_clonal[1], obj$run_parameters$prior_lims_clonal[2], n = 1)
    prec_overdispersion = rgamma(shape = 3,rate = 1, n = 1)

    
    weights = gtools::rdirichlet(n = 1, alpha = 1/ rep(1,(obj$run_parameters$K + theo_clonal_num(k, range = F))))
      
    ## Clonal component
    
    bm_1 = theo_peaks * obj$run_parameters$number_of_trials_clonal_mean
    
    bm_2 = obj$run_parameters$number_of_trials_clonal_mean - bm_1
    
    betas_clone_mean = rbeta(shape1 = bm_1, shape2 = bm_2, n = length(bm_1))
    
    betas_clone_n_samples = rlnorm(n = length(bm_1), log(prior_overdispersion), 1/prec_overdispersion)
    
    
    if(obj$run_parameters$K > 0){
      
      adj_ccf = (subclonal_ccf * obj$run_parameters$purity) / (2 * (1 - obj$run_parameters$purity) + theo_clonal_tot(k) * obj$run_parameters$purity)
      k_means = runif(n = obj$run_parameters$K, ((subclonal_ccf + 0.01) * obj$run_parameters$purity) / (2 * (1 - obj$run_parameters$purity) + theo_clonal_tot(k) * obj$run_parameters$purity),
                      ((subclonal_ccf + 0.01) * obj$run_parameters$purity) / (2 * (1 - obj$run_parameters$purity) + theo_clonal_tot(k) * obj$run_parameters$purity))
      
      if (obj$run_parameters$subclonal_prior == "Moyal"){
        scale_subclonal = rgamma(n = obj$run_parameters$K, 10,0.1)
        subclone_mean = rtruncmoyal(n = obj$run_parameters$K, 
                                    loc = k_means - 1/scale_subclonal * (-digamma(1) + log(2)), scale = 1/scale_subclonal, 
                                    lower = obj$model_parameters[[k]]$tail_scale - 1e-5,upper = min(theo_peaks))
      } else {
        
        num_trials_subclonal = runif( obj$run_parameters$K, obj$run_parameters$prior_lims_k[1], obj$run_parameters$prior_lims_k[2])
        subclone_mean = rbeta(n = obj$run_parameters$K,shape1 = k_means * num_trials_subclonal,shape2 =  (1 - k_means) * num_trials_subclonal)
      }
    }
    
    if (obj$run_parameters$tail > 0){
      
      ## TAIL distributions
      
      tail_probs = gtools::rdirichlet(n = 1,alpha = c(1 , 1 + obj$run_parameters$K))
      alpha_precision = rgamma(1, obj$run_parameters$alpha_precision_concentration, obj$run_parameters$alpha_precision_rate)
      alpha = rlnorm(meanlog = log(alpha_prior),sdlog =  1 / alpha_precision, n = obj$run_parameters$tail)
      
      if( obj$run_parameters$truncated_pareto){
        if( obj$run_parameters$K > 0 & obj$run_parameters$multi_tail){
          multitails_weights = gtools::rdirichlet(n = 1,alpha = rep(1 , 1 + obj$run_parameters$K)) 
          tcm = theo_peaks - max(adj_ccf)
          tcm[tcm < obj$model_parameters[[k]]$tail_scale] = obj$model_parameters[[k]]$tail_scale
          tcm = c(tcm, adccf)
          U = tcm
        }else{
          multitails_weights = NA
          U = min(theo_peaks)
        }
      }else{
        multitails_weights = NA
        U = 1
      }
    }
      
      ## Clonal likelihood
      
      data_l = length(obj$model_parameters[[k]]$cluster_assignments)
      DP = obj$data %>% filter(karyotype == k, !is.na(cluster)) %>% pull(DP)
      beta = lapply(1:length(betas_clone_mean), FUN =  function(i) VGAM::rbetabinom.ab(shape1 = betas_clone_mean[i] * betas_clone_n_samples[i],
                                                                                       shape2 =  (1 - betas_clone_mean[i]) * betas_clone_n_samples[i], 
                                 size = DP, n = data_l)) %>% do.call(cbind,.)
      ## Subclonal likelihood
      
      if ( obj$run_parameters$K > 0){
        has_subclones = TRUE
          subclonal_lk = lapply(1:length(subclone_mean), 
                                FUN =  function(i) rbinom(prob = subclone_mean[i],
                                                          n = data_l, size = DP)) %>% do.call(cbind,.)
      }else{
        has_subclones = FALSE
      }
        
      
      if (obj$run_parameters$tail == 1){
        has_tail = TRUE
        multi_penalty = 0
        if (obj$run_parameters$K > 0 & obj$run_parameters$truncated_pareto & obj$run_parameters$multi_tail) {
          clonal_prop = 1 - max(subclonal_ccf)
          sub_ccf = subclonal_ccf
          theo_weights = c(clonal_prop, sub_ccf)
          theo_weights = theo_weights / sum(theo_weights)
          multi_penalty = log(length(NV)) * sqrt((theo_weights - multitails_weights)**2)
        }
          
        
        pareto = VGAM::rtruncpareto( n = data_l, shape = alpha, lower = obj$model_parameter[[k]]$tail_scale, upper =  U)
        pareto_bin = rbinom(prob = pareto, 
                            size = DP, n = data_l)
        
      } 
      if (has_tail & has_subclones){
        
        all_samples = cbind(beta, subclonal_lk)
        all_samples = cbind(all_samples, pareto_bin)
        probs = c(weights * tail_probs[2], tail_probs[1])
        idx = sample.int(size = nrow(all_samples),n = ncol(all_samples), prob = probs, replace = T)
        final_sample = lapply(1:nrow(all_samples), function(i) all_samples[i,idx[i]]) %>%
          do.call(c,.) %>% unname()
      
      }
      
      if (has_tail &  !has_subclones) {
        all_samples = cbind(beta, pareto_bin)
        probs = c(weights * tail_probs[1], tail_probs[2])
        idx = sample.int(size = nrow(all_samples),n = ncol(all_samples), prob = probs, replace = T)
        final_sample = lapply(1:nrow(all_samples), function(i) all_samples[i,idx[i]]) %>%
          do.call(c,.) %>% unname()
      }
        
      
      if ( !has_tail & has_subclones){
        
        all_samples = cbind(beta, subclonal_lk)
        probs = weights 
        idx = sample.int(size = nrow(all_samples),n = ncol(all_samples), prob = probs, replace = T)
        final_sample = lapply(1:nrow(all_samples), function(i) all_samples[i,idx[i]]) %>%
          do.call(c,.) %>% unname()
      }

      
      if ( !has_tail &  !has_subclones){
        
        all_samples = beta
        probs = weights
        idx = sample.int(size = nrow(all_samples),n = ncol(all_samples), prob = probs, replace = T)
        final_sample = lapply(1:nrow(all_samples), function(i) all_samples[i,idx[i]]) %>%
          do.call(c,.) %>% unname()

      }
      res[[k]] <-final_sample
  }
  return(res)
}
  
  


theo_clonal_num <- function(kr, range = T){
  kr = stringr::str_split(kr, pattern = ":")[[1]]
  kr = as.integer(kr)
  if(range)
    return(1:max(kr))
  return(max(kr))
}

theo_clonal_tot <- function(kr) {
  kr = stringr::str_split(kr, pattern = ":")[[1]]
  kr = as.integer(kr)
  return(sum(kr))
}

scale_pareto <- function(VAF, min_vaf_scale_tail){
  NBINS = 100
  hist_v = hist(VAF, breaks = seq.default(0,1,length.out = NBINS),plot = F)
  vals = hist_v$breaks
  
  idx = which(hist_v$counts[1:(length(hist_v$counts)-1)] > hist_v$counts[2:length(hist_v$counts)])
  
  best_scale = vals[idx[1]]
  
  if (best_scale > max_vaf)
    return(min(VAF) - 1e-10)
  else
    return(best_scale - 1e-10)
}


log_sum_exp <- function(M){
  c = apply(M, 1, max)
  return (c + log(apply(exp(args - M), 1, sum)))
}



posterior_predictive_checks <- function(obj, samples = 100){
  
  samples <- easypar::run(posterior_predictive_checks_aux,lapply(1:samples, list),parallel = F,export = "obj")
  return(samples)
  
  
}


posterior_predictive_checks_aux <- function(i){
  
  res <- vector(length = length(names(obj$model_parameters)), mode = "list")
  names(res) <- names(obj$model_parameters)
  
 
  
  for(k in names(obj$model_parameters)){
    
    data_l = length(obj$model_parameters[[k]]$cluster_assignments)
    
    
    alpha <- obj$model_parameters[[k]]$tail_shape
    
    
    if(obj$run_parameters$K > 0){
      subclonal_ccf = obj$model_parameters[[k]]$ccf_subclones
    }
    
    theo_peaks = (theo_clonal_num(k) * obj$run_parameters$purity - 1e-9) / (2 * (1 - obj$run_parameters$purity) + theo_clonal_tot(k) * obj$run_parameters$purity)
    
    
    prior_overdispersion = obj$model_parameters[[k]]$beta_concentration1 + obj$model_parameters[[k]]$beta_concentration2
    prec_overdispersion = 1/obj$model_parameters[[k]]$dispersion_noise
    
    
    weights = obj$model_parameters[[k]]$mixture_probs
    
    ## Clonal component
    
    bm_1 =  obj$model_parameters[[k]]$beta_concentration1
    
    bm_2 = obj$model_parameters[[k]]$beta_concentration2
    
    betas_clone_mean_p = bm_1/(bm_1 + bm_2)
    
    betas_clone_n_samples = rlnorm(n = length(bm_1), log(prior_overdispersion), 1/prec_overdispersion)
    
    betas_clone_mean <- lapply(1:length(betas_clone_mean_p), function(x) rbeta(shape1 = betas_clone_mean_p[x] * betas_clone_n_samples[x],
                                                                               shape2 =  (1 - betas_clone_mean_p[x]) * betas_clone_n_samples[x], 
                                                                               n = data_l))
    
    if(obj$run_parameters$K > 0){
      
      adj_ccf = (subclonal_ccf * obj$run_parameters$purity) / (2 * (1 - obj$run_parameters$purity) + theo_clonal_tot(k) * obj$run_parameters$purity)
      k_means = runif(n = obj$run_parameters$K, ((subclonal_ccf + 0.01) * obj$run_parameters$purity) / (2 * (1 - obj$run_parameters$purity) + theo_clonal_tot(k) * obj$run_parameters$purity),
                      ((subclonal_ccf + 0.01) * obj$run_parameters$purity) / (2 * (1 - obj$run_parameters$purity) + theo_clonal_tot(k) * obj$run_parameters$purity))
      
      if (obj$run_parameters$subclonal_prior == "Moyal"){
        scale_subclonal = obj$model_parameters[[k]]$scale_subclonal
        subclone_mean = rtruncmoyal(n = obj$run_parameters$K, 
                                    loc = k_means - 1/scale_subclonal * (-digamma(1) + log(2)), scale = 1/scale_subclonal, 
                                    lower = obj$model_parameters[[k]]$tail_scale - 1e-5,upper = min(theo_peaks))
      } else {
        
        num_trials_subclonal = obj$model_parameters[[k]]$n_trials_subclonal
        subclone_mean = rbeta(n = obj$run_parameters$K, k_means * num_trials_subclonal, (1 - k_means) * num_trials_subclonal)
      }
    }
    
    if (obj$run_parameters$tail > 0){
      
      ## TAIL distributions
      
      alpha_precision = 1/obj$model_parameters[[k]]$tail_noise
      alpha = rlnorm(meanlog = log(alpha),sdlog =  1 / alpha_precision, n = obj$run_parameters$tail)
      
      if( obj$run_parameters$truncated_pareto){
        if( obj$run_parameters$K > 0 & obj$run_parameters$multi_tail){
          multitails_weights = obj$model_parameters[[k]]$multi_tail_weights
          tcm = theo_peaks - max(adj_ccf)
          tcm[tcm < obj$model_parameters[[k]]$tail_scale] = obj$model_parameters[[k]]$tail_scale
          tcm = c(tcm, adccf)
          U = tcm
        }else{
          multitails_weights = NA
          U = min(theo_peaks)
        }
      }else{
        multitails_weights = NA
        U = 1
      }
    }
    
    ## Clonal likelihood
    
    DP = obj$data %>% filter(karyotype == k, !is.na(cluster)) %>% pull(DP)
    beta = lapply(1:length(betas_clone_mean), FUN =  function(i) rbinom(prob  = betas_clone_mean[[i]], 
                                                                                     size = DP, n = data_l)) %>% do.call(cbind,.)
    ## Subclonal likelihood
    
    if ( obj$run_parameters$K > 0){
      has_subclones = TRUE
      subclonal_lk = lapply(1:length(subclone_mean), 
                            FUN =  function(i) rbinom(prob = subclone_mean[i],
                                                      n = data_l, size = DP)) %>% do.call(cbind,.)
    }else{
      has_subclones = FALSE
    }
    
    
    if (obj$run_parameters$tail == 1){
      has_tail = TRUE
      multi_penalty = 0
      if (obj$run_parameters$K > 0 & obj$run_parameters$truncated_pareto & obj$run_parameters$multi_tail) {
        clonal_prop = 1 - max(subclonal_ccf)
        sub_ccf = subclonal_ccf
        theo_weights = c(clonal_prop, sub_ccf)
        theo_weights = theo_weights / sum(theo_weights)
        multi_penalty = log(length(NV)) * sqrt((theo_weights - multitails_weights)**2)
      }
      
      
      pareto = VGAM::rtruncpareto( n = data_l, shape = alpha, lower = obj$model_parameter[[k]]$tail_scale, upper =  U)
      pareto_bin = rbinom(prob = pareto, 
                          size = DP, n = data_l)
      
    } 
    if (has_tail & has_subclones){
      
      all_samples = cbind(pareto_bin, subclonal_lk)
      all_samples = cbind(all_samples, beta)
      probs = weights
      idx = sample.int(size = nrow(all_samples),n = ncol(all_samples), prob = probs, replace = T)
      final_sample = lapply(1:nrow(all_samples), function(i) all_samples[i,idx[i]]) %>%
        do.call(c,.) %>% unname()
      
    }
    
    if (has_tail &  !has_subclones) {
      all_samples = cbind(pareto_bin, beta)
      probs = weights 
      idx = sample.int(size = nrow(all_samples),n = ncol(all_samples), prob = probs, replace = T)
      final_sample = lapply(1:nrow(all_samples), function(i) all_samples[i,idx[i]]) %>%
        do.call(c,.) %>% unname()
    }
    
    
    if ( !has_tail & has_subclones){
      
      all_samples = cbind(subclonal_lk, beta)
      probs = weights 
      idx = sample.int(size = nrow(all_samples),n = ncol(all_samples), prob = probs, replace = T)
      final_sample = lapply(1:nrow(all_samples), function(i) all_samples[i,idx[i]]) %>%
        do.call(c,.) %>% unname()
    }
    
    
    if ( !has_tail &  !has_subclones){
      
      all_samples = beta
      probs = weights
      idx = sample.int(size = nrow(all_samples),n = ncol(all_samples), prob = probs, replace = T)
      final_sample = lapply(1:nrow(all_samples), function(i) all_samples[i,idx[i]]) %>%
        do.call(c,.) %>% unname()
      
    }
    res[[k]] <-final_sample
  }
  return(res)
}


calculate_distance_ecdf <- function(obj, samples, used_mutations, karyo_collapse = c("none", "sum", "mean", "max")){
  
  data_real <- obj$data %>% mutate(mutation_id = paste(chr,from,to, sep = ":")) %>% filter(mutation_id %in% used_mutations) %>% select(NV,DP, karyotype)
  data_real <- split(data_real, data_real$karyotype, drop = T)
  data_real <- lapply(data_real, function(x) x[,1])
  
  dist_vec = vector(length = 2L)
  names(dist_vec) <- names(obj$model_parameters)
  
  for(k in names(obj$model_parameters)) {
    x <- seq.default(0,max(data_real[[k]]) + 10,by = 1)
    
    predictive_samples <- lapply(samples, function(x) x[[k]]) %>% do.call("cbind",.)
    
    avg_post <- apply(predictive_samples,1, mean)
    dist_k <- ks.test(data_real[[k]], avg_post)$statistic
    
    dist_vec[k] <- dist_k
    
  }
  
  if(karyo_collapse == "none"){
    return(dist_vec)
  }else if (karyo_collapse == "mean") {
    return(mean(dist_vec))
  }else {
    return(do.call(karyo_collapse,as.list(dist_vec)))
  }
  
}

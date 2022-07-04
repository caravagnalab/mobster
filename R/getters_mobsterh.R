###################################################
# Auxiliary functions to extract information from #
# a fit with MOBSTERh                             #
###################################################

# Parameters getters for the tail
get_pareto = function(x) {
  used = x$model_parameters %>% names

  if(!has_tail(x)) return(NULL)

  one_Beta = c("1:0", "1:1")
  two_Beta = c("2:0", "2:1", "2:2")

  df = NULL
  for (k in used){
    beta_mean_k <- x$model_parameters[[k]]$beta_concentration1 / (x$model_parameters[[k]]$beta_concentration2 + x$model_parameters[[k]]$beta_concentration1)
    min_beta_clonal_mean <- ifelse(test = k %in% one_Beta,
                                   beta_mean_k,
                                   min(beta_mean_k))
    df = df %>%
    bind_rows(
      data.frame(
        karyotype = k,
        scale = x$model_parameters[[k]]$tail_scale,
        shape = x$model_parameters[[k]]$tail_shape,
        mixing = x$model_parameters[[k]]$mixture_probs[1],
        shape_noise = x$model_parameters[[k]]$tail_noise,
        location = ifelse(is_truncated(x),
                          min_beta_clonal_mean,
                          1 - 1e-5),
        cluster = "Tail"
      )
    )
  }

  df
}

# Parameters getters for the Betas
get_beta = function(x) {
  used = x$model_parameters %>% names

  one_Beta = c("1:0", "1:1")
  two_Beta = c("2:0", "2:1", "2:2")

  df = NULL
  for (k in used)
  {
    betas = x$model_parameters[[k]] %>% names
    betas = startsWith(betas, 'beta')

    tibble_beta = x$model_parameters[[k]][betas] %>% as_tibble()
    colnames(tibble_beta) = c("a", "b")

    tibble_beta$karyotype = k
    cluster <-  vector(length = nrow(tibble_beta))
    if(k %in% one_Beta){
      cluster[1] = "C1"
      if(mobster:::has_tail(x))
        tibble_beta$mixing = x$model_parameters[[k]]$mixture_probs[2]
      else
        tibble_beta$mixing = x$model_parameters[[k]]$mixture_probs[1]
    } else {
      cluster[1:2] = c("C1", "C2")
      if(mobster:::has_tail(x))
        tibble_beta$mixing = x$model_parameters[[k]]$mixture_probs[2:3]
      else
        tibble_beta$mixing = x$model_parameters[[k]]$mixture_probs[1:2]
    }

    tibble_beta$cluster <-  cluster


    df = df %>% bind_rows(tibble_beta)
  }

  df
}

get_ccf_subclones <- function(x) {
  
  
  if(!mobster:::has_subclones(x)) return(NULL)
  
  used = x$model_parameters %>% names
  

  k = used[1]

  tibble_ccf <- x$model_parameters[[k]]$ccf_subclones %>% as_tibble()
  colnames(tibble_ccf) <- "CCF"
  tibble_ccf$cluster <-  paste0("S", 1:nrow(tibble_ccf))
  
  return(tibble_ccf)
}

get_moyal_sub <- function(x) {
  if(!mobster:::has_subclones(x)) return(NULL)
  used = x$model_parameters %>% names
  
  df = NULL
  
  for(k in used) {
    tibble_moyal <- x$model_parameters[[k]]$loc_subclones %>% as_tibble()
    colnames(tibble_moyal) <- "location"
    tibble_moyal$scale <- x$model_parameters[[k]]$scale_subclonal
    tibble_moyal$karyotype <- k
    tibble_moyal$cluster <-  paste0("S", 1:nrow(tibble_moyal))
    L <- length(x$model_parameters[[k]]$mixture_probs)
    K <- nrow(tibble_moyal)
    tibble_moyal$mixing <- x$model_parameters[[k]]$mixture_probs[(L-K+1):L]
    
    df = df %>% bind_rows(tibble_moyal)
    
  }
  
  return(df)
  
}


get_beta_sub <- function(x) {
  if(!mobster:::has_subclones(x)) return(NULL)
  used = x$model_parameters %>% names
  
  df = NULL
  
  for(k in used) {
    
    tibble_beta <- (x$model_parameters[[k]]$loc_subclones * 
                      x$model_parameters[[k]]$n_trials_subclonal) %>% as_tibble()
    colnames(tibble_beta) <- "a"
    tibble_beta$b <- x$model_parameters[[k]]$n_trials_subclonal - tibble_beta$a 
    
    tibble_beta$karyotype = k
    tibble_beta$cluster <-  paste0("S", 1:nrow(tibble_beta))
    
    L <- length(x$model_parameters[[k]]$mixture_probs)
    K <- nrow(tibble_beta)
    tibble_beta$mixing <- x$model_parameters[[k]]$mixture_probs[(L-K+1):L]
    
    df = df %>% bind_rows(tibble_beta)
    
  }
  
  return(df)
  
}


get_mixture_weights <- function(x){

  return(clonality_interpreter(x))

}

get_assignment_probs <- function(x, cutoff = 0){

  used = x$model_parameters %>% names
  clusts <-  clonality_interpreter(x)
  df = NULL
  for(k in used){
    k_probs <- x$model_parameters[[k]]$cluster_probs %>% t
    k_probs <-  k_probs %>%  as.data.frame()
    k_probs$cluste <-  apply(k_probs, 1, function(kk) {
      if(max(kk) < cutoff) return(NA)
      else return(which.max(kk))
    })
    k_probs <-  k_probs  %>%  arrange(-cluste, across()) %>% select(-cluste) %>% as.matrix()
    colnames(k_probs) <- clusts %>% filter(karyotype == !!k) %>% pull(cluster)
    k_probs <- k_probs %>%  reshape2::melt() %>% as_tibble()
    colnames(k_probs) <-  c("Point", "Cluster", "Value")
    k_probs$Karyotype <-  k
    to_filt <-  k_probs %>%  group_by(Point) %>%  summarize(Max = max(Value))
    k_probs <-  dplyr::full_join(to_filt, k_probs, by = "Point")
    k_probs <-  k_probs %>%  mutate(Cluster = if_else(Max > cutoff, as.character(Cluster), "NA"))
    df = df %>%  bind_rows(k_probs)

  }

  return(df)
}



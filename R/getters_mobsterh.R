###################################################
# Auxiliary functions to extract information from #
# a fit with MOBSTERh                             #
###################################################

# Parameters getters for the tail
get_pareto = function(x) {
  used = x$model_parameters %>% names

  df = NULL
  for (k in used)
    df = df %>%
    bind_rows(
      data.frame(
        karyotype = k,
        scale = x$model_parameters[[k]]$tail_scale,
        shape = x$model_parameters[[k]]$tail_shape,
        mixing = x$model_parameters[[k]]$mixture_probs[1],
        shape_noise = x$model_parameters[[k]]$tail_noise,
        location = if(is_truncated(x), min(x$model_parameters[[k]]$), NULL)
        cluster = "Tail"
      )
    )

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

    if(has_tail(x))
      tibble_beta$mixing = x$model_parameters[[k]]$mixture_probs[-1]
    else
      tibble_beta$mixing = x$model_parameters[[k]]$mixture_probs
    tibble_beta$karyotype = k
    cluster <-  vector(length = nrow(tibble_beta))
    if(k %in% one_Beta){
      cluster[1] = "C1"
      if(nrow(tibble_beta) > 1)
        cluster[2:nrow(tibble_beta)] <-  paste0("S", 1:(nrow(tibble_beta)-1))
    } else {
      cluster[1:2] = c("C1", "C2")
      if(nrow(tibble_beta) > 2)
        cluster[3:nrow(tibble_beta)] <-  paste0("S", 1:(nrow(tibble_beta)-2))
    }

    tibble_beta$cluster <-  cluster


    df = df %>% bind_rows(tibble_beta)
  }

  df
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



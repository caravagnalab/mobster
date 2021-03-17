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
        cluster = "Tail"
      )
    )

  df
}

# Parameters getters for the Betas
get_beta = function(x) {
  used = x$model_parameters %>% names

  df = NULL
  for (k in used)
  {
    betas = x$model_parameters[[k]] %>% names
    betas = startsWith(betas, 'beta')

    tibble_beta = x$model_parameters[[k]][betas] %>% as_tibble()
    colnames(tibble_beta) = c("a", "b")

    tibble_beta$mixing = x$model_parameters[[k]]$mixture_probs[-1]
    tibble_beta$karyotype = k
    tibble_beta$cluster = paste0("Beta", 1:nrow(tibble_beta))

    df = df %>% bind_rows(tibble_beta)
  }

  df
}

# Interpret clonality from a MOBSTERh fit
clonality_interpreter = function(x)
{
  tail_params = get_pareto(x) %>%
    mutate(what = "Neutral") %>%
    select(karyotype,
           cluster,
           what)

  one_Beta = c("1:0", "1:1")
  two_Beta = c("2:0", "2:1", "2:2")

  Beta_params = get_beta(x) %>% mutate(what = case_when(
    (karyotype %in% one_Beta) & (cluster == "Beta1") ~ "Clonal",
    (karyotype %in% one_Beta) & (cluster != "Beta1") ~ "Subclone",
    (karyotype %in% two_Beta) &
      (cluster %in% c("Beta1", "Beta2")) ~ "Clonal",
    (karyotype %in% two_Beta) &
      !(cluster %in% c("Beta1", "Beta2")) ~ "Subclone"
  )) %>%
    select(karyotype,
           cluster,
           what)

  return(tail_params %>%
           bind_rows(Beta_params))
}

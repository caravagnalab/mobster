# Interpret clonality from a MOBSTERh fit
clonality_interpreter = function(x)
{



  Beta_params = get_beta(x) %>% mutate(what = dplyr::case_when(
    (substr(cluster,1,1) == "C") ~ "Clone",
    (substr(cluster,1,1) == "S") ~ "Subclone")) %>%
    select(karyotype,
           cluster,
           what, mixing)

  if(has_tail(x)){
    tail_params = get_pareto(x) %>%
      mutate(what = "Neutral") %>%
      select(karyotype,
             cluster,
             what, mixing)

    return(tail_params %>%
             bind_rows(Beta_params))

  } else {
    return(Beta_params)
  }


}

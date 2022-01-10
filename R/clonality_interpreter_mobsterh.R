# Interpret clonality from a MOBSTERh fit
clonality_interpreter = function(x)
{


  Beta_params = get_beta(x) %>% mutate(what = "Clone") %>%
    select(karyotype,
           cluster,
           what, mixing)
  
  if(has_subclones(x)) {
    if(is_moyal(x)) {
      subclonal_params = get_moyal_sub(x) %>% mutate(what = "Subclone") %>%
        select(karyotype,
               cluster,
               what, mixing)
    } else {
      subclonal_params = get_beta_sub(x) %>% mutate(what = "Subclone") %>%
        select(karyotype,
               cluster,
               what, mixing)
    }
    
    Beta_params <- Beta_params %>% bind_rows(subclonal_params)
  }

  if(has_tail(x)){
    tail_params = get_pareto(x) %>%
      mutate(what = "Neutral") %>%
      select(karyotype,
             cluster,
             what, mixing)

    Beta_params <- Beta_params %>%
             bind_rows(tail_params)
  }

  return(Beta_params)

}

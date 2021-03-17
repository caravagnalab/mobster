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

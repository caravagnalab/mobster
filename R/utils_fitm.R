
format_data_mobsterm <- function(x) {

  x_unique <- x %>%
    distinct(sample_id, .keep_all = TRUE)
  
  sample_names <- x_unique$sample_id
  purity <- x_unique$purity
  karyotype <- x_unique$karyotype
  
  x <- x %>%
    dplyr::select(NV, DP, mutation_id, sample_id, purity, karyotype)
  
  # From long to wide format
  x_wide <- x %>%
    pivot_wider(
      names_from = sample_id,
      values_from = c(NV, DP, karyotype, purity),
      names_sep = "_" # separator between variable and sample_id
    )
  
  # Remove rows with NA
  x_wide = na.omit(x_wide)
  
  # Extract NV and DP for python code
  NV <- x_wide %>%
    select(starts_with("NV_")) %>%
    as.matrix()
  
  DP <- x_wide %>%
    select(starts_with("DP_")) %>%
    as.matrix()
  
  mutation_id = x_wide %>%
    select(mutation_id) %>% 
    as.matrix()
  
  # Add a "cluster" column
  x_wide$cluster <-  NULL
  
  
  return(list(
    x_wide = x_wide,
    NV = NV,
    DP = DP,
    purity = purity,
    mutation_id = mutation_id,
    karyotype = karyotype
  ))
  
}


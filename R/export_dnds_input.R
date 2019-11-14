extract_dnds_input = function(x, genes)
{
  message("TODO document this function")
  
  mobster:::is_mobster_fit(x)
  
  if("gene" %in% colnames(x$data))
  {
    pio::pioStr("Found 'gene' column in the data, using that.")
    pio::pioDisp(x %>% filter(gene %in% genes))
  }
  
  pio::pioStr("Missing 'gene' column in the data, annotating that.")
  
  x = x %>% 
    annotate_genes_from_locations() %>%
    filter(gene %in% genes) 
  
  x %>% pio::pioDisp
  
  return(x)
}

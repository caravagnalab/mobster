extract_dnds_input = function(x, genes_list = NULL)
{
  message("TODO document this function")
  
  mobster:::is_mobster_fit(x)
  
  if("gene" %in% colnames(x$data))
  {
    pio::pioStr("Found 'gene' column in the data, using that.", suffix = '\n')
  }
  else
  {
    pio::pioStr("Missing 'gene' column in the data, annotating that.", suffix = '\n')
  
    x = x %>% annotate_genes_from_locations
  }
  
  if(!is.null(genes_list))
    x = x %>% dplyr::filter(gene %in% genes_list) 
  
  pio::pioStr("Annotated data", ifelse(is.null(genes_list), "All genes", "Subset to custom genes list"), suffix = '\n')
  
  x %>% pio::pioDisp()                                                                                                                                                                                       
  
  return(x)
}

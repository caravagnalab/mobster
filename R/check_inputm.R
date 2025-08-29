# Error checking inputs
check_inputm = function(
    x,
    K_list,
    max_iter,
    seed_list, 
    lr,
    par_threshold, 
    loss_threshold)
{
  
  stopifnot(all(sapply(K_list, function(k) k >= 0))) # Check all K values are positive
  
  if (!is.list(K_list)) K_list = as.list(K_list) # Check K values is a list
  if (!is.list(seed_list)) seed_list = as.list(seed_list) # Check seed values is a list
  
  if(lr > 0.05){
    cli::cli_alert_warning("You have selected a relatively high learning rate, consider that such a choice can cause instabilities.")
  }
  
  stopifnot(max_iter > 0)
  stopifnot(par_threshold > 0)
  stopifnot(loss_threshold > 0)
  
  # Here check that karyotype is unique per sample
  n_samples <- length(unique(x$sample_id))
  n_pairs   <- nrow(unique(x[c("sample_id", "karyotype")]))
  
  if(n_samples != n_pairs){
    stop(
      "Some sample_id values are associated with multiple karyotypes. Each sample_id must have exactly one karyotype."
    )
  }
  
}
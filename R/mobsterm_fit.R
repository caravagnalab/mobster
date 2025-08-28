
#' Fit a multivariate model with MOBSTERm
#'
#' @description This function fits a multivariate version of the MOBSTER model implemented in \code{mobster_fit}. 
#' 
#' @param x Input data.frame (or tibble). The input data.frame should have at least 
#' 6 coloumns named as mutation_id, NV (number of variant), DP (depth), sample_id, purity and karyotype.
#' @param K_list A vector with the number of mixture components to use. All values of \code{K_list} must be positive
#' and strictly greater than 0.
#' @param max_iter Maximum number of steps for a variational inference fit. 
#' If parameters and loss convergence is not achieved before these steps, the fit is interrupted.
#' @param seed A vector with the number of seeds to test for the fit.
#' @param lr Learning rate used by the optimizer.
#' @param par_threshold Tolerance for parameter convergence. As ELBO oscillations are common in gradient based VI, we will monitor the convergence of all the parameters in the model,
#' the inference stops when (abs(new-old) / abs(old)) < par_threshold for 200 consecutive iterations, for all the parameters.
#' @param loss_threshold Tolerance for loss convergence. As ELBO oscillations are common in gradient based VI, we will monitor the convergence of the loss in the model,
#' the inference stops when (abs(new_loss-old_loss) / abs(old_loss)) < loss_threshold for 200 consecutive iterations.
#'
#' @return TBD
#' @export
#'
#' @examples
mobster_fit = function(x,
                       K_list=c(2, 8),
                       max_iter=1500,
                       seed_list=c(123, 1234), 
                       lr=0.01,
                       par_threshold = 0.005, 
                       loss_threshold = 0.01
                       )
{
  print(x)
  pio::pioHdr(paste0("MOBSTERm fit"))
  cat('\n')
  
  
  # Tibble
  if (is.matrix(x) | is.data.frame(x))
  {
    # Check columns
    if (!all(c("NV", "DP", "karyotype", "mutation_id", "sample_id", "purity") %in% colnames(x)))
      stop(
        "Please provide a data.frame with the following columns: mutation_id, NV, DP, sample_id, purity and karyotype."
      )
    
    cli::cli_alert_warning("Using input purity {.field {purity}}")
    
    can_work = TRUE
  }
  
  if (!can_work) {
    stop(
      "Input must be any a data.frame."
    )
  }
  
  if (is.null(x))
    return(NULL)
  
  # Check for basic input requirements
  # Here check data type and also that karyotype is unique for each sample
  mobster:::check_inputm(
    x,
    K_list,
    max_iter,
    seed_list, 
    lr,
    par_threshold, 
    loss_threshold
  )
  
  # Then: format_data_mobsterm
  # From long to wide etc etc
  formatted_data <-
    mobster:::format_data_mobsterm(x)
  
  x_wide = formatted_data$x_wid
  NV = formatted_data$NV
  DP = formatted_data$DP
  purity = formatted_data$purity
  mutation_id = formatted_data$mutation_id
  karyotype = formatted_data$karyotype
  
  # Inference here
}
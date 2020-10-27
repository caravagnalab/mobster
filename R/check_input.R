# Error checking inputs
check_input = function(x, K, samples, init, tail, epsilon, maxIter, fit.type, seed, model.selection, trace)
{
  # x 
  stopifnot(is.matrix(x) | is.data.frame(x))
  stopifnot(ncol(x) > 0)
  stopifnot('VAF' %in% colnames(x) | all(is.numeric(x$VAF)))
  
  # numerical errors
  stopifnot(all(x$VAF > 0)) 
  stopifnot(all(x$VAF < 1)) 

  # other params
  stopifnot(all(sapply(K, function(k) k >= 1 & k <= nrow(x))))
  stopifnot(samples >= 1)

  stopifnot(
    init %in%
      c('random', 'peaks')
  )
  
  stopifnot(
    tail %in% c(T, F) | is.null(tail)
  )

  stopifnot(maxIter > 0)
  stopifnot(!is.na(epsilon))
  stopifnot(epsilon > 0)
  
  stopifnot(
    fit.type %in%
      c('MM', 'MLW')
  )
  

  stopifnot(!is.null(seed) | seed > 0)

  stopifnot(model.selection %in% c('ICL', 'reICL', 'BIC', 'AIC', 'NLL'))

  stopifnot(
    trace %in% c(T, F) | is.null(tail)
  )
  
  # Check input column names, something should be reserved
  fixed_names = c("cluster", "Tail", paste0("C", 1:100))
  fixed_names = fixed_names[which(fixed_names %in% colnames(x))]
  
  if(length(fixed_names) > 0){
    stop("There are some reserved names in the input data that cannot be used, please remove or rename columns: ", paste0(fixed_names, collapse = ', '))
  }
  
}

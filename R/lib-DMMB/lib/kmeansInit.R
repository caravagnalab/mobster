#' Function for initializing parameters using 'kmeans' algorithm.
#' 
kmeansInit <- function(X, K, w, s.dot, conds, lambdas, eqProp=FALSE){
  N     <- NROW(X)                          # Number of objects
  cl    <- kmeans(X/w, K, nstart=25)        # Use Kmeans with random starts
  C.n   <- cl$cluster                       # Get the mixture components
  if (eqProp){
    pi.c <- rep(1/K, K)                     # Equal mixing proportions
  }else{
    pi.c  <- as.vector(table(C.n)/N)        # Mixing proportions using to 'kmeans'
  }
  
  C.mat <- matrix(0, nrow=N, ncol=K)        # Create matrix of cluster assignments
  for(i in 1:N){
    C.mat[i, C.n[i]] <- 1
  }
  
  # Call function for calculating the unknown parameters lambda
  lambdas <- calcLambdas(X=X, 
                         w=w, 
                         s.dot=s.dot, 
                         conds=conds, 
                         postResp=C.mat, 
                         lambdas=lambdas)
  
  # Wrap all the parameters in a list
  initParams          <- list()
  initParams$pi.c     <- pi.c 
  initParams$lambdas  <- lambdas
  
  return(initParams)
}

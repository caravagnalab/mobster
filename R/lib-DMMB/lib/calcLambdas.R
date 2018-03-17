#' Function for calculating the unknown parameters l, which
#' correspond to the clustering parameters that define the 
#' profiles of the gene in cluster k across biological 
#' condition d. Thus, the dimension of the matrix is d x K . 
#' 
#' l_{jk} = (\sum_{i}t_{ik}X_{ij.}) / (s_{j.}\sum_{i}t_{ik}X_{i..})
#' 
calcLambdas <- function(X, w, s.dot, conds, postResp, lambdas){
  d       <- NROW(lambdas)                # Number of conditions
  K       <- NCOL(lambdas)                # Number of clusters
  rhsDen  <- colSums(postResp * w)        # Calculate RHS of denominator
  for(j in 1:d) {
    denom <- s.dot[j] * rhsDen            # Multiply with LHS to get denominator value
                                          # Compute \sum_{i}\sum_{l} X_{ijl} = X.j.
    X.j.  <- rowSums(as.matrix(X[,which(conds == (unique(conds))[j])]))
                                          # Compute the numerator
    num   <- colSums(postResp * matrix(rep(X.j., K), ncol=K))
    lambdas[j,] <- num / denom            # Get the parameter value for condition d
  }
  return(lambdas)
}

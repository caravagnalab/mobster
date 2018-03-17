#' Performs EM algorithm for Binomial Mixture Models (BMMs).
#' 
#' Notation and algorithm follows Bishop's book Ch.9, "Pattern Recognition 
#' and Machina Learning". BUT, it computes the Negative Log Likelihood (NLL):
#' - ln p(X|p,pi) = - Sum_{N}(ln(Sum_{K}(pi_k * Bin(x_n|r, p_k))))
#' 
#' This method works in different modes, depending on the parameters given. 
#' If no 'theta' parameter is given, then it initializes the parameters using
#' 'kmeans' algorithm. 
##===============================================================================

bmm.EM <- function(X, r, K=2, theta, epsilon=1e-6, maxIter=2000, isLog=TRUE, isDebug=FALSE){
  
  N         <- length(X)                        # Length of the dataset
  post.resp <- matrix(0, nrow=N, ncol=K)        # Hold responsibilities
  pdf.w     <- matrix(0, nrow=N, ncol=K)        # Hold weighted PDFs
  all.NLL   <- vector(mode="numeric")           # Hold NLL for all EM iterations
  NLL       <- 1e+40                            # Initialize Negative Log Likelihood
  
  # If 'theta' parameter is empty, we initialize parameters using 'kmeans'
  if (missing(theta)){
    cl      <- kmeans(X/r, K, nstart=25)        # Use Kmeans with random starts
    C.n     <- cl$cluster                       # Get the mixture components
    p       <- as.vector(cl$centers)            # Mean for each cluster
    pi.c    <- as.vector(table(C.n)/NROW(X))    # Mixing proportions
  }else{
    p       <- theta$p                          # Mean for each cluster
    pi.c    <- theta$pi.c                       # Mixing proportions
  }
  
  if (isDebug){
    cat("Initial values:\n")
    cat("Probability:", p, "\n")
    cat("Mixing proportions:", pi.c, "\n")
    cat('Mixture components', table(C.n), '\n')
    cat("Initial NLL:", NLL, "\n")
    cat('logTransform', isLog, '\n')
  }
  
  ##=========================================
  # Run Expectation Maximization  algorithm #
  ##=========================================
  for (i in 1:maxIter){                         # Loop until convergence
    prevNLL  <- NLL                             # Store NLL to check for convergence
    
    ##===================
    #       E-Step      #
    ##===================
    if (!isLog){
      # Calculate weighted PDF of each cluster for each data point
      for (k in 1:K){
        pdf.w[,k] <- pi.c[k] * dbinom(X, size=r, prob=p[k])
      }
      Z           <- rowSums(pdf.w)             # Normalization constant
      post.resp   <- pdf.w / Z                  # Get responsibilites by normalization
      NLL         <- -sum(log(Z))               # Evaluate the NLL
    }
    else
    {
      
      # print('-- pdf.w')
      for (k in 1:K) {
        pdf.w[,k] <- log(pi.c[k]) + dbinom(X, size=r, prob=p[k], log=TRUE)
       }
      # Calculate probabilities using the logSumExp trick for numerical stability
      Z           <- apply(pdf.w, 1, logSumExp)
      post.resp   <- pdf.w - Z
      post.resp   <- apply(post.resp, 2, exp)   # Exponentiate to get actual probabilities
      NLL         <- -sum(Z)                    # Evaluate the NLL
    }
    
    
    
    ##===================
    #       M-Step      #
    ##===================
    N.k   <- colSums(post.resp)                 # Sum of responsibilities for each cluster
    pi.c  <- N.k/N                              # Update mixing proportions for each cluster
    for (k in 1:K){
      p[k]  <- post.resp[,k]%*%X/post.resp[,k]%*%r     # Update probabilities    
    NLL.Diff  <- prevNLL - NLL                  # Compute NLL difference after ith iteration
    if (NLL.Diff < 0){
      message("Negative log likelihood increases - Something is wrong!\n")
      message("Finishing EM...!")
      break
    }
    if (isDebug){
      cat("i:", i, "\t")
      cat("NLL:", NLL, "\t\t")
      cat("NLL-diff:", NLL.Diff, "\n")
    }
    all.NLL   <- c(all.NLL, NLL)                # Keep all NLL in a vector  

    # cat(i, ' -- -ll', NLL, ', eps = ', epsilon, ' @ ', NLL.Diff < epsilon, '\n')
    }
    
    
    if (NLL.Diff < epsilon){                    # Check for convergence.
      break
    }
  } #End of Expectation Maximization loop.
  
  # message("Total iterations: ", i, "\n")
  if (i == maxIter){
    message("Warning: EM did not converge with given maximum iterations (NLL.Diff = ", 
    		NLL.Diff,  
    		", epsilon = ", 
    		epsilon,")!\n\n")
  }
  
  # Add names to the estimated variables for clarity
  names(pi.c)   <- paste("Clust", 1:K)
  names(p)      <- paste("Clust", 1:K)
  
  # Cluster labels of each data point. Each data point is assigned to the cluster
  # with the highest posterior responsibility.
  labels <- unlist(apply(post.resp, 1, function(x) which(x == max(x, na.rm = TRUE))[1]))
  
  ##===========================
  # Perform model selection   #
  ##===========================
  numParams <- (K-1) + K          # Total number of parameters i.e. pi.c + p
  
  BIC <- 2*NLL + numParams*log(N) # BIC = -2*ln(L) + params*ln(N)
  AIC <- 2*NLL + 2*numParams      # AIC = -2*ln(L) + 2*params
  
  entropy <- -sum(post.resp * log(post.resp), na.rm=TRUE)
  ICL <- BIC + entropy            # Integrated Complete Likelihood criterion
  
  ##======================
  # Create a BMM object  #
  ##======================
  results           <- list()
  results$X         <- X
  results$r         <- r
  results$K         <- K
  results$N         <- N
  results$postResp  <- post.resp
  results$labels    <- labels
  results$pi.c      <- pi.c
  results$p         <- p
  results$NLL       <- NLL
  results$all.NLL   <- all.NLL
  results$BIC       <- BIC
  results$AIC       <- AIC
  results$ICL       <- ICL
  
  class(results) <- "BMM"
  return(results)
}

# library(R.utils)
# sourceDirectory("lib", modifiedOnly=FALSE)
# load('PD3851a.RData')
# 
# 
# data = data.p
#   data$total_reads = apply(data, 1, overall_read_counts)
#   data$mutated_reads = apply(data, 1, mutated_read_counts)
#   
#   result = bmm.EM(
#     X = data$mutated_reads,
#     r = matrix(data$total_reads, ncol = 1),
#     K = 2,
#     epsilon = 1e-8,
#     maxIter = 100,
#     isDebug = F
# )


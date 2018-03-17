#' Function for computing the library size normalization factors using
#' 4 different methods. These normalization factors s are estimated from
#' the data and are considered to be fixed parameters in the Poisson
#' mixture model.
#' 
#' Possible methods: 
#' 1. TC    : simply normalize by total counts
#' 2. MED   : Compute the median of the data and normalize
#' 3. DESEQ : Median ration normalization developed by DESeq package
#' 4. TMM   : Trimmed Mean of M-values normalization used in edgeR package
#'              To use this method, the edgeR package should be loaded first.
#'
normFactors <- function(X, libSize=FALSE, libType="TC"){
  q <- NCOL(X)                    # Number of variables
  s <- rep(NA, q)                 # Normalized library size for each variable
  if (libSize == FALSE) {
    s <- rep(1, q)                # If lib.size is not considered set s=1 on each variable
  }else{
    if (libType == "TC"){        # 'TC' case
      s <- colSums(X) / sum(X)
    }else if (libType == "MED"){ # 'MED' case
      s <- apply(X, MARGIN=2, FUN=median) / sum(apply(X, MARGIN=2, FUN=median))
    }else if (libType == "DESEQ"){ # 'DESEQ' case
      ## Code from DESeq function 'estimateSizeFactorsForMatrix'
      loggeomeans <- rowMeans(log(X))
      s <- apply(X, MARGIN=2, FUN=function(x) 
        exp(median((log(x)-loggeomeans)[is.finite(loggeomeans)])))
      s <- s / sum(s)
    }else if (libType == "TMM"){   # 'TMM' case
      f <- calcNormFactors(X, method = "TMM")
      s <- colSums(X)*f / sum(colSums(X)*f)
    }
  }
  return(s)
}

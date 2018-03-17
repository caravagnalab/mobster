#' This function simulates data from a Poisson mixture model, as described by 
#' Rau et al. (2011). Data are simulated with varying expression level (wi) 
#' for 4 clusters. Clusters may be simulated with "high" or "low" separation, 
#' and three different options are available for the library size setting: 
#' "equal", "A", and "B", as described by Rau et al. (2011).
#' 
#' References: Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011). 
#' Clustering highthroughput sequencing data with Poisson mixture models. Inria Research 
#' Report 7786. Available at http://hal.inria.fr/inria-00638082.

poisMixSim <- function(N=2000, libsize, separation){
  
  if(length(libsize) > 1)
    stop(paste(sQuote("libsize"), "must be equal to one of the following:", dQuote("A"), "or",
               dQuote("B"), "or", dQuote("equal"))) 
  if(libsize != "A" & libsize != "B" & libsize != "equal")
    stop(paste(sQuote("libsize"), "must be equal to one of the following:", dQuote("A"), "or",
               dQuote("B"), "or", dQuote("equal"))) 
  if(length(separation) > 1)
    stop(paste(sQuote("separation"), "must be equal to one of the following:", 
               dQuote("high"), "or", dQuote("low")))
  if(separation != "high" & separation != "low")
    stop(paste(sQuote("separation"), "must be equal to one of the following:", 
               dQuote("high"), "or", dQuote("low")))
  if(length(N) > 1)
    stop(paste(sQuote("N"), "must be a positive integer"))
  if(N <= 0)
    stop(paste(sQuote("N"), "must be a positive integer"))
  if(round(N) != N)
    stop(paste(sQuote("N"), "must be a positive integer"))
  
  ## libsize <- c("equal", "A", "B")
  if(libsize == "A" | libsize == "equal") {
    conds <- c(1, rep(2,4), rep(3, 3))
    mean.expr <- 1640
    s.norm <- c(0.112, 0.156, 0.071, 0.248, 0.165, 0.014, 0.028, 0.206)
  }else if(libsize == "B") {
    conds <- c(rep(1, 4), rep(2,2))
    mean.expr <- 1521
    s.norm <- c(0.096, 0.084, 0.253, 0.205, 0.224, 0.138)
  }
  
  r       <- table(conds)
  q       <- length(conds)          # All replicates in all conditions
  d       <- length(unique(conds))  # Number of conditions
  K.true  <- 4                      # True number of clusters
  w       <- round(rexp(N, 1/mean.expr)) # Generate data from Exp
  s.true  <- ifelse(rep(libsize, q)=="equal", rep(1/q, q), s.norm)
  s.dot.true <- rep(NA, d)          # Calculate normalized library size for each condition
  for(j in 1:d){
    s.dot.true[j] <- sum(s.true[which(conds == (unique(conds))[j])])
  }
  
  ##=================================
  #     Choosing lambda values      #
  ##=================================
  lambda.true <- matrix(NA, nrow = d, ncol = K.true)
  
  if(separation == "high"){       # High separation
    tmp <- cbind(c(1,3,5), c(5,1,3), c(3,5,1), c(5,3,1))
  }else if(separation == "low") { # Low separation
    tmp <- cbind(c(1,3,5), c(2,4,4), c(1,5,4), c(2,5,3))
  }
  
  if(libsize == "B"){             # Keep only d conditions
    tmp <- tmp[1:d,]
  }
  
  ## Choosing lambda values so that colSums(s.dot.true * lambda.true) = 1
  for(k in 1:K.true){
    lambda.norm     <- tmp[,k]/sum(tmp[,k])   # Normalized lambda values on each cluster k
    lambda.true[,k] <- lambda.norm/s.dot.true # Divide by Sj. to get appropriate lambdas
  }
  pi.true <- c(.10, .20, .30, .40)            # Mixing proportions
  
  ##===============================
  #      Simulating data          #
  ##===============================
  X <- matrix(NA, nrow=N, ncol=q)
  label.true <- rep(NA, N)
  
  for(i in 1:N){
    label.true[i] <- sample(1:K.true, prob=pi.true, size=1)
    lambda.tmp    <- rep(lambda.true[,label.true[i]], times=r)
    X[i,]         <- rmultinom(n=1, size=w[i], prob=s.true*lambda.tmp)
  }
  
  ## Remove rows with all zeros
  if(min(rowSums(X) == 0)){
    X <- X[-which(rowSums(X) == 0),]
    label.true <- label.true[-which(rowSums(X) == 0),]
    w <- w[-which(rowSums(X) == 0),]
  }
  
  simulation        <- list()
  simulation$N      <- N
  simulation$X      <- X
  simulation$w      <- w
  simulation$pi.c   <- pi.true
  simulation$conds  <- conds
  simulation$labels <- label.true
  simulation$lambda <- lambda.true
  simulation$s.true <- s.true
  simulation$s.dot.true <- s.dot.true
  
  return(simulation)
}
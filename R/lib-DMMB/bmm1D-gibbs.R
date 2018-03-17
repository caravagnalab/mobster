bmm1D.gibbs <-  function(X, r, K=2, N.Sims=10000, burnin=5000, params, 
                         ord.constr=FALSE, stephens=FALSE, logl=TRUE){
  ##================================================================
  # Function which uses Gibbs sampling so as to find the posterior #
  # of the Hierarchical Dirichlet Finite Mixture Model with        #
  # Beta prior on probability , using one dimensional dataset.     #
  ##================================================================
  
  N             <- length(X)                  # Length of the dataset
  post.resp     <- matrix(0, nrow=N, ncol=K)  # Posterior responsibilities
  pdf.w         <- matrix(0, nrow=N, ncol=K)  # PDF of each point on each cluster k
  C.n           <- matrix(0, nrow=N, ncol=K)  # Mixture components
  C.matrix      <- matrix(0, nrow=N, ncol=K)  # Total Mixture components
  NLL           <- vector(mode="numeric")     # Hold NLL for all MCMC iterations
  sum.x.k       <- vector(length=K)           # Sample sum for each cluster k 
  diff.k        <- vector(length=K)           # Sum of difference for each cluster k 
  p.draws       <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Success prob of each trial
  a.draws       <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Success prob of each trial
  b.draws       <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Success prob of each trial
  pi.draws      <- matrix(0, nrow=N.Sims-burnin, ncol=K) # Mixing Proportions
  if (stephens) # Use Stephens algorithm for relabelling MCMC outputs
    postRespArr   <- array(0, dim=c(N.Sims-burnin, N, K)) # Post resp for each MCMC run
  
  ##===============================================
  # If 'params' not defined initialize parameters #
  ##===============================================
  if (missing(params)){
    Binom       <- list()
    cl          <- kmeans(X/r, K, nstart = 25)  # Use Kmeans with random starts
    C           <- cl$cluster                   # Get the mixture components
    pi.cur      <- as.vector(table(C)/NROW(X))  # Mixing proportions
    dir.a       <- rep(1/K, K)                  # Dirichlet concentration parameter
    #dir.a       <- rep(1e-5, K)                  # Dirichlet concentration parameter
    Binom$p     <- as.vector(cl$centers)        # Binomial probability for each cluster
#    Binom$Beta  <- list(a=1, b=1)               # Initialize Beta hyperparameters
    Binom$Beta$a  <- rep(1, K)               # Initialize Beta hyperparameters
    Binom$Beta$b  <- rep(1, K)               # Initialize Beta hyperparameters
    
       
    if (ord.constr){
      ##=============================================================
      # Ordering constraint:                                        #
      #   Order all the values according to the mean parameter 'mu' #
      ##=============================================================
      Order     <- order(Binom$p)
      pi.cur    <- pi.cur[Order]
      Binom$p   <- Binom$p[Order]
    }
  }else{
    Binom       <- params$Binom
    pi.cur      <- params$pi.cur
    dir.a       <- params$dir.a
  }
  
  my.breaks = c(seq(1, N.Sims, by = N.Sims/20), N.Sims+1) # avoid overflow
  ticks = 0
  
  cat('Gibbs Sampler [. = 5%]: ')
  for (t in 1:N.Sims){  # Start Gibbs sampling
    if(t > my.breaks[ticks+1]){
      cat('.')
      ticks = ticks+1
      
      }
    
    # Compute responsibilities
    res         <- compute.resp(X, r, pdf.w, K, Binom, pi.cur, logl)
    post.resp   <- res$post.resp
    # Draw mixture components for ith simulation
    C.n         <- c.n.update(N, K, post.resp)
    # Calculate component counts of each cluster
    N.k         <- colSums(C.n)
    # Update mixing proportions using new cluster component counts
    pi.cur      <- pi.update(dir.a, N.k)
    for (k in 1:K){ 
      if (N.k[k] == 0){ # If there are no objects in the cluster
        sum.x.k[k]  <- 0
        diff.k[k]   <- 0
      }else{ # Else, calculate sum and sum of difference for each cluster
        sum.x.k[k]  <- sum(X[C.n[,k] == 1])
        diff.k[k]   <- sum(r[C.n[,k] == 1] - X[C.n[,k] == 1])
      }
    }
    # Update posterior mean
    obj_out <- p.update(K, Binom, sum.x.k, diff.k)
    Binom$p <- obj_out$p.posterior
    alpha.n <- obj_out$alpha.n
    beta.n <- obj_out$beta.n
    
    if (ord.constr){
      ##=============================================================
      # Ordering constraint:                                        #
      #   Order all the values according to the mean parameter 'mu' #
      ##=============================================================
      Order     <- order(Binom$p)
      pi.cur    <- pi.cur[Order]
      Binom$p   <- Binom$p[Order]
    }
    
    # Keep only the simulations after the burned in period has passed
    if (t > burnin){
      NLL                     <- c(NLL, res$NLL)
      C.matrix                <- C.matrix + C.n
      pi.draws[t - burnin,]   <- pi.cur
      p.draws[t - burnin,]    <- Binom$p
      a.draws[t - burnin,]    <- alpha.n
      b.draws[t - burnin,]    <- beta.n
      if (stephens) # Use Stephens algorithm for relabelling MCMC outputs
        postRespArr[t-burnin, , ] <- post.resp
    }
  }
  
  cat(' [OK]\n')
  
  # Object to keep input data
  dat         <- NULL
  dat$X       <- X
  dat$r       <- r
  dat$K       <- K
  dat$N       <- N
  dat$N.Sims  <- N.Sims
  dat$burnin  <- burnin
  
  # Object to hold all the MCMC draws
  draws       <- NULL
  draws$pi    <- pi.draws
  draws$p     <- p.draws
  draws$a     <- a.draws
  draws$b     <- b.draws
  
  # Object to hold the summaries for the parameters
  summary     <- NULL
  summary$pi  <- apply(pi.draws, 2, mean)   # Expected value of mix. prop.
  summary$p   <- apply(p.draws, 2, mean)    # Expected value of mean
  summary$a   <- apply(a.draws, 2, mean)    # Expected value of mean
  summary$b   <- apply(b.draws, 2, mean)    # Expected value of mean
  summary$C   <- C.matrix / (N.Sims-burnin) # Convert C.matrix to probs
  summary$NLL <- NLL
  if (stephens) # Use Stephens algorithm for relabelling MCMC outputs
    summary$PR <- postRespArr
  
  # Object to hold the credible intervals for the parameters
  cred.interv     <- NULL
  cred.interv$pi  <- apply(pi.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  cred.interv$p   <- apply(p.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  cred.interv$a   <- apply(a.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  cred.interv$b   <- apply(b.draws, 2, quantile, prob=c(0.025, 0.5, 0.975))
  
  return(list(dat=dat, draws=draws, summary=summary, cred.interv=cred.interv))
}

# Compute the responsibilities
compute.resp <- function(X, r, pdf.w, K, Binom, pi.cur, logl){
  if (logl){
    for (k in 1:K){ # Calculate the PDF of each cluster for each data point
      pdf.w[,k] <- log(pi.cur[k]) + dbinom(X, size=r, prob=Binom$p[k], log=TRUE)
    }
    # Calculate probabilities using the logSumExp trick for numerical stability
    Z           <- apply(pdf.w, 1, logSumExp)
    post.resp   <- pdf.w - Z
    post.resp   <- apply(post.resp, 2, exp)   # Exponentiate to get actual probabilities
    NLL         <- -sum(Z)                    # Evaluate the NLL
  }else{
    for (k in 1:K){ # Calculate the PDF of each cluster for each data point
      pdf.w[,k] <- pi.cur[k] * dbinom(X, size=r, prob=Binom$p[k])
    }
    Z           <- rowSums(pdf.w)             # Normalization constant
    post.resp   <- pdf.w / Z                  # Get responsibilites by normalization
    NLL         <- -sum(log(Z))               # Evaluate the NLL
  }
  return(list(post.resp=post.resp, NLL=NLL))
}

# Update the mixture components 
c.n.update <- function(N, K, post.resp){
  c.i.draw <- matrix(0, nrow=N, ncol=K)
  for (i in 1:N){ # Sample one point from a multinomial i.e. ~ Discrete
    c.i.draw[i,] = rmultinom(1, 1, post.resp[i,])
  }
  return(c.i.draw)
}

# Update the mixing proportions
pi.update <- function(dir.a, N.k){
  a.n <- dir.a + N.k
  return(as.vector(rdirichlet(n=1, alpha=a.n))) # Sample from Dirichlet
}

# Update the posterior mean
p.update <- function(K, Binom, sum.x.k, diff.k){
  p.posterior = alpha.n = beta.n <- vector(length=K)
  for (k in 1:K){
    alpha.n[k]   <- Binom$Beta$a[k] + sum.x.k[k]
    beta.n[k]    <- Binom$Beta$b[k] + diff.k[k]
    p.posterior[k] <- rbeta(1, shape1=alpha.n[k], shape2=beta.n[k])
  }
  return(list(p.posterior = p.posterior, alpha.n = alpha.n, beta.n = beta.n))
}
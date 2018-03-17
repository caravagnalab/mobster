##=============================================
# Main script that will be called to read and #
# generate data                               #
##============================================= 

##=============================================
# Generate data:  mixture of 1D Gaussians     #
##=============================================
gen.gaussian <- function(N=500, K=2, pi.c=c(0.6,0.4), mus=c(10,20), stds=c(2,3)){
  # Sample the N components according to their mixing proportions
  components <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Generate N random samples from a 1D-Gaussian, with the corresponding weights
  samples    <- rnorm(n=N,mean=mus[components],sd=stds[components])
  
  #plot(density(samples),main="Density Estimate of the Mixture Model")
  return (samples)
}

##=============================================
# Generate data:  mixture of MV Gaussians     #
##=============================================
gen.MV.gaussian <- function(N=500, K=2, pi.c=c(0.6,0.4), mus=matrix(c(1.5,4,4,1.5),2,2), S=list()){
  if (missing(S)){
    S[[1]] <- matrix(c(1,0,0,1),2,2)
    S[[2]] <- matrix(c(1,.7,.7,1),2,2)
  }
  # Sample the N components according to their mixing proportions
  components <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Rows of the mean matrix are clusters and columns are the dimensionality of the data
  D         <- NCOL(mus)
  # Generate N random samples from a MV-Gaussian, with the corresponding weights
  samples   <- matrix(0, nrow=N, ncol=D)
  for (i in 1:N)
    samples[i,]    <- rmvnorm(1,mean=mus[components[i],], sigma=S[[components[i]]])
  return (samples)
}

##=============================================
# Generate data: mixture of 1D Poissons       #
##=============================================
gen.poisson <- function(N=500, K=2, pi.c=c(0.6,0.4), lambdas=c(8,22)){
  # Sample the N components according to their mixing proportions
  components <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Generate N random samples from a Gaussian, with the corresponding weights
  samples    <- rpois(n=N,lambda=lambdas[components])
  
  return (samples)
}

##=============================================
# Generate data:  mixture of MV Poissons      #
##=============================================
gen.MV.poisson <- function(N=500, K=2, pi.c=c(0.6,0.4), lambdas=matrix(c(16,10,7,30),2,2)){
  # Sample the N components according to their mixing proportions
  components  <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Rows of the mean matrix are clusters and columns are the dimensionality of the data
  D           <- NCOL(lambdas)
  # Generate N random samples from a MV-Gaussian, with the corresponding weights
  samples     <- matrix(0, nrow=N, ncol=D)
  for (i in 1:N)
    samples[i,]    <- rpois(2,lambda=lambdas[components[i],])
  return (samples)
}

##============================================
# Generate data: mixture of 1D Binomials     #
##============================================
gen.binomial <- function(N=500, K=2, pi.c=c(0.6,0.4), p=c(0.3,0.8), r){
  # Sample the N components according to their mixing proportions
  components <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Generate N random samples from a Gaussian, with the corresponding weights
  samples    <- rbinom(n=N, size=r, prob=p[components])
  return (samples)
}

##=============================================
# Generate data:  mixture of MV Poissons      #
##=============================================
gen.MV.binomial <- function(N=500, K=2, pi.c=c(0.6,0.4), p=matrix(c(0.1,0.8,0.5,0.2),2,2), r){
  # Sample the N components according to their mixing proportions
  components  <- sample(1:K, prob=pi.c, size=N, replace=TRUE)
  # Rows of the mean matrix are clusters and columns are the dimensionality of the data
  D           <- NCOL(p)
  # Generate N random samples from a MV-Gaussian, with the corresponding weights
  samples     <- matrix(0, nrow=N, ncol=D)
  for (i in 1:N)
    samples[i,]    <- rbinom(2, size=r[i,], prob=p[components[i],])
  return (samples)
}

##===================================================
# Generate data: Mixture of 3 methylation profiles  #
##===================================================
gen.meth.data <- function(N=300, pi.c=c(0.45, 0.35, 0.2), maxL=25, minLoc=-100, maxLoc=100){
  # Create a list to store for each methylation region its corresponding data
  X       <- list()
  # For each of the N objects
  for (i in 1:N){
    # L is the number of Cs found in the ith region
    L <- rbinom(n=1, size=maxL, prob=.8)
    X[[i]] <- matrix(0, nrow=3, ncol=L)
    # Randomly sample locations for the Cs and scale them, so the data lie in the (-2,2) region
    X[[i]][1,] <- sort(sample(minLoc:maxLoc,L))/max(abs(minLoc), abs(maxLoc))
    
    if (i < N * pi.c[1]){   # First methylation profile
      lb <- round(L/2.5)
      mb <- round(L/5)
      
      X[[i]][2,1:lb] <- rbinom(n=lb, size=20, prob=.9)
      repeat{
        X[[i]][3,1:lb] <- rbinom(n=lb, size=14, prob=.9)
        if(all(X[[i]][2,1:lb] > X[[i]][3,1:lb]))
          break
      }
      
      X[[i]][2,(lb+1):(lb+mb)] <- rbinom(n=mb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1):(lb+mb)] <- rbinom(n=mb, size=7, prob=.9)
        if (all(X[[i]][2,(lb+1):(lb+mb)] > X[[i]][3,(lb+1):(lb+mb)]))  
          break
      }
      
      X[[i]][2,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=3, prob=.8)
        if (all(X[[i]][2,(lb+1+mb):L] > X[[i]][3,(lb+1+mb):L]))
          break
      }
    }else if (i < (N * pi.c[2] + N * pi.c[1])){     # Second methylation profile
      lb <- round(L/2.5)
      mb <- round(L/4)
      
      X[[i]][2,1:lb] <- rbinom(n=lb, size=20, prob=.9)
      repeat{
        X[[i]][3,1:lb] <- rbinom(n=lb, size=3, prob=.8)
        if(all(X[[i]][2,1:lb] > X[[i]][3,1:lb]))
          break
      }
      
      X[[i]][2,(lb+1):(lb+mb)] <- rbinom(n=mb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1):(lb+mb)] <- rbinom(n=mb, size=7, prob=.9)
        if (all(X[[i]][2,(lb+1):(lb+mb)] > X[[i]][3,(lb+1):(lb+mb)]))
          break
      }
      
      X[[i]][2,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=14, prob=.9)
        if (all(X[[i]][2,(lb+1+mb):L] > X[[i]][3,(lb+1+mb):L]))
          break
      }
    }else{                  # Third methylation profile
      lb <- round(L/3)
      mb <- round(L/3)
      
      X[[i]][2,1:lb] <- rbinom(n=lb, size=20, prob=.9)
      repeat{
        X[[i]][3,1:lb] <- rbinom(n=lb, size=4, prob=.9)
        if(all(X[[i]][2,1:lb] > X[[i]][3,1:lb]))
          break
      }
      
      X[[i]][2,(lb+1):(lb+mb)] <- rbinom(n=mb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1):(lb+mb)] <- rbinom(n=mb, size=14, prob=.9)
        if (all(X[[i]][2,(lb+1):(lb+mb)] > X[[i]][3,(lb+1):(lb+mb)]))  
          break
      }
      
      X[[i]][2,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=4, prob=.9)
        if (all(X[[i]][2,(lb+1+mb):L] > X[[i]][3,(lb+1+mb):L]))
          break
      }
    }
  }
  return(X)
}


##===================================================
# Generate data: Mixture of 3 methylation profiles  #
##===================================================
gen.meth.data2 <- function(N=300, pi.c=c(0.45, 0.35, 0.2), maxL=25, minLoc=-100, maxLoc=100){
  # Create a list to store for each methylation region its corresponding data
  X       <- list()
  # For each of the N objects
  for (i in 1:N){
    # L is the number of Cs found in the ith region
    L <- rbinom(n=1, size=maxL, prob=.8)
    X[[i]] <- matrix(0, nrow=3, ncol=L)
    # Randomly sample locations for the Cs and scale them, so the data lie in the (-2,2) region
    X[[i]][1,] <- sort(sample(minLoc:maxLoc,L))/max(abs(minLoc), abs(maxLoc))
    
    if (i < N * pi.c[1]){   # First methylation profile
      lb <- round(L/3)
      
      X[[i]][2,1:lb] <- rbinom(n=lb, size=20, prob=.9)
      repeat{
        X[[i]][3,1:lb] <- rbinom(n=lb, size=14, prob=.9)
        if(all(X[[i]][2,1:lb] > X[[i]][3,1:lb]))
          break
      }
      
      X[[i]][2,(lb+1):L] <- rbinom(n=L-lb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1):L] <- rbinom(n=L-lb, size=3, prob=.9)
        if (all(X[[i]][2,(lb+1):L] > X[[i]][3,(lb+1):L]))  
          break
      }
    }else if (i < (N * pi.c[2] + N * pi.c[1])){     # Second methylation profile
      lb <- round(L/1.5)
      
      X[[i]][2,1:lb] <- rbinom(n=lb, size=20, prob=.9)
      repeat{
        X[[i]][3,1:lb] <- rbinom(n=lb, size=3, prob=.8)
        if(all(X[[i]][2,1:lb] > X[[i]][3,1:lb]))
          break
      }
      
      X[[i]][2,(lb+1):L] <- rbinom(n=L-lb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1):L] <- rbinom(n=L-lb, size=14, prob=.9)
        if (all(X[[i]][2,(lb+1):L] > X[[i]][3,(lb+1):L]))
          break
      }
    }else{                  # Third methylation profile
      lb <- round(L/3)
      mb <- round(L/3)
      
      X[[i]][2,1:lb] <- rbinom(n=lb, size=20, prob=.9)
      repeat{
        X[[i]][3,1:lb] <- rbinom(n=lb, size=3, prob=.9)
        if(all(X[[i]][2,1:lb] > X[[i]][3,1:lb]))
          break
      }
      
      X[[i]][2,(lb+1):(lb+mb)] <- rbinom(n=mb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1):(lb+mb)] <- rbinom(n=mb, size=14, prob=.9)
        if (all(X[[i]][2,(lb+1):(lb+mb)] > X[[i]][3,(lb+1):(lb+mb)]))  
          break
      }
      
      X[[i]][2,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=3, prob=.9)
        if (all(X[[i]][2,(lb+1+mb):L] > X[[i]][3,(lb+1+mb):L]))
          break
      }
    }
  }
  return(X)
}


##===================================================
# Generate data: Mixture of 3 methylation profiles  #
##===================================================
gen.meth.data3 <- function(N=300, pi.c=c(0.45, 0.35, 0.2), maxL=25, minLoc=-100, maxLoc=100){
  # Create a list to store for each methylation region its corresponding data
  X       <- list()
  # For each of the N objects
  for (i in 1:N){
    # L is the number of Cs found in the ith region
    L <- rbinom(n=1, size=maxL, prob=.8)
    X[[i]] <- matrix(0, nrow=3, ncol=L)
    # Randomly sample locations for the Cs and scale them, so the data lie in the (-1,1) region
    X[[i]][1,] <- sort(sample(minLoc:maxLoc,L))/max(abs(minLoc), abs(maxLoc))
    
    if (i < N * pi.c[1]){   # First methylation profile
      lb <- round(L/4)
      
      X[[i]][2,1:lb] <- rbinom(n=lb, size=20, prob=.9)
      repeat{
        X[[i]][3,1:lb] <- rbinom(n=lb, size=14, prob=.9)
        if(all(X[[i]][2,1:lb] > X[[i]][3,1:lb]))
          break
      }
      
      X[[i]][2,(lb+1):L] <- rbinom(n=L-lb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1):L] <- rbinom(n=L-lb, size=2, prob=.9)
        if (all(X[[i]][2,(lb+1):L] > X[[i]][3,(lb+1):L]))  
          break
      }
    }else if (i < (N * pi.c[2] + N * pi.c[1])){     # Second methylation profile
      lb <- round(L/1.5)
      
      X[[i]][2,1:lb] <- rbinom(n=lb, size=20, prob=.9)
      repeat{
        X[[i]][3,1:lb] <- rbinom(n=lb, size=2, prob=.8)
        if(all(X[[i]][2,1:lb] > X[[i]][3,1:lb]))
          break
      }
      
      X[[i]][2,(lb+1):L] <- rbinom(n=L-lb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1):L] <- rbinom(n=L-lb, size=14, prob=.9)
        if (all(X[[i]][2,(lb+1):L] > X[[i]][3,(lb+1):L]))
          break
      }
    }else{                  # Third methylation profile
      lb <- round(L/2.5)
      mb <- round(L/3.5)
      
      X[[i]][2,1:lb] <- rbinom(n=lb, size=20, prob=.9)
      repeat{
        X[[i]][3,1:lb] <- rbinom(n=lb, size=2, prob=.9)
        if(all(X[[i]][2,1:lb] > X[[i]][3,1:lb]))
          break
      }
      
      X[[i]][2,(lb+1):(lb+mb)] <- rbinom(n=mb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1):(lb+mb)] <- rbinom(n=mb, size=14, prob=.9)
        if (all(X[[i]][2,(lb+1):(lb+mb)] > X[[i]][3,(lb+1):(lb+mb)]))  
          break
      }
      
      X[[i]][2,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=20, prob=.9)
      repeat{
        X[[i]][3,(lb+1+mb):L] <- rbinom(n=L-mb-lb, size=2, prob=.9)
        if (all(X[[i]][2,(lb+1+mb):L] > X[[i]][3,(lb+1+mb):L]))
          break
      }
    }
  }
  return(X)
}


##=====================================
# Read data from Old Faithful Dataset #
##=====================================
read.faithful1D <- function(dim=2){
  X <-  faithful[,dim]
}

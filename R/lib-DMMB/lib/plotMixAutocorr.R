plotMixAutocorr <- function(mcmc.obj, nams){
  # Plot autocorrelation functions
  K     <- mcmc.obj$dat$K           # Number of clusters
  pars  <- length(mcmc.obj$draws)   # Number of parameters
  if (missing(nams)){
    nams  <- names(mcmc.obj$draws)  # Store names for showing in plots
  }
  
  for (i in 1:pars){            # Plot each parameter for ...
    par(mfrow=c(1, K))          # Set plot frames
    for (j in 1:K) {            # each cluster k
      if (is.list(mcmc.obj$draws[[i]])){
        if (length(mcmc.obj$draws[[i]]) == 0){
          next
        }else{
          for (m in 1:length(mcmc.obj$draws[[i]]))
            acf(mcmc.obj$draws[[i]][[m]][,j], lag.max=100, col=2, main=paste(nams[i],j,sep="_"))
        }
      }else if (is.matrix(mcmc.obj$draws[[i]])){
        acf(mcmc.obj$draws[[i]][,j], lag.max=100, col=2, main=paste(nams[i],j,sep="_"))
      }else{
        acf(mcmc.obj$draws[[i]], lag.max=100, col=2, main=paste(nams[i],j,sep="_"))
      }
    }
  }
}

plotMixTrace <- function(mcmc.obj, nams){
  # Plot traces of sampled parameters during MCMC
  K     <- mcmc.obj$dat$K           # Number of clusters
  pars  <- length(mcmc.obj$draws)   # Number of parameters
  if (missing(nams)){
    nams  <- names(mcmc.obj$draws)  # Store names for showing in plots
  }
  par(mfrow=c(pars, K))             # Set plot frames
  
  for (i in 1:pars){            # Plot each parameter for ...
    par(mfrow=c(1, K))          # Set plot frames
    for (j in 1:K) {            # each cluster k
      if (is.list(mcmc.obj$draws[[i]])){
        if (length(mcmc.obj$draws[[i]]) == 0){
          next
        }else{
          for (m in 1:length(mcmc.obj$draws[[i]]))
            plot(mcmc.obj$draws[[i]][[m]][,j], main=paste(nams[i],j,sep="_"), type="l")
        }
      }else if (is.matrix(mcmc.obj$draws[[i]])){
        plot(mcmc.obj$draws[[i]][,j], main=paste(nams[i],j,sep="_"), type="l")
      }else{
        plot(mcmc.obj$draws[[i]], main=paste(nams[i],j,sep="_"), type="l")
      }
    }
  }
}

# Functions to prepare Binomial fit data
overall_read_counts = function(x) {
	return(
		as.numeric(x['expt_a_count']) + 
		as.numeric(x['expt_c_count']) + 
		as.numeric(x['expt_t_count']) + 
		as.numeric(x['expt_g_count']))
}

mutated_read_counts = function(x){
	return(round(as.numeric(x['mut_allele_proportion']) * overall_read_counts(x)))
}

# Wrappers for EM
runEM = function(data, k, epsilon = 1e-3, maxIter = 2000, synthetic = FALSE) # one-run of the EM
{
	if(! synthetic)
	{
		data$total_reads = apply(data, 1, overall_read_counts)
		data$mutated_reads = apply(data, 1, mutated_read_counts)
  	}
  
  result = bmm.EM(
	  X = data$mutated_reads,
	  r = matrix(data$total_reads, ncol = 1),
	  K = k,
	  epsilon = epsilon,
	  maxIter = maxIter
  )
	
	return(result)
}

# N runs of the EM, with a table created as output and some plot to show the likelihood trends
sample.BIC.EM = function(data, K, N, epsilon = 1e-3, maxIter = 2000,  synthetic = FALSE, parallel = FALSE)
{
  BIC.values = matrix(0, ncol = length(K), nrow = N)
  DRES = NULL
  
  for(k in K)
  { 
    RES = NULL
    
    if(!parallel)
    {
      cat('K = ',k, ' ')
      for(i in 1:N) {
        cat('.')
        result = runEM(data, k, epsilon, maxIter, synthetic)
        BIC.values[i, k] = result$BIC
        RES = append(RES, list(result))  
      }
      cat('[OK]\n') 
    }
    else
    {
      require(parallel)
      require(doParallel)
      
      # set the number of cores to be used (80% of available ones)
      cores = as.integer(.8 * (detectCores() - 1))
      if (cores < 1) cores = 1
      
      # setup the parallelization to perform the bootstrap
      cl = makeCluster(cores)
      registerDoParallel(cl)
      
      message('Using parallel with ', cores, "/", detectCores(), ' cores')
      
      r = foreach(num = 1:N, .export = ls(globalenv())) %dopar% 
      {
        result = runEM(data, k, epsilon, maxIter, synthetic)
      }
      
      RES = r
      
      DRES = append(DRES, list(RES))  
      
      for(i in 1:N) 
        BIC.values[i, k] = RES[[i]]$BIC
    
      stopCluster(cl)
      }
  }
  
  colnames(BIC.values) = paste('K=', K)
  rownames(BIC.values) = 1:N
  
  return(list(BIC.values=BIC.values, samples = DRES))
}




############### GIBBS WRAPPERS



# Wrap Gibbs sampler
runGibbs = function(data, k, N.Sims = 20000,  burnin = 13000, synthetic = FALSE, ...)
{
	if(! synthetic)
	{
  	data$total_reads = apply(data, 1, overall_read_counts)
  	data$mutated_reads = apply(data, 1, mutated_read_counts)
	}
	  
  result = bmm1D.gibbs(
	  X = data$mutated_reads,
	  r = data$total_reads, # this sampler here gets a vector, not a matrix as the EM one
	  K = k,
	  N.Sims = N.Sims,
	  burnin = burnin,
	  ...
  )
  
  message(
    paste('\n(Average) Mixture proportions\n', paste(round(result$summary$pi, 2), collapse = ','),
        '\n(Average) Mixture parameters\n', paste(round(result$summary$p, 2), collapse = ','),
        '\n(Average) Posterior a, b for Beta(a,b)\n', 
        paste(round(result$summary$a, 2), collapse = ','), '\n',
        paste(round(result$summary$b, 2), collapse = ','),
        '\n(Average) NLL:', round(mean(result$summary$NLL), 2)
      )
  )
	
	return(result)
}

sample.proportions.Gibbs = function(data, K, N, concentration = 1e-5,  N.Sims = 20000,  burnin = 13000,   synthetic = FALSE, parallel = FALSE)
{
  create.initialCond = function()
  {
    if(! synthetic)
    {
      data$total_reads = apply(data, 1, overall_read_counts)
      data$mutated_reads = apply(data, 1, mutated_read_counts)
    }

    Binom       <- list()
    cl          <- kmeans(data$mutated_reads/data$total_reads, K, nstart = 1)   # Use Kmeans with no starts
    C           <- cl$cluster                   # Get the mixture components
    pi.cur      <- as.vector(table(C)/NROW(data$mutated_reads))  # Mixing proportions
    dir.a       <- rep(concentration, K)                  # Dirichlet concentration parameter
    Binom$p     <- as.vector(cl$centers)        # Binomial probability for each cluster
    Binom$Beta  <- list(a=1, b=1)               # Initialize Beta hyperparameters
    
    return(list(Binom=Binom, pi.cur=pi.cur, dir.a=dir.a)) 
  }
  
  proportion.values = matrix(0, ncol = K, nrow = N)
  colnames(proportion.values) = paste('pi_', 1:ncol(proportion.values))
  rownames(proportion.values) = 1:N
  
  cat('Sampling', N, 'chains with', N.Sims, 'steps.\n')
  
  RES = NULL
  
  if(!parallel) {
    for(i in 1:N) 
    {
      result = runGibbs(data = data,  k = K, N.Sims = N.Sims, burnin = burnin, synthetic = synthetic, params = create.initialCond())
      RES = append(RES, list(result))  
      
      proportion.values[i, ] = sort(result$summary$pi, decreasing = TRUE)
    }
  }
  else
  {
    require(parallel)
    require(doParallel)

    # set the number of cores to be used (80% of available ones)
    cores = as.integer(.8 * (detectCores() - 1))
    if (cores < 1) cores = 1
    
    # setup the parallelization to perform the bootstrap
    cl = makeCluster(cores)
    registerDoParallel(cl)

    message('Using parallel with ', cores, "/", detectCores(), ' cores')
    
    r = foreach(num = 1:N, 
                .packages = "MCMCpack",
                .export = ls(globalenv())) %dopar% {
      
        result = runGibbs(data = data,  k = K, N.Sims = N.Sims, burnin = burnin, synthetic = synthetic, params = create.initialCond())
      }
  
    # recollect after dopar  
    RES = r
    for(i in 1:N) 
      proportion.values[i, ] = sort(RES[[i]]$summary$pi)
    
    stopCluster(cl)
    print(proportion.values)
    
  }
  
  # We take the median of the proportion, if label switching happens rarely, removing outliers is good
  require(robustbase)
  medians = colMedians(proportion.values)
  averages = colMeans(proportion.values)
  
  # Visualize the logLik (NLL)
  require(RColorBrewer)
  cols = colorRampPalette(brewer.pal(n = 9, 'Set1'))(K) 

  for(i in 1:K) 
  {
    #dens = density(proportion.values[, i])
    #dens$y = dens$y/max(dens$y)
    #if(i == 1) plot(dens, col = cols[i], xlim = c(0,1))    
    #else lines(dens, col = cols[i], xlim = c(0,1))
    hist(proportion.values[, i], main = paste("pi_", i, ': median ', round(medians[i], 2)), col = cols[i])
  }

  return(list(proportion.values = proportion.values, medians = medians, averages = averages, samples = RES))

}
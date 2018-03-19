# Optimal BIC scores as a function of K. With visualization and 2 selection criteria that resembles the ones used for Gap statistics.
optimalBIC =  function(BIC.values, coeff = 0.05, plot = TRUE, logscale = TRUE)
{
  plot.BIC.values = BIC.values
  if(logscale) {
    plot.BIC.values = log(plot.BIC.values)
    message('Boxplot rescaled in log space')
  }
  means = colMeans(-BIC.values)
  
	global.optimal = which(means == max(means))
	cat('Means:', round(means,2), '\n')
	cat('[Global Max] K =', global.optimal, 'with average BIC', means[global.optimal])

	dderiv = sapply(1:(length(means)-1), function(i) means[i+1]-means[i])
	
	largest.gap = max(dderiv) * coeff
	deriv.optimal = min(which(dderiv < largest.gap, arr.ind = TRUE))
	
	cat('\n\nDiscrete Derivatives:', round(dderiv, 2), 
	    '\n', coeff * 100, '% of max =',largest.gap,
	    '\n[Derivative Max] K =', deriv.optimal , 
	    'with average BIC', -BIC.values[deriv.optimal])

  diff = round(abs(means[deriv.optimal] - means[global.optimal]), 2)
	cat("\n\nAbsolute difference among average BIC scores:", diff) 
	
	if(plot)	
	{	
		boxplot(-plot.BIC.values, 
			sub = paste('Glob. Max K=', global.optimal, ', Der. Max K=', deriv.optimal,', diff =', diff, sep = ''), 
			main = 'BIC scores', 
			xlab = 'K',
			ylab = ifelse(logscale, 'log(BIC)', 'BIC'),
		)
	}		
	return(list(deriv.optimal = deriv.optimal, global.optimal = global.optimal))
}
	
plot.NLL.optimalK = function(results, collapse = FALSE)
{
  # Visualize the logLik (NLL)
  if(collapse) par(mfrow = c(max(K)/3, max(K)/3), mar = c(2,2,2,2))
  
  require(RColorBrewer)
  K = ncol(results$BIC.values)
  N = length(results$samples)/K
  
  print(K)
  print(N)
  
  cols = colorRampPalette(brewer.pal(n = 9, 'Set1'))(N)
  for(i in 1:(N)) 
  {
    
    if(i == 1) plot(results$samples[[i * K]]$all.NLL, type = 'l', main = paste("K =", k), col = cols[i], ylab = 'NLL', xlab = 'Step')
    else lines(results$samples[[i * K]]$all.NLL, col = cols[i])
  }
  
}
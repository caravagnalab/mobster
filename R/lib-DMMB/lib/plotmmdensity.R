# Binomial MM density plot
plotmmdensity =  function(gibbs,  trials = 100, palette = 'Set1', is.norm = TRUE)
{
	# average of mixing proportions and parameters of the mixture
	pi = gibbs$summary $pi 
	p = gibbs$summary $p 	
	K = gibbs$dat$K

	if(! is.norm)  pi = pi /sum(pi)

	# Trials, 1 to ...
	x = 1:trials

	# Simulate a Binomial sampling with "trials", and report the probability of having "x" successes
	# according to the mixture
	density = NULL
	mixture = 0.0
	for(i in 1:K)
	{
		# Calculate the estimated density (rescaled for the proprortions)
	  mixture.points = dbinom(x, size = trials, prob = gibbs$summary$p[i]) * gibbs$summary$pi[i]
		density = append(density, list(mixture.points))
		mixture <- mixture + mixture.points
	}
	
	# plot the overall mixture, and the densities
	require(RColorBrewer)
	colors = colorRampPalette(brewer.pal(9, palette))(K)
	  
	plot(mixture, 
	     col = 'darkred', 
	     type = 'l', 
	     lwd = 2, 
	     ylab = 'mixture density',
	     xlab = 'number of successes',
	     main = paste('Probability of success for N = ', trials, 'trials'))
	
	legend(
	  'topright', 
	  paste('p=', round(gibbs$summary$p, 2), '(prop.', round(gibbs$summary$pi, 2),')'), 
	  col = colors, lty = 1, bg = 'grey96')

	for(i in 1:K) lines(x, density[[i]], col= colors[i], lty = 2)
	#lines(x, mixture, col = 'darkred', lwd = 2)
}

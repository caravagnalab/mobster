# Gaussian KDE of posterior density for univariate Dirichlet MM with Binomial likelihood.
# gibbs -- structure returned by Andreas' Gibbs sample 
# kde.density.smooth -- density smoothing (suggested <1)
labposterior =  function(gibbs,  is.norm = TRUE, collapse.plots = FALSE)
{
	pi = gibbs$draws$pi #(mixing proportions, usually written as \pi_j)
	p = gibbs$draws$p 	# (parameter of the mixture)
	K = gibbs$dat$K

	if(! is.norm)  pi = t(apply(pi, 1, function(x) x/sum(x)))

	# # Test for label switching -- plot all groups
	# bp.pi.cols = sample.coda[ ,pi.cols]
	# bp.kappa.cols = sample.coda[ , kappa.cols]
	# bp.kappa.cols = t(apply(bp.kappa.cols, 1, function(x) x/sum(x)))

	if(collapse.plots) par(mfrow = c(K, 2))
	for(i in 1:K)
	{
	  if(!collapse.plots) layout(matrix(c(1, 1, 1, 2), nrow = 1))
	    
		ymax = max(max(pi[, i]), p[, i])		
		avg = round(mean(pi[, i]), 2)
		
		
		plot(1:nrow(pi), pi[, i],  
			type = 'l', 
			main = paste('Cluster k =', i, ' (avg pi = ', avg, ')'), 
			ylim = c(0, ymax), 
			xlab = 'MCMC step',
			ylab = '',
			col = 'red')
		lines(1:nrow(p), p[, i], col = 'darkblue')
				
		# legend('topright', c(paste('p_', i), paste('pi_', i)), col = c('red', 'darkblue'), lty = 1, bg = 'grey96')
		legend('topright', 
			expression(pi[k], p[k]), 
			col = c('red', 'darkblue'), lty = 1, bg = 'grey96')
  

		d <- density(pi[, i], adjust = 0.05)
		plot(d, main="KDE", col = "black", bty = 'n', ylab = '', xlab = '', yaxt = 'n')
		
	}
	
		
	# # check a bit of this density	
		
	# par(mfrow = c(4,1))

	# index = 1
	# pi = sample.coda[index,pi.cols]
	# k = sample.coda[index,kappa.cols] / sum(sample.coda[index,kappa.cols])

	# plot(pi, col = 'red', pch = 19, main = 'Mixture parameter (Binomial: p ~ success) -- one Gibbs sample')
	# lines(x = 1:length(pi), y = pi, col = 'red', pch=16, type = 'h')
	# points(pi)

	# plot(k, col = 'lightblue', pch = 19, main = 'Mixture proportions (Categorical: pi_k) -- one Gibbs sample')
	# lines(x = 1:length(pi), y = k, col = 'lightblue', pch=16, type = 'h')
	# points(k)

	# # this density estimation does not depend on from/to, but rather on adjust
	# # 	- adjust: the bandwidth used is actually adjust*bw. This makes it easy to specify values like ‘half the default’ bandwidth

	# dens =  density(
		# c(sample.coda[1,pi.cols]), 
		# from = 0, to = 1,
		# weights = c(sample.coda[1,kappa.cols]) / sum(c(sample.coda[1,kappa.cols])), adjust = density.smooth)
	
	# py = dens$y
	# px = dens$x
		
	# plot(dens, 'KDE with Gaussian and adjust = 0.05')
	# points(px, py, pch = 19, col = 'orange')	
	# lines(px,py)

	# dens =  density(
		# c(sample.coda[1,pi.cols]), 
		# from = 0, to = 1,
		# weights = c(sample.coda[1,kappa.cols]) / sum(c(sample.coda[1,kappa.cols])), adjust = 1)
	
	# py = dens$y
	# px = dens$x
		
	# plot(dens,  'KDE with Gaussian and adjust = 1')
	# points(px, py, pch = 19, col = 'orange')	
	# lines(px,py)
	
	# dev.copy2pdf(file='Plot-KDE.pdf')

	
	# for (i in 1:dim(sample.coda)[1]) 
	# {
		# posterior.samples[,i] <- density(
			# c(sample.coda[i,pi.cols]), 
			# weights=c(sample.coda[i,kappa.cols]) / sum(c(sample.coda[i,kappa.cols])), adjust=density.smooth, from=0,to=1)$y
	# }
	# ncol(posterior.samples)
	# nrow(posterior.samples)
	# head(posterior.samples)
	
	# plot(posterior.samples[, 1]) 
	# for (i in 1: ncol(posterior.samples))
		# lines(posterior.samples[, i]) 


}

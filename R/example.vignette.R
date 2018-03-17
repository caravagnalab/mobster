# # library(dbpmm)
# #
# B1 = rbeta(100, 50, 50) # clonal cluster
# TAIL = sads::rpareto(400, shape = 2, scale = 0.05) # tail
# TAIL = TAIL[TAIL < 1 ] # remove the few outliers
#
# X = c(B1, TAIL)
# hist(X, breaks = seq(0,1,0.01))
# #
# # # Fit
# # fit = dbpmm.fit(X,
# #                 K = 1:2,
# #                 samples = 10,
# #                 tail = c(TRUE, FALSE),
# #                 epsilon = 1e-10,
# #                 maxIter = 2000,
# #                 fit.type = 'MM')
# # plot(fit)
# # fit$tail
#
#
# find_peaks <- function (x, m = 3){
#   shape <- diff(sign(diff(x, na.pad = FALSE)))
#   pks <- sapply(which(shape < 0), FUN = function(i){
#     z <- i - m + 1
#     z <- ifelse(z > 0, z, 1)
#     w <- i + m + 1
#     w <- ifelse(w < length(x), w, length(x))
#     if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
#   })
#   pks <- unlist(pks)
#   pks
# }
#
# myPeaks = function(X, K, pareto.shape = list(min.val = 0.01, max.val = 5)) {
#
#   # Compute KDE
#   h = hist(X, breaks = seq(0, 1, 0.01), freq = F)
#   maxD = max(h$density)
#
#   # Detect peaks
#   peaks = find_peaks(h$density, 1)
#   x.peaks = (peaks * 0.01)
#
#   # store only peaks above 0.1
#   peaks = peaks[x.peaks > 0.1]
#   x.peaks = x.peaks[x.peaks > 0.1]
#   peakValues = h$density[peaks]
#
#   # Cluster their x-coordinates
#   clus = kmeans(x.peaks, K, nstart = 100)
#
#   # centers = clus$centers[clus$cluster, ]
#   # distances =  sqrt((peaks - centers)^2)
#   #
#   # outliers = order(distances, decreasing=T)[1:5]
#   # # print(outliers)
#   # print(clus$centers)
#   # points(x.peaks[outliers], peakValues[outliers], pch = 8, col = 'black')
#
#   # Colors to spot them
#   col = RColorBrewer::brewer.pal(8, 'Set2')
#
#   # Print the peaks, coloured by cluster assignment
#   df = data.frame(peaks = x.peaks, peakValues, c = clus$cluster)
#   lapply(split(df, f = df$c), function(w)
#          points(w$peaks, w$peakValues, col = col[w$c], pch = 19))
#
#   # A line for the mean
#   abline(v = as.vector(clus$centers), lty = 2, col = col)
#
#   # Rectangles
#   col = ggplot2::alpha(col, 0.5)
#   for(i in 1:K){
#     d = split(df, f = df$c)[[i]]
#     x.left = min(d$peaks)
#     rect(
#       min(d$peaks), 0,
#       max(d$peaks),
#       maxD,
#       col = col[i], density = 23, border = NA)
#   }
#
#   # Random Beta
#   # Sampler for Beta parameters -- reject samples where the mean/variance lead to negative
#   # parameters (a and b are required to be strictly positive)
#   bsamp = function(m) {
#     repeat{
#       v = runif(n = 1, min = 0.001, max = 0.25)
#       p = c(dbpmm:::.estBetaParams(m, v), mean=m, var=v)
#       if(all(p > 0)) return(p)
#     }
#   }
#
#
#   domain = seq(0, 1, 0.01)
#   for(i in 1:K) {
#     par = bsamp(as.vector(clus$centers)[i])
#     lines(domain, dbeta(domain, par$alpha, par$beta), col = col[i])
#   }
#
#   shape = runif(1, min = pareto.shape$min.val, max = pareto.shape$max.val)
#   scale = min(X) - 1e-9
#   lines(domain, sads::dpareto(domain, shape, scale), col = 'red')
#
# }
#
# timon <- readRDS("~/Documents/GitHub/dbpmm/R/simulation_results_theide_blinded.rds")
# myPeaks = Vectorize(myPeaks, vectorize.args = 'K')
#
# w= names(timon)[4]
#
# lapply(names(timon),
#        function(w){
#           X = timon[[w]]$data
#           X = X$alt/X$depth
#
#           pdf(paste(w, '.pdf', sep = ''))
#           myPeaks(X, 1:5)
#           dev.off()
# })
#
# myPeaks(X, 3)
# rect(0.2, 0, .4, 15, col = 'red', density = 23, border = NA)
#
# clus$centers
#
#
#
#
# # m = tapply(distances, clus$cluster,mean)
# # divide each distance by the mean for its cluster:
# # d <- distances/(m[clus$cluster])
# # d = d[order(d, decreasing=TRUE)][1:5]
# # abline(h = d, lty = 2, col = 'red')
#

x=fit$best

dom = seq(0, 1, 0.01)
points = mobster:::template_density(x)
points_y = sapply(points, function(x) x$y)
colnames(points_y) = names(points)

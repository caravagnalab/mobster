BetaParetoMM.mselection = function(X, K, init = 'random', tail = TRUE, epsilon = 1e-6, maxIter = 10000, is_verbose = FALSE, fit.type = 'MM', restarts = 10, parallel = FALSE, cores.ratio =.8, file = NA)
{
  SCORES = NULL
  FIT = NULL

  if(!is.na(file)) pdf(file, width = 8, height = 8)

  # cl = NULL
  # if(parallel) cl = setup_parallel(cores.ratio = cores.ratio)

  for(k in K)
  {
    fit = BetaParetoMM.fit(X = X, K = k, init = init, tail = tail, epsilon = epsilon, maxIter = maxIter, fit.type = fit.type, is_verbose = is_verbose, parallel = parallel, cores.ratio = cores.ratio, restarts = restarts)
    SCORES = append(SCORES, list(fit$scores))
    FIT = append(FIT, list(fit))

    names(SCORES)[length(SCORES)] =
      names(FIT)[length(FIT)] = paste('K =', k)

    if(!is.na(file)) plot(fit, cex = 1)
  }

  # if(parallel) stop_parallel(cl)

  SCORES = Reduce(rbind, SCORES)
  rownames(SCORES) = paste('K =', K)

  cat(bgRed('\n\nBEST MODEL'), red('\n--------------------------------------------------------------\n'))
  print(SCORES)
  cat('\n')
  for(i in 1:ncol(SCORES))
    cat(bgRed(colnames(SCORES)[i]), rownames(SCORES)[which.min(SCORES[,i])], '\n')
  cat(red('--------------------------------------------------------------\n'))

  if(!is.na(file))
  {
    cur = par()$mfrow
  par(mfrow = c(3,2))
  for(i in 1:ncol(SCORES)) {
    plot(K, SCORES[,i], type = 'l', lty = 2, xlab = 'K', ylab = '', col = 'gray')
    points(K, SCORES[,i], pch = 19, col = 'black')
    m = which.min(SCORES[,i])
    points(K[m], SCORES[m,i], pch = 19, col = 'red')

    title(colnames(SCORES)[i])
  }
  par(mfrow = cur)

 dev.off()
  }


  return(list(SCORES = SCORES, FIT = FIT))
}


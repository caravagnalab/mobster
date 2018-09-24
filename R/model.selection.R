# BetaParetoMM.mselection = function(X, K, init = 'random', tail = TRUE, epsilon = 1e-6, maxIter = 10000, is_verbose = FALSE, fit.type = 'MM', restarts = 10, parallel = FALSE, cores.ratio =.8, file = NA)
# {
#   SCORES = NULL
#   FIT = NULL
#
#   if(!is.na(file)) pdf(file, width = 8, height = 8)
#
#   # cl = NULL
#   # if(parallel) cl = setup_parallel(cores.ratio = cores.ratio)
#
#   for(k in K)
#   {
#     fit = BetaParetoMM.fit(X = X, K = k, init = init, tail = tail, epsilon = epsilon, maxIter = maxIter, fit.type = fit.type, is_verbose = is_verbose, parallel = parallel, cores.ratio = cores.ratio, restarts = restarts)
#     SCORES = append(SCORES, list(fit$scores))
#     FIT = append(FIT, list(fit))
#
#     names(SCORES)[length(SCORES)] =
#       names(FIT)[length(FIT)] = paste('K =', k)
#
#     if(!is.na(file)) plot(fit, cex = 1)
#   }
#
#   # if(parallel) stop_parallel(cl)
#
#   SCORES = Reduce(rbind, SCORES)
#   rownames(SCORES) = paste('K =', K)
#
#   cat(bgRed('\n\nBEST MODEL'), red('\n--------------------------------------------------------------\n'))
#   print(SCORES)
#   cat('\n')
#   for(i in 1:ncol(SCORES))
#     cat(bgRed(colnames(SCORES)[i]), rownames(SCORES)[which.min(SCORES[,i])], '\n')
#   cat(red('--------------------------------------------------------------\n'))
#
#   if(!is.na(file))
#   {
#     cur = par()$mfrow
#   par(mfrow = c(3,2))
#   for(i in 1:ncol(SCORES)) {
#     plot(K, SCORES[,i], type = 'l', lty = 2, xlab = 'K', ylab = '', col = 'gray')
#     points(K, SCORES[,i], pch = 19, col = 'black')
#     m = which.min(SCORES[,i])
#     points(K[m], SCORES[m,i], pch = 19, col = 'red')
#
#     title(colnames(SCORES)[i])
#   }
#   par(mfrow = cur)
#
#  dev.off()
#   }
#
#
#   return(list(SCORES = SCORES, FIT = FIT))
# }

#' Title
#'
#' @param results
#' @param scores.suitable
#'
#' @return
#' @export
#'
#' @examples
model_selection = function(results, scores.suitable = c('ICL', 'BIC', 'AIC', 'NLL', 'reICL'), silent = TRUE)
{
  ####################################################################################################
  # New temporary code to update the implementation -- support for "reduced Enrtropy ICL"
  # move this code in the EM function (main fit)
  reduced_entropy = function(z_nk, labels)
  {
    # z_nk = mb$runs[[case]]$z_nk
    # labels = mb$runs[[case]]$labels

    # Cannonical entropy
    entropy = -sum(z_nk * log(z_nk), na.rm = TRUE)

    # Here we take the MAP estimates of z_nk, and select only entries that are assigned
    # to a Beta component (i.e. those with arg_max != Tail)
    cz_nk = z_nk[labels != 'Tail', 2:ncol(z_nk), drop = FALSE]

    # This is un-normalized -- we compute the empirical normalizing constant (C)
    C = rowSums(cz_nk)
    for(i in 1:nrow(cz_nk)) cz_nk[i, ] = cz_nk[i, ]/C[i]

    # The reduced entropy is the entropy of this distribution
    rentropy = -sum(cz_nk * log(cz_nk), na.rm = TRUE)

    # cat("Entropy ", entropy, " [ Reduced ", rentropy,']\n')
    rentropy
  }

  for(i in 1:length(results)) {
    z_nk = results[[i]]$z_nk
    labels = results[[i]]$data$cluster

    rentropy = reduced_entropy(z_nk, labels)
    results[[i]]$scores$reICL = results[[i]]$scores$BIC + rentropy
  }
  ####################################################################################################


  # Exatrct table with scores
  tab = lapply(results, function(w) {
    val = w$scores
    val$K = w$Kbeta
    val$tail = w$fit.tail
    val
  })

  tab = Reduce(rbind, tab)

  model.ids = c("K", 'tail')

  model.selection = NULL
  model.rank = NULL

  for(s in scores.suitable) {


    new.tab = tab[order(tab[, s], decreasing = FALSE), ]
    new.runs = results[as.integer(rownames(new.tab))]
    new.best = new.runs[[1]]

    reorder.fit = list(fits.table = new.tab, runs = new.runs, best = new.best)
    model.selection = append(model.selection, list(reorder.fit))

    model.rank = rbind(model.rank, new.tab[1, model.ids])
  }
  names(model.selection) = scores.suitable
  rownames(model.rank) = scores.suitable

  model.rank$score = rownames(model.rank)

  if(!silent)
  {
    pio::pioTit("Model selection table")
    pio::pioDisp(model.rank)
  }
  
  return(list(model.rank = model.rank, model.selection = model.selection))
}




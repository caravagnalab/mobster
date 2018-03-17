
# BetaParetoMM.comparator = function(X, K, init, tail, epsilon, maxIter, is_verbose, fit.type, restarts, parallel, cores.ratio, file)
# {
#   BPMM = BetaParetoMM.mselection(X = X, K = K, init = init, tail = TRUE,
#                                  epsilon = epsilon, maxIter = maxIter, fit.type = fit.type,
#                                  file = paste(file, 'BPMM-fit', '.pdf', sep = ''),
#                                  restarts = restarts, cores.ratio = cores.ratio, parallel = parallel)
#   save(BPMM, file = paste(file, '-BPMM-fit', '.RData', sep = ''))
#
#   BMM = BetaParetoMM.mselection(X = X, K = K, init = init, tail = FALSE,
#                                 epsilon = epsilon, maxIter = maxIter, fit.type = fit.type,
#                                 file = paste(file, 'BMM-fit', '.pdf', sep = ''),
#                                 restarts = restarts, cores.ratio = cores.ratio, parallel = parallel)
#   save(BPMM, file = paste(file, '-BMM-fit', '.RData', sep = ''))
#
#   cat(bgRed('BPMM\n'))
#   print(BPMM$SCORES)
#
#   cat(bgRed('BMM\n'))
#   print(BMM$SCORES)
#
#   mBPMM = min(BPMM$SCORES[, 'ICL'])
#   mBMM = min(BMM$SCORES[, 'ICL'])
#
#   best = BPMM$FIT
#   scores = BPMM$SCORES
#   best.minScore = mBPMM
#
#   if(mBPMM > mBMM) {
#     best = BMM$FIT
#     scores = BMM$SCORES
#     best.minScore = mBMM
#   }
#
#   best = best[[ rownames(scores)[which(scores[, 'ICL'] == best.minScore)] ]]
#
#   pdf(paste(file, '-best-fit.pdf', sep = ''))
#   plot(best)
#   dev.off()
#
# }



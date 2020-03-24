model_selection = function(x, scores.suitable = c('ICL', 'BIC', 'AIC', 'NLL', 'reICL'), silent = TRUE)
{
  results = x
  
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




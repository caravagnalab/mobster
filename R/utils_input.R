is_list_mobster_fits = function(x)
{
  nOK = all(c('fits.table', 'runs', 'best', 'model.selection') %in% names(x))
  lOK = is.list(x)
  mF = all(sapply(x$runs, class) == 'dbpmm')
  mbF = (class(x$best) == 'dbpmm')
  
  all(nOK, lOK, mF, mbF)
}



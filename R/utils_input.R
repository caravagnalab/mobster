is_mobster_fit = function(x)
{
  if (!inherits(x, 'dbpmm'))
  {
    cat(
      crayon::red(
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nThe input object is not a MOBSTER fit (see below).\n"
      )
    )
    
    cat(
      crayon::red(
        "If you have run \"mobster_fit\" to fit your data, try to access the object named \"best\" (e.g., x$best).\n"
      )
    )
    
    cat(
      crayon::red(
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      )
    )
    
    cat("\nstr() for the input object.\n")
    print(str(x))
    cat("\n\n")
    
    stop("Wrong parameter, a fit was instead expected.")
  }
}

is_mobster_input_dataset = function(x)
{
  merr = function(s) {
    cat(
      crayon::red(
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
        "The input object is not a MOBSTER valid input.\n",
        s, '\n',
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
      )
    )
    
    stop("Bad MOBSTER input.")
  }
  
  if (!is.data.frame(x)) merr("The input is not of class data.frame.")
  if (ncol(x) == 0) merr("The input does not have any column.")
  if (!('VAF' %in% colnames(x))) merr("The input does not have a column named VAF.")
  if (!(all(is.numeric(x$VAF)))) merr("The column named VAF is not numeric.")
  if (any(is.na(x$VAF)) ) merr("The column named VAF contains NA values.")
}



is_list_mobster_fits = function(x)
{
  nOK = all(c('fits.table', 'runs', 'best', 'model.selection') %in% names(x))
  lOK = is.list(x)
  mF = all(sapply(x$runs, class) == 'dbpmm')
  mbF = (class(x$best) == 'dbpmm')
  
  all(nOK, lOK, mF, mbF)
}

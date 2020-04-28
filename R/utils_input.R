is_mobster_fit = function(x)
{
  if (!inherits(x, 'dbpmm'))
  {
    cli::cli_alert_danger(
      "The input object is not a MOBSTER fit; if you have run \"mobster_fit\" access the object named \"best\" (e.g., x$best)."
    )
    
    stop("Wrong parameter, a fit was instead expected.")
  }
}

is_mobster_input_dataset = function(x)
{
  merr = function(s) {
    cli::cli_alert_danger(
      "The input object is not a MOBSTER fit; {.value {s}}."
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
  merr = function(s) {
      cli::cli_alert_danger(
        "The input object is not a list of MOBSTER fits; {.value {s}}."
      )
      
    stop("Bad MOBSTER input (list of fits).")
  }
  
  if(!is.list(x)) merr("Input isn ot a list.")
  if(!all(c('fits.table', 'runs', 'best', 'model.selection') %in% names(x))) merr("Missing names, there should be 'fits.table', 'runs', 'best', 'model.selection' at least.")
  if(!all(sapply(x$runs, class) == 'dbpmm')) merr("Runs objects are not MOBSTER fits.")
  if(class(x$best) != 'dbpmm') merr("Best is not MOBSTER a fits.")
}


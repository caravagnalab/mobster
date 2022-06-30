m_ok = function(m)
{
  paste("{crayon::green(clisymbols::symbol$tick)}", m)
}

m_wrn = function(m)
{
  paste("{crayon::yellow('!')}", m)
}

m_err = function(m)
{
  paste("{crayon::red(clisymbols::symbol$cross)}", m)
}

m_inf = function(m)
{
  paste("{crayon::white(clisymbols::symbol$info)}", m)
}

m_txt = function(m, symbol = 'clisymbols::symbol$pointer')
{
  paste('{', symbol, '}', m)
}

crash_ifnotinstalled = function(packages)
{
  list_installed = installed.packages()
  found = packages %in% list_installed
  
  missing = packages[which(!found)]
  
  for(m in missing) {
    cli::cli_alert_danger("Package {.field {m}} is not installed, and computation cannot be carried out.")
  }
  
  if(any(!found)) stop("Missing required package.")
}

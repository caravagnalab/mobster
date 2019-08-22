.onLoad <- function(libname, pkgname)
{
  requirements = c('tidyverse', 'pio', 'crayon')
  
  suppressMessages(sapply(requirements, require, character.only = TRUE))
  
  invisible()
}
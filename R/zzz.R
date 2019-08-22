.onLoad <- function(libname, pkgname) 
{
  requirements = c('tidyverse', 'pio', 'crayon')
  
  sapply(requirements, function(x) { suppressMessages(require(x)) })
  
  invisible()
}
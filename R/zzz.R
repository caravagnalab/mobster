.onLoad <- function(libname, pkgname)
{
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  requirements = c('tidyverse', 'pio', 'crayon', 'easypar', 'ggpubr', 'sads')
  
  suppressMessages(sapply(requirements, require, character.only = TRUE))
  
  # =-=-=-=-=-=-
  # Package options
  # =-=-=-=-=-=-
  options(pio.string_fg_colour = crayon::bgYellow$black)
  
  # =-=-=-=-=-=-
  # Header
  # =-=-=-=-=-=-
  
  mobster_welcome_message =  getOption('mobster_welcome_message', default = TRUE)
  
  if(mobster_welcome_message)
  {
    pio::pioHdr('MOBSTER - Model-based clustering in cancer')
    pio::pioStr("Author : ", "Giulio Caravagna <gcaravagn@gmail.com>", suffix = '\n')
    pio::pioStr("GitHub : ", "caravagn/mobster", suffix = '\n')
    
    cat(
      "\n > ",
      crayon::blue("[https://https://caravagn.github.io/evoverse]"), 
      "Use package",
      crayon::green("\"evoverse\""),
      "to implement cancer evolution analyses from multi-sample data, integrating MOBSTER, VIBER and CNAqc packages.\n\n")

    options(mobster_welcome_message = FALSE) 
  }
  
  invisible()
}
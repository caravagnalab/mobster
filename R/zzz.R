.onLoad <- function(libname, pkgname)
{
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  
  requirements = c('tidyverse', 'pio', 'crayon')
  
  suppressMessages(sapply(requirements, require, character.only = TRUE))
  
  # =-=-=-=-=-=-
  # Package options
  # =-=-=-=-=-=-
  options(pio.string_fg_colour = crayon::bgYellow$black)
  
  # =-=-=-=-=-=-
  # Header
  # =-=-=-=-=-=-
  pio::pioHdr('MOBSTER - Model-based clustering in cancer')
  pio::pioStr("Author : ", "Giulio Caravagna <gcaravagn@gmail.com>", suffix = '\n')
  pio::pioStr("GitHub : ", "caravagn/mobster", suffix = '\n')
  
  cat(
    "\n > ",
    crayon::blue("[https://github.com/caravagn/mvMOBSTER]"), 
    "See package",
    crayon::green("\"mvmobster\""),
    "for support with multi-region sequencing analyses.\n\n")
  
  
  
  
  invisible()
}
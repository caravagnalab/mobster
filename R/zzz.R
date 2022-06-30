.onLoad <- function(libname, pkgname)
{
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  # requirements = c('tidyverse', 'pio', 'crayon', 'easypar', 'ggpubr', 'sads')
  # 
  # # dplyr,
  # # tidyr,
  # # sads,
  # # ggplot2,
  # # ggpubr,
  # # cowplot,
  # # cli,
  # # crayon,
  # # easypar,
  # # pio,
  # # ctree
  # 
  # suppressMessages(sapply(requirements, require, character.only = TRUE))
  
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
    # pio::pioHdr('MOBSTER - Model-based clustering in cancer')
    # pio::pioStr("Author : ", "Giulio Caravagna <gcaravagn@gmail.com>", suffix = '\n')
    # pio::pioStr("GitHub : ", "caravagn/mobster", suffix = '\n')
    # pio::pioStr("   WWW : ", "https://caravagn.github.io/mobster/", suffix = '\n')
    # 
    # 
    # cat(
    #   "\n > MOBSTER is part of the", crayon::green("\"evoverse\""), 
    #   crayon::blue("[https://bit.ly/2orn94e]"),
    #   "- a collection of packages to implement Cancer Evolution analyses from cancer sequencing data.\n"
    #   )
    
    pk = 'mobster'
    pk_l = 'Model-based clustering in cancer'
    www = "https://caravagn.github.io/mobster/"
    em = "gcaravagn@gmail.com"
    
    cli::cli_text("{crayon::green(clisymbols::symbol$tick)} Loading {.field {pk}}, {.emph \'{pk_l}\'}. Support : {.url { www}}.")
    
    options(mobster_welcome_message = FALSE) 
  }
  
  invisible()
}
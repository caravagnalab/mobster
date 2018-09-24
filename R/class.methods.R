#' Summaries for an object of class dbpmm is like a print.
#'
#' @param x the obj of class dbpmm
#' @param ...
#'
#' @return see print
#' @export
#'
#' @examples something..
summary.dbpmm = function(x, ...) { print.dbpmm(x, ...) }

#' Print an object of class dbpmm.
#'
#' @param x an object of class dbpmm..
#' @param ...
#'
#' @return nothing.
#' @export
#' @import crayon
#'
#' @examples something..
print.dbpmm = function(x, ...)
{
  cat(crayon::bgYellow(crayon::black("[ MOBSTER ]")),
      'N =', x$N, crayon::cyan("samples with"),
      'K =', x$Kbeta, crayon::cyan("Beta clusters, fit by"),
      crayon::yellow(x$fit.type), crayon::cyan('in'), length(x$all.NLL), crayon::cyan("steps"),
      ifelse(x$status, crayon::green('[CONVERGED]'), crayon::red('[NON CONVERGED]')), "\n")

  # pio::pioTit("   Beta components   ")
  # 
  # print(
  #   x$Clusters %>%
  #   dplyr::select(-init.value) %>%
  #   dplyr::filter(cluster != 'Tail', type =='Mean' | type == 'Variance') %>%
  #   dplyr::mutate(fit.value = formatC(fit.value, format = "e", digits = 2)) %>%
  #   spread(key = type, value = fit.value) 
  # )
  #   
  # ####################### Tail
  # pio::pioTit("   Pareto components   ")
  # 
  # if(x$fit.tail) cat(crayon::red('Off\n'))
  # else
  #   print(
  #     x$Clusters %>%
  #     dplyr::select(-init.value) %>%
  #     dplyr::filter(cluster == 'Tail', type =='Shape' | type == 'Scale') %>%
  #     dplyr::mutate(fit.value = formatC(fit.value, format = "e", digits = 2)) %>%
  #     spread(key = type, value = fit.value) 
  #   )
    
  # ####################### Pi
  # pio::pioTit("   Mixing proportions   ")
  # 
  # print(
  #   x$Clusters %>%
  #     dplyr::select(-init.value) %>%
  #     dplyr::filter(type == 'Mixing proportion') %>%
  #     dplyr::mutate(fit.value = formatC(fit.value, format = "e", digits = 2)) %>%
  #     spread(key = type, value = fit.value) 
  # )
  # 
  # cat(crayon::cyan("Clusters dimension"))
  # 
  # clus.size = table(x$labels)
  # clus.size = clus.size[order(clus.size)]
  # 
  # ####################### Scores
  # pio::pioTit("   Scores   ")
  # 
  # print(clus.size)
  # tibble::as.tibble(x$scores)
  
  
  # Pareto = x$Clusters %>%
  #   dplyr::select(-init.value) %>%
  #   dplyr::filter(cluster == 'Tail', type == 'Shape' | type == 'Scale') %>%
  #   dplyr::mutate(fit.value = formatC(fit.value, digits = 2)) 
  # 
  # Betas = x$Clusters %>%
  #   dplyr::select(-init.value) %>%
  #   dplyr::filter(cluster != 'Tail', type == 'Mean' | type == 'Variance') %>%
  #   dplyr::mutate(fit.value = formatC(fit.value, format = "e", digits = 2)) 
  # 
  # Pi = x$Clusters %>%
  #   dplyr::select(-init.value) %>%
  #   dplyr::filter(type == 'Mixing proportion') %>%
  #   dplyr::mutate(fit.value = formatC(fit.value, digits = 4)) 
  #   dplyr::mutate_all(as.character)
  # 

  # Betas =
  #   x$Clusters %>%
  #     dplyr::select(-init.value) %>%
  #     dplyr::filter(cluster != 'Tail', type == 'Mean' | type == 'Variance' | type == 'Mixing proportion') %>%
  #     dplyr::mutate(fit.value = formatC(fit.value, format = "e", digits = 2)) %>% 
  #     spread(key = type, value = fit.value)   
  #   
  # Betas = Betas[, c('Mean', 'Variance', 'Mixing proportion')]
  # 
  # Pareto =
  #   x$Clusters %>%
  #   dplyr::select(-init.value) %>%
  #   dplyr::filter(cluster == 'Tail', type == 'Shape' | type == 'Scale' | type == 'Mixing proportion') %>%
  #   dplyr::mutate(fit.value = formatC(fit.value, format = "e", digits = 2)) %>% 
  #   spread(key = type, value = fit.value)   
  # 
  # Pareto = Pareto[, c('Scale', 'Shape', 'Mixing proportion')]
  # 
  # ####################### Pi
  # pio::pioTit("  Beta components   ")
  # print(Betas)
  # 
  # pio::pioTit("  Pareto component  ")
  # print(Pareto)
  # 
  # pio::pioTit("  Hard clusterting  ")
  # 
  # clus.size = table(x$data$cluster)
  # clus.size = clus.size[order(clus.size)]
  # print(clus.size)
  # 
  # ####################### Scores
  # pio::pioTit("       Scores      ")
  # 
  # print(tibble::as.tibble(x$scores))
  # 
  Betas =
    x$Clusters %>%
    dplyr::select(-init.value) %>%
    dplyr::filter(cluster != 'Tail', type == 'Mean' | type == 'Variance' | type == 'Mixing proportion') %>%
    dplyr::mutate(fit.value = formatC(fit.value, format = "e", digits = 2))

  Pareto =
    x$Clusters %>%
    dplyr::select(-init.value) %>%
    dplyr::filter(cluster == 'Tail', type == 'Shape' | type == 'Scale' | type == 'Mixing proportion') %>%
    dplyr::mutate(fit.value = formatC(fit.value, format = "e", digits = 2))


  ####################### Pi
  clus.size = table(x$data$cluster)
  clus.size = clus.size[order(clus.size)]

  print(clus.size)
  cat('\n')
  
  cat(crayon::black(crayon::bgYellow("  Components (fit)  \n")))
  print(as.data.frame(dplyr::bind_rows(Betas, Pareto)), row.names = FALSE)
  
  ####################### Scores
  cat(crayon::black(crayon::bgYellow("\n  Scores (model selection)  \n")))
  print(x$scores, row.names = FALSE)
  
  # print(tibble::as.tibble(x$scores))
}

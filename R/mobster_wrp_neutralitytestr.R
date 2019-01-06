#' Title
#'
#' @param x
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
mobster_wrp_neutralitytest = function(x, samples = x$samples, ql = 0.05, qr = 0.95)
{
  if (is.null(x$fit.MOBSTER))
    stop("Missing MOBSTER fits from this object?")
  
  tests = lapply(samples,
                 function(sample)
                 {
                   pio::pioTit(paste0(
                     sample,
                     ": neutralitytest from package neutralitytestr (Williams et al.)"
                   ))
                   
                   this.fit = x$fit.MOBSTER[[sample]]$best
                   tail.size = this.fit$Clusters %>%
                     filter(cluster == 'Tail', type == 'Mixing proportion') %>%
                     pull(fit.value)
                   
                   if (tail.size == 0) {
                     message("\nNo tail detected, cannot run test.")
                     
                     return(NULL)
                   }
                   
                   points = this.fit$data %>%
                     filter(cluster == 'Tail') %>% pull(VAF)
                   
                   # Integration range
                   # fmin = round(min(points), 4)
                   # fmax = round(max(points), 4)
                   # 5 an 95% quantiles..
                   fmin = round(quantile(points, ql), 2)
                   fmax = round(quantile(points, qr), 2)  
                   
                   lsq = neutralitytestr::neutralitytest(points, fmin = fmin, fmax = fmax)
                   lsq$N_tail = length(points)
                   
                   pio::pioStr("Lower quantile =", ql)
                   pio::pioStr("Upper quantile =", qr)
                   
                   pio::pioStr("Integration range =",
                               paste0('[', fmin, ', ', fmax, ']\n'))
                   
                   pio::pioStr("Tail size = ",
                               paste0(length(points), '\n'))
                   
                   print(lsq)
                   
                   return(lsq)
                 })
  names(tests) = samples
  
  # Mutation rate from non-null tests
  pio::pioTit(paste0("Mutation rate (linear combination)"))
  
  tests.null =  sapply(tests, is.null)
  
  mr = sapply(tests[!tests.null], function(w) w$mutation.rate)
  
  # Normalized tail size gives the weight of the estimate
  coeff = sapply(tests[!tests.null], function(w) w$N_tail)
  coeff = coeff/sum(coeff)
  
  mutation.rate = mr %*% coeff
  
  pio::pioStr("Mutation rate per sample.", '', suffix = '\n')
  pio::pioDisp(mr)
  
  pio::pioStr("Weight (normalised tail size).", '', suffix = '\n')
  pio::pioDisp(coeff)
  
  pio::pioStr("Linear combination = ", mutation.rate, suffix = '\n')
  
  # Ploits
  
  pio::pioTit(paste0("Computing plots"))
  
  plots = lapply(tests,
                 function(w)
                 {
                   if (is.null(w))
                     return(ggplot() + labs(title = paste0(w, ": No tail detected.")))
                   else
                     ggpubr::annotate_figure(
                       neutralitytestr::plot_all(w),
                       top = paste0(
                         "Neutralitytest from package neutralitytestr (Williams et al.)"
                       )
                     )
                 })
  
  figure = ggarrange(
    plotlist = plots,
    nrow = length(samples),
    ncol = 1,
    labels = samples
  )
  
  return(list(neutralitytest = tests, plot = figure))
}

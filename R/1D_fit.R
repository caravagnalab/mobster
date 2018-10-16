


#' MOBSTER fitting function
#'
#' @param X Data, a tibble with at least a VAF column  with values in (0,1).
#' @param K A vector with the number of Beta components to use. A single  Pareto component is used to model
#' tails, see also "tail". All values of K must be positive and strictly greater than 0.
#' @param n Number of fits that should be attempted for each value of K. This value >1 makes senso only if
#' one uses randomized initial conditions, of course.
#' @param init Initial values for the paremeters. By using "ranodm" mean and variance for each Beta component
#' are randomply sampled in the interval (0,1). PARETO??? If this is a list with K-dimensional vectors named 'a' and 'b',
#' and scalar 'shape' and 'scale' are provided, these parameters will be used.
#' @param tail Set it FALSE to fit the data without the Pareto tail.
#' @param epsilon Tolerance for convergency estimation. For MLE fit this is compared to the differential of the NLL.
#' for MM fit to the largest differential among the mixing proportions (pi).
#' @param maxIter Maximum number of steps of the fit. If convergency is not achieved before these steps, the fit is interrupted.
#' @param is_verbose Verbose output. This also shows a sort of progress bar for the fitting.
#' @param fit.type A string that determines the type of fit. "MLE" is the Maximum Likelihood Estimate of the Beta paraneters,
#' while "MM" is Momemnt Matching. MLE is numerical, and thus slower.
#' @param parallel If TRUE, package parallel is used to run all fits in parallel.
#' @param cores.ratio Ratio (0,1) of the total number of cores that should be used if parallel = TRUE
#' @param file.dump If this is not NA then it should be a string. The method will save all fits in an RData named accordingly
#' and will plot all the fits in a PDF file.
#' @param seed Seed for the random numbers generator
#' @param top Number of top fits to return, ranked by ICL.
#' @param annotation Subtitle annotation, if plotting.
#' @param model.selection Criterion to pick the best model -- one of ICL, BIC, AIC, NLL or entropy. We advise to use only ICL and BIC.
#'
#' @return An object of class "dbpmm" which has methods for print and visualization.
#' @export
#'
#' @import crayon
#' @import parallel
#' @import doParallel
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import magrittr
#'
#' @examples will make some
mobster_fit = function(x,
                       K = 1:3,
                       n = 10,
                       init = 'peaks',
                       tail = c(TRUE, FALSE),
                       epsilon = 1e-10,
                       maxIter = 2000,
                       is_verbose = FALSE,
                       fit.type = 'MM',
                       parallel = FALSE,
                       cores.ratio = .8,
                       file.dump = NA,
                       seed = 12345,
                       annotation = NULL,
                       model.selection = 'reICL',
                       trace = FALSE)
{

  samples = n
  
  # Check for basic input requirements
  stopifnot(is.numeric(samples))
  
  stopifnot(is.data.frame(x))
  stopifnot(ncol(x) > 0) 
  stopifnot('VAF' %in% colnames(x))
   
  X = tibble::as.tibble(x) 

  ###################### Initializations
  set.seed(seed)
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  
  # Storage variables
  best = obj = runs = NULL
  
  # Configurations that will be used for model selection 
  tests = expand.grid(K, 1:samples, tail, stringsAsFactors = FALSE)

  colnames(tests) = c('K', 'Run', 'tail')
  tests = tests[order(tests$tail, tests$K), ]
  
  ntests = nrow(tests)
  
  ###################### Print headers
  pio::pioHdr("MOBSTER fit")

  cat(
    cyan('\n\tDump:'),
    blue(file.dump),
    cyan('\t Annotation:'),
    blue(annotation),
    '\n'
  )



  cat(
    cyan('\n\t-'),
    'N =',
    nrow(x),
    cyan('samples with'),
    'K =',
    K,
    cyan('Beta; Pareto Type-I power-law tail:'),
    ifelse(tail, green('ON'), red('OFF'))
  )
  cat(
    cyan('\n\t- Fit   : '),
    yellow(fit.type),
    cyan('for'),
    maxIter,
    cyan('max. steps with'),
    ifelse(all(is.character(init)), init, 'custom'),
    cyan('initialization;'),
    yellow('\u03B5 ='),
    epsilon,
    cyan('. Model selection: '),
    yellow(model.selection)
  )
  cat(
    cyan('\n\t- Runs  : '),
    samples,
    ' x ',
    length(K),
    ' x ',
    length(tail),
    '=',
    yellow(ntests),
    cyan(' with '),
    ifelse(
      parallel,
      paste(green('PARALLEL'), '[', cores.ratio * 100, '% of cores ]'),
      red('SERIAL')
    ),
    cyan(' output is '),
    ifelse(is_verbose, green('VERBOSE'), red('SILENT')),
    '\n'
  )

  # Input Checks
  if (any(x$VAF == 0)) {
    cat(crayon::red('\n[VAFs = 0] setting them to 1e-9 to avoid numerical errors'))
    
    X$VAF[X$VAF == 0] = 1e-9
  }

  if (any(X == 1)) {
    cat(crayon::red('\n[VAFs = 1] setting them to 1-1e-9 to avoid numerical errors'))
    X$VAF[X$VAF == 1] = 1 - 1e-9
  }
  # END: Input Checks
  
  if (!parallel)
  {
    # formatter
    fmt = function(x, n){ sprintf(paste0('%-', n,'s'), x) }
    fmt = Vectorize(fmt, vectorize.args = 'x')
    
    # Header
    tkS = 15
    
    pio::pioTit(
      paste(fmt(c(paste0('| Run (tot. ', ntests, ')'), '| K ', '| Tail', "| Best", "| Time (m)"), tkS), collapse = '')
    )
    
    for (r in 1:ntests)
    {
      LOCAL.TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

      flush.console()

      # Print values
      s.run = cyan(fmt(paste0('| ', r), tkS))

      s.mod = paste0('| ', tests[r, 'K'])
      s.mod = yellow(fmt(s.mod, tkS))
      
      s.tail = ifelse(tests[r, 'tail'], 
                      green(fmt("| YES", tkS)), 
                      red(fmt("| NO", tkS)))

      cat(paste0('\n', s.run, s.mod, s.tail))
      
      # Compute fit
      obj = .dbpmm.EM(
        x,
        K = tests[r, 'K'],
        init = init,
        tail = tests[r, 'tail'],
        epsilon = epsilon,
        maxIter = maxIter,
        is_verbose = is_verbose,
        fit.type = fit.type,
        trace = trace
      )
      
      runs[[r]] = obj

      # Best score
      s.best = paste0("| ", round(obj$scores[, model.selection], 2))
      s.best = fmt(s.best, tkS)
      
      if (is.null(best) ||
          (!is.null(best) &&
           obj$scores[, model.selection] < best$scores[, model.selection]))
      {
        best = obj

        s.best = bgGreen(s.best)      
      }
      else s.best = white(s.best)   
      

      
      # timing
      LOCAL.TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"),
                            LOCAL.TIME,
                            units = "mins")
      
      s.LOCAL.TIME = as.numeric(LOCAL.TIME, "mins")
      if(s.LOCAL.TIME < 1) s.LOCAL.TIME = "< 1"
      else s.LOCAL.TIME = round(s.LOCAL.TIME, 2)


      s.time = fmt(paste0('| ', sprintf("%-4s", s.LOCAL.TIME)), tkS)
      s.time = blue(s.time)
      
      cat(paste0(s.best, s.time))
    }

  }
  else
  {
    # Setup clusters for parallel computing
    cl = .setup_parallel(cores.ratio = cores.ratio)

    # perform parallel inferences
    r = foreach(
      num = 1:ntests,
      .packages = "crayon",
      .export = ls(globalenv(), all.names = TRUE)
    ) %dopar%
    {
      obj = .dbpmm.EM(
        x,
        K = tests[num, 'K'],
        init = init,
        tail = tests[num, 'tail'],
        epsilon = epsilon,
        maxIter = maxIter,
        is_verbose = is_verbose,
        fit.type = fit.type,
        trace = trace
      )
    }

    runs = r
    .stop_parallel(cl)

    # Get best output a posteriori
    best = NULL
    for (i in 1:ntests) {
      obj = runs[[i]]

      if (is.null(best) ||
          (!is.null(best) &&
           obj$scores[, model.selection] < best$scores[, model.selection]))
        best = obj
    }
  }

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  cat(bold("\n\nCOMPLETED."), round(TIME, 2), cyan('mins'), '\n')

  #### SUBSET TOP FITS
  cat("\n", bold("TOP SCORES:", model.selection, "\n\n"))

  # Get all scores
  tests = cbind(tests, Reduce(rbind, lapply(runs, function(w)
    w$scores)))

  # clean up some repeated results -- show unique fits
  tests$Run = NULL
  scores.columns = colnames(runs[[1]]$scores)
  tests[, scores.columns] = apply(tests[, scores.columns], 2, round, digits = 2)

  # print(tibble::as.tibble(tests))
  # print(!duplicated(tests))

  runs = runs[!duplicated(tests)] # remove duplicated entries..
  tests = tests[!duplicated(tests), ] # remove duplicated entries..
  rownames(tests) = NULL

  # Model selection -- this will be returned later ..
  model = model_selection(runs, scores.suitable = model.selection)
  model = model$model.selection[[model.selection]]
  model$model.selection = model.selection

  print(tibble::as.tibble(model$model.rank))
  
  # at most store 'top' model fits
  # ntests = nrow(tests)
  # if(ntests > top)
  # {
  #   tests = tests[1:top, ]
  #   runs = runs[1:top]
  # }
  # print(tests)

  ###### SHOW BEST FIT
  cat(bold("\n BEST FIT\n\n"))
  print.dbpmm(model$best)

  ###### DUMP DATA TO DISK AS REQUIRED
  if (!is.na(file.dump))
  {
    cat(bold("\n DUMP\n\n"))

    save(model, file = paste(file.dump, '-allRuns.RData', sep = ''))
    cat(yellow('\n- Results : '),
        paste(file.dump, '-allRuns.RData', sep = ''))

    plot_report_MOBSTER(
      model,
      title = 'MOBSTER top fit',
      cex = 1,
      TOP = 5,
      palette = 'Set1',
      boxplot.ICL.range = NULL
    )
    ggsave(paste0(file.dump, "-MOBSTER_fit.pdf"),
           width = 10,
           height = 10)
    
    cat(yellow('\n   - Fits : '),
        paste(file.dump, '-MOBSTER_fit.pdf', sep = ''))
  }

  return(model)
}



.dbpmm.EM <-
  function(X,
           K = 3,
           init = 'peaks',
           tail = TRUE,
           epsilon = 1e-10,
           maxIter = 1000,
           is_verbose = FALSE,
           fit.type = 'MM',
           trace = FALSE) 
  {
    stopifnot(fit.type %in% c('MLE', 'MM'))
    stopifnot(tail | K > 0)

    ##=============================##
    # Create a BetaParetoMM object  #
    ##=============================##
    fit           = list()
    class(fit) <- "dbpmm"

    fit$data      = X
    fit$Call      = match.call()

    fit$fit.type  = fit.type # MLE or MM

    fit$fit.tail  = tail 
    fit$Kbeta     = K
    fit$K         = K + 1 # One extra pareto Mode -- specially assigned to slot #1

    fit$N         = nrow(X)   # Number of samples
    fit$z_nk      = matrix(0, nrow = fit$N, ncol = fit$K)        # Responsibilities
    fit$N.k       = NULL       # Clustering assignments

    fit$pdf.w     = matrix(0, nrow = fit$N, ncol = fit$K)        # Weighted PDFs
    fit$all.NLL   = vector(mode = "numeric")           # Hold NLL for all EM iterations
    fit$NLL       = 1e+40                            # Initialize Negative Log Likelihood

    fit$labels    = NULL # hard clustering assignments
    fit$trace     = NULL # trace for the fit 


    # Names of components 
    names.BetaC = paste('C', 1:K, sep = '')
    names.ParetoC = 'Tail'

    # Compute initial conditions
    fit$Clusters  = .initializer(X$VAF, K = fit$Kbeta, tail = tail, init = init)
    fit$Clusters$fit.value = fit$Clusters$init.value
    
    # Extract Beta values
    fit$a = .params_Beta(fit)$a
    fit$b = .params_Beta(fit)$b
    
    names(fit$a) = names(fit$b) = names.BetaC
    
    # Extract Tail values
    fit$shape = fit$scale = NA
    if(fit$fit.tail)
    {
      fit$shape = .params_Pareto(fit)$Shape
      fit$scale = .params_Pareto(fit)$Scale
    }
    
    # Extract Mixin Prop
    fit$pi = .params_Pi(fit)

    fit$scores = NULL

    logX = log(X$VAF)

    ##=========================================
    # Run Expectation Maximization  algorithm #
    ##=========================================
    for (i in 1:maxIter) {
      ################## 
      # Convergence:
      # - NLL with MLE
      # - pi with MM 
      prevNLL = fit$NLL               
      prevpi = fit$pi
      
      ################## 
      # Store trace to visualize fit
      # current =  data.frame(step = i, NLL = fit$NLL, 
      
      if(trace){
        step.density = fit$Clusters
        step.density$step = i
        fit$trace = dplyr::bind_rows(fit$trace, step.density)
      }

      ##===================
      #       E-Step      #
      ##===================
      # When pi(Pareto) --> 0 the MLE fit for shape --> 0/0 = NaN and thus at some point we get NaN here
      for (k in 1:fit$K)
        fit$pdf.w[, k] = ddbpmm(fit,
                                data = fit$data$VAF,
                                components = k,
                                a = fit$a, 
                                b = fit$b,
                                pi = fit$pi,
                                shape = fit$shape,
                                scale = fit$scale,
                                log = TRUE)

      # Calculate probabilities using the logSumExp trick for numerical stability
      Z          = apply(fit$pdf.w, 1, .log_sum_exp)
      fit$z_nk   = fit$pdf.w - Z
      fit$z_nk   = apply(fit$z_nk, 2, exp)    # Exponentiate to get actual probabilities
      fit$NLL    = -sum(Z)  # Evaluate the NLL
      fit$all.NLL   <-
        c(fit$all.NLL, fit$NLL)    # Keep all NLL in a vector

      if (any(is.infinite(fit$z_nk)))
        stop('Error? All latent variables (z_nk) are Infinite.')


      ##===================
      #       M-Step      #
      ##===================
      N.k   <-
        colSums(fit$z_nk)            # Sum of responsibilities for each cluster
      
      fit$pi  <-
        N.k / fit$N                  # Update mixing proportions for each cluster

      names(fit$pi) = c(names.ParetoC, names.BetaC)
      fit = .set_params_Pi(fit, fit$pi)
      
      # PARETO: NUMERICAL VERSION NOT REQUIRED
      # Pfun = NLLParetoMix(X, z_nk, pi, scale)
      # fit = mle(Pfun, start = list(shape = as.numeric(shape)))
      # shape = coef(fit)['shape']

      # PARETO: analytical MLE
      fit$shape = as.numeric(-1 * (sum(fit$z_nk[, 1])) / (fit$z_nk[, 1] %*% (log(fit$scale) - logX)))
      fit = .set_params_Pareto(fit, fit$shape, fit$scale)
      
      # BETA: numerical MLE or analytical MM
      for (k in 2:fit$K)
      {
        if (fit.type == 'MLE')
          # MLE
        {
          # Compute a functional of the negative logLik, which we minimize
          MLE.fit = stats4::mle(
            minuslogl = .NLLBetaMix(k, X$VAF, fit$z_nk, fit$pi),
            start = list(a = fit$a[k - 1], b = fit$b[k - 1])
          )

          fit$a[k - 1] = as.numeric(stats4::coef(MLE.fit)['a'])
          fit$b[k - 1] = as.numeric(stats4::coef(MLE.fit)['b'])
        }
        else
          # Moments Matching
        {
          mean = as.numeric((fit$z_nk[, k] %*% X$VAF) / (fit$N * fit$pi[k]))
          var = as.numeric((fit$z_nk[, k] %*% ((X$VAF - mean) ** 2)) / (fit$N * fit$pi[k]))

          if (is.na(mean) & is.na(var)) {
            warning('Possible singularity in one Beta component a/b --> Inf.')
          }
          else {
            par = .estBetaParams(mean, var)
            fit$a[k - 1] = par$a
            fit$b[k - 1] = par$b
          }
          
          names(fit$a) = names(fit$b) = names.BetaC
          fit = .set_params_Beta(fit, fit$a, fit$b)
        }
      }

      ## Convergency test
      if (.stoppingCriterion(i,
                             prevNLL,
                             fit$NLL,
                             prevpi,
                             fit$pi,
                             fit.type,
                             epsilon,
                             is_verbose))
        break


    } #End of Expectation Maximization loop.

    if(is_verbose) cat("EM!")
    
    # Status of convergence: TRUE/ FALSE
    fit$status = (i < maxIter)

    ########### Add names to the estimated variables for clarity, populate summary matrices
    names(fit$pi) =  colnames(fit$z_nk) = colnames(fit$pdf.w) = c(names.ParetoC, names.BetaC)
    names(fit$a) = names(fit$b) = names.BetaC
    
    # Update fit table
    fit = .set_params_Beta(fit, fit$a, fit$b)
    fit = .set_params_Pareto(fit, fit$shape, fit$scale)
    fit = .set_params_Pi(fit, fit$pi)
    
    # Cluster labels of each data point. Each data point is assigned to the cluster
    # with the highest posterior responsibility. Hard assignment.
    fit$data$cluster =  unlist(apply(fit$z_nk, 1,
                               function(x) {
                                 names(fit$pi)[which(x == max(x, na.rm = TRUE))[1]]
                               }))

    # Summary numbers
    fit$N.k = rep(0, fit$K)
    names(fit$N.k) = names(fit$pi)
    obFreq = table(fit$data$cluster)

    fit$N.k[names(obFreq)] = obFreq
    
    
    ##==============================
    # Scores for model selection   #
    ##==============================
    fit$scores = latent_vars_scores(
      latent_vars(fit), # Extract latent variables
      fit$K, 
      fit$fit.tail, 
      fit$data$cluster)

    if (is_verbose) {
      print.dbpmm(fit)
      cat('\n')
    }

    return(fit)
  }




#' #' Title
#' #'
#' #' @param X
#' #' @param output.folder
#' #' @param fit.tail
#' #' @param second.fit
#' #' @param do.plots
#' #' @param annotation
#' #' @param dbpmm.fit.params
#' #' @param DP.fit.params
#' #' @param vbdbmm.fit.params
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @import crayon
#' #' @import parallel
#' #' @import doParallel
#' #' @import ggplot2
#' #' @import DPpackage
#' #' @import vbdbmm
#' #' @import BMix
#' #'
#' #' @examples will make some
#' dbpmm.2steps.fit = function(
#'   X,
#'   output.folder = '.',
#'   fit.tail = TRUE,
#'   second.fit = c('DP', 'vbdbmm', 'BMix'),
#'   do.plots = TRUE,
#'   annotation = 'dbpmm.tumour',
#'   dbpmm.fit.params = list(
#'     K = 1,
#'     samples = 1,
#'     init = 'peaks',
#'     tail = TRUE,
#'     epsilon = 1e-10,
#'     maxIter = 6000,
#'     is_verbose = FALSE,
#'     fit.type = 'MM',
#'     parallel = F,
#'     cores.ratio = .8,
#'     file.dump = NA,
#'     seed = 12345,
#'     top = 10,
#'     annotation = NULL
#'   ),
#'   DP.fit.params = list(
#'     alpha_0 = 1e-4,
#'     a1 = 1,
#'     b1 = 1,
#'     ngrid = 100,
#'     nburn = 5000,
#'     nsave = 10000,
#'     nskip = 3,
#'     ndisplay = 100,
#'     state = NULL,
#'     status = TRUE
#'   ),
#'   vbdbmm.fit.params = list(
#'     K = 5,
#'     alpha_0 = 1e-8,
#'     a_0 = 1,
#'     b_0 = 1,
#'     max_iter = NA,
#'     epsilon_conv = 1e-10,
#'     restarts = 10,
#'     parallel = TRUE,
#'     silent = TRUE
#'   ),
#'   bmix.fit.params = list(
#'     K.Binomials = 0:2,
#'     K.BetaBinomials = 0:2,
#'     epsilon = 1e-8,
#'     samples = 10
#'   )
#' )
#' {
#'   if(output.folder != '.') {
#'     current.dir = getwd()
#'     dir.create(output.folder)
#'     setwd(output.folder)
#'   }
#'
#'   if(do.plots) pdf('plot.pdf')
#'
#'
#'   stopifnot(ncol(X) >= 3)
#'   stopifnot(is.data.frame(X))
#'
#'   X = X[, 1:3]
#'   colnames(X) == c("NV", "DP", "VAF")
#'
#'   pio::pioTit("Input data")
#'   pio::pioDisp(X)
#'
#'   # X$VAF = X[, 1]/X[, 2]
#'
#'   fits.table = NULL
#'   fits.tail.table = NULL
#'
#'   # First, if required fit a tail. If you do, then remove tail mutations
#'   best.fit = NULL
#'   if (fit.tail)
#'   {
#'     cat(bgBlue(white('\n\n ===== Tail-detection via Dirichlet Beta-Pareto Mixtures ===== \n\n')))
#'
#'     tail.fit = dbpmm.fit(
#'       X = X$VAF,
#'       K = dbpmm.fit.params$K,
#'       samples = dbpmm.fit.params$samples,
#'       init = dbpmm.fit.params$init,
#'       tail = dbpmm.fit.params$tail,
#'       epsilon = dbpmm.fit.params$epsilon,
#'       maxIter = dbpmm.fit.params$maxIter,
#'       is_verbose = dbpmm.fit.params$is_verbose,
#'       fit.type = dbpmm.fit.params$fit.type,
#'       parallel = dbpmm.fit.params$parallel,
#'       cores.ratio = dbpmm.fit.params$cores.ratio,
#'       file.dump = dbpmm.fit.params$file.dump,
#'       seed = dbpmm.fit.params$seed,
#'       top = dbpmm.fit.params$top,
#'       annotation = NA
#'     )
#'
#'     # Best fit
#'     best.fit = tail.fit[[1]]
#'
#'     fits.tail.table = best.fit$beta[c('a', 'b', 'mean'), , drop = FALSE]
#'
#'     if(do.plots) plot(best.fit, silent = FALSE)
#'
#'     save(best.fit, file = 'dbpmm.fit.RData')
#'
#'     # remove tail mutations for the best fit, if the tail is used
#'     if(all(is.na(best.fit$tail)))
#'     {
#'       cat(red('\n\n ===== Best dbpmm fit does not have a tail ===== \n\n'))
#'       save(X, file = paste('Input-postMOBSTER.RData', sep = ''))
#'     }
#'     else
#'     {
#'       curX = nrow(X)
#'
#'       which.tail = which(best.fit$labels != 'Tail')
#'
#'       X = X[which.tail, ]
#'       save(X, file = paste('Input-postMOBSTER.RData', sep = ''))
#'
#'       cat(green('\n\n ===== Best dbpmm fit with tail, removed',
#'                 (curX-nrow(X)),
#'                 'observations; reduction to',
#'                 round(nrow(X)/curX, 2) * 100,
#'                 '% ===== \n\n'))
#'
#'       stopifnot(ncol(X) > 0)
#'     }
#'   }
#'
#'   fit.vbdbmm = fit.DP = fit.bmix = NULL
#'
#'   #### DP package
#'   if('DP' %in% second.fit)
#'   {
#'     cat(bgBlue(white('\n\n ===== Clustering via Dirichlet Process Binomial Mixtures ===== \n\n')))
#'
#'     # prior -- pointwise or Bayesian
#'     prior = NULL
#'     if(length(DP.fit.params$alpha_0) == 1) prior = list(alpha = DP.fit.params$alpha_0, a1 = DP.fit.params$a1, b1 = DP.fit.params$b1)
#'     if(length(DP.fit.params$alpha_0) == 2) prior = list(a0 = DP.fit.params$alpha_0[1], b0 = DP.fit.params$alpha_0[2], a1 = DP.fit.params$a1, b1 = DP.fit.params$b1)
#'
#'     # fittin algorithm from DPpackage
#'     fit.DP = DPpackage::DPbetabinom(
#'       y = X,
#'       ngrid = DP.fit.params$ngrid,
#'       prior = prior,
#'       mcmc = list(nburn = DP.fit.params$nburn,
#'                   nsave = DP.fit.params$nsave,
#'                   nskip = DP.fit.params$nskip,
#'                   ndisplay = DP.fit.params$ndisplay),
#'       state = DP.fit.params$state,
#'       status = DP.fit.params$status
#'     )
#'
#'     fits.table = append(fits.table, list(.extract.DP.fit(fit.DP)))
#'     names(fits.table)[length(fits.table)] = 'DP'
#'
#'     save(fit.DP, file = paste('DP.fit.RData'))
#'   }
#'
#'   #### vbdbmm package
#'   if('vbdbmm' %in% second.fit)
#'   {
#'     cat(bgBlue(white('\n\n ===== Clustering via Variational Bayes Binomial Mixtures ===== \n\n')))
#'
#'     fit.vbdbmm = vbdbmm::vb_bmm1D_fit(
#'       X = X,
#'       K = vbdbmm.fit.params$K,
#'       alpha_0 = vbdbmm.fit.params$alpha_0,
#'       a_0 = vbdbmm.fit.params$a_0,
#'       b_0 = vbdbmm.fit.params$b_0,
#'       max_iter = vbdbmm.fit.params$max_iter,
#'       epsilon_conv = vbdbmm.fit.params$epsilon_conv,
#'       restarts = vbdbmm.fit.params$restarts,
#'       parallel = vbdbmm.fit.params$parallel,
#'       silent = vbdbmm.fit.params$silent
#'     )
#'
#'     fits.table = append(fits.table, list(.extract.vbdbmm.fit(fit.vbdbmm)))
#'     names(fits.table)[length(fits.table)] = 'vbdbmm'
#'
#'     if(do.plots) vbdbmm::vb_bmm_summary(fit.vbdbmm)
#'
#'     save(fit.vbdbmm, file = paste('vbdbmm.fit.RData'))
#'   }
#'
#'   #### BMix package
#'   if('BMix' %in% second.fit)
#'   {
#'     require(BMix)
#'
#'     cat(bgBlue(white('\n\n ===== Clustering via Dirichlet Process Binomial Mixtures ===== \n\n')))
#'
#'     fit.bmix = BMix::bmixfit(
#'       X,
#'       K.Binomials = bmix.fit.params$K.Binomials,
#'       K.BetaBinomials = bmix.fit.params$K.BetaBinomials,
#'       epsilon = bmix.fit.params$epsilon,
#'       samples = bmix.fit.params$samples
#'     )
#'
#'     fits.table = append(fits.table, list(.extract.Bmix.fit(fit.bmix)))
#'     names(fits.table)[length(fits.table)] = 'Bmix'
#'
#'     if(do.plots) BMix:::plot.bmix(fit.bmix, data = X, coverage = mean(X[, 2]))
#'
#'     save(fit.bmix, file = paste('BMix.fit.RData'))
#'   }
#'
#'   # second.fit = c('DP', 'vbdbmm', 'BMix'),
#'
#'
#'   cat(bgBlue(white('\n\n ===== MOBSTER ===== \n\n')))
#'
#'   cat(bgBlue(white('TAIL')), '\n')
#'   if(!all(is.null(best.fit$tail))) print(data.frame(best.fit$tail, row.names = 'TAIL'))
#'   else print("NONE")
#'
#'   cat(bgBlue(white('\nBeta components')), '\n')
#'   print(fits.tail.table)
#'
#'   cat(bgBlue(white('\nCLONES')), '\n')
#'   fits.table = Reduce(rbind, fits.table)
#'   rownames(fits.table) = second.fit
#'   print(fits.table)
#'
#'   save(fits.table, file = 'fits.table.RData')
#'
#'   if(output.folder != '.') setwd(current.dir)
#'   if(do.plots) dev.off()
#'
#'   fits = NULL
#'   fits$MOBSTER = best.fit
#'   fits$vbdbmm = fit.vbdbmm
#'   fits$DP = fit.DP
#'   fits$BMix = fit.bmix
#'
#'   fits$table = fits.table
#'
#'   fits
#' }

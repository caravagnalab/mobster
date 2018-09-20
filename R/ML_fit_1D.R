
#' dbpmm.fit
#'
#' @param X Data, a vector of observations in (0,1).
#' @param K A vector with the number of Beta components to use. A single  Pareto component is used to model
#' tails, see also "tail". All values of K must be positive and strictly greater than 0.
#' @param samples Number of fits that should be attempted for each value of K. This value >1 makes senso only if
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
#'
#' @examples will make some
dbpmm.fit = function(X,
                     K = 1:3,
                     samples = 10,
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
                     # top = length(K) * samples * length(tail),
                     annotation = NULL,
                     model.selection = 'reICL'
)
{
  best = obj = runs = NULL
  stopifnot(is.numeric(samples))
  set.seed(seed)

  cat(bgYellow(black(" [ MOBSTER ] ")),
      yellow(' Finite Dirichelt Mixture Models with Beta and Pareto mixtures (univariate)\n'))

  cat(cyan('\n\tDump:'), blue(file.dump), cyan('\t Annotation:'), blue(annotation),'\n')


  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  tests = expand.grid(K, 1:samples, tail, stringsAsFactors = FALSE)
  tests = tests[order(tests[, 1]), ]

  colnames(tests) = c('K', 'Run', 'tail')
  ntests = nrow(tests)

  cat(cyan('\n\t-'), 'N =',  length(X), cyan('samples with'), 'K =', K, cyan('Beta; Pareto Type-I power-law tail:'), ifelse(tail, green('ON'), red('OFF')))
  cat(cyan('\n\t- Fit   : '), yellow(fit.type), cyan('for'), maxIter, cyan('max. steps with'),
      ifelse(all(is.character(init)), init, 'custom'), cyan('initialization;'),  yellow('\u03B5 ='), epsilon, cyan('. Model selection: '), yellow(model.selection))
  cat(cyan('\n\t- Runs  : '), samples, ' x ', length(K), ' x ', length(tail), '=', yellow(ntests), cyan(' with '),
      ifelse(parallel, paste(green('PARALLEL'), '[', cores.ratio * 100, '% of cores ]'), red('SERIAL')),
      cyan(' output is '),
      ifelse(is_verbose, green('VERBOSE'), red('SILENT')), '\n'
  )

  # Input Checks
  if(any(X == 0)) {
    cat(crayon::red('\n[VAFs = 0] setting them to 1e-9 to avoid numerical errors'))
    X[X == 0] = 1e-9
  }

  if(any(X == 1)) {
    cat(crayon::red('\n[VAFs = 1] setting them to 1-1e-9 to avoid numerical errors'))
    X[X == 1] = 1-1e-9
  }
  # END: Input Checks

  if(!parallel)
  {
    for(r in 1:ntests)
    {
      LOCAL.TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

      flush.console()

      cat(
        bgBlue(
          '\n**** Run',
          sprintf('%3s /', r), sprintf('%3s', ntests), ' : '),
        yellow(
          sprintf("K = %-10s", paste0(tests[r, 'K'], ifelse(tests[r, 'tail'], " + Tail", "")))))

      if(is_verbose) cat('\n')

      obj = .dbpmm.EM(X, K = tests[r, 'K'], init = init, tail = tests[r, 'tail'], epsilon = epsilon, maxIter = maxIter, is_verbose = is_verbose, fit.type = fit.type)
      runs[[r]] = obj

      if (is.null(best) || (!is.null(best) && obj$scores[, model.selection] < best$scores[, model.selection]))
      {
        best = obj

        cat(green(" New best", model.selection))
      }
      else cat(red('    Worst', model.selection))

      LOCAL.TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), LOCAL.TIME, units = "mins")
      cat(' Min(s): ',
          blue(
            sprintf("%-8s", round(LOCAL.TIME, 2))
          ),
          ';', model.selection, ' = ', obj$scores[, model.selection])
    }

  }
  else
  {
    # Setup clusters for parallel computing
    cl = .setup_parallel(cores.ratio = cores.ratio)

    # perform parallel inferences
    r = foreach(num = 1:ntests, .packages = "crayon", .export = ls(globalenv(), all.names = TRUE)) %dopar%
    {
      obj = .dbpmm.EM(X, K = tests[num, 'K'], init = init, tail = tests[num, 'tail'], epsilon = epsilon, maxIter = maxIter, is_verbose = is_verbose, fit.type = fit.type)
    }

    runs = r
    .stop_parallel(cl)

    # Get best output a posteriori
    best = NULL
    for(i in 1:ntests) {
      obj = runs[[i]]

      if (is.null(best) || (!is.null(best) && obj$scores[, model.selection] < best$scores[, model.selection])) best = obj
    }
  }

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

  cat(bold("\n\nCOMPLETED."), blue(round(TIME, 2)), cyan('mins'), '\n')

  #### SUBSET TOP FITS
  cat(bold("\n", model.selection, "TOP OF ALL DISTINCT FITS\nn"))

  # Get all scores
  tests = cbind(tests, Reduce(rbind, lapply(runs, function(w) w$scores)))

  # clean up some repeated results -- show unique fits
  tests$Run = NULL
  scores.columns = colnames(runs[[1]]$scores)
  tests[, scores.columns] = apply(tests[, scores.columns], 2, round, digits = 2)

  print(tests)
  print(!duplicated(tests))

  runs = runs[!duplicated(tests)] # remove duplicated entries..
  tests = tests[!duplicated(tests), ] # remove duplicated entries..
  rownames(tests) = NULL

  # Model selection -- this will be returned later ..
  model = model_selection(runs, scores.suitable = model.selection)
  model = model$model.selection[[model.selection]]
  model$model.selection = model.selection

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

  if(!is.na(file.dump))
  {
    cat(bold("\n DUMP\n\n"))

    save(model, file = paste(file.dump, '-allRuns.RData', sep = ''))
    cat(yellow('\n- Results : '), paste(file.dump, '-allRuns.RData', sep = ''))

    plot_report_MOBSTER(model,
                        title = 'MOBSTER top fit', cex = 1, TOP = 5,
                        palette = 'Set1', boxplot.ICL.range = NULL)
    ggsave(paste0(file.dump, "-MOBSTER_fit.pdf"), width = 10, height = 10)
    cat(yellow('\n   - Fits : '), paste(file.dump, '-MOBSTER_fit.pdf', sep = ''))

    # pdf(paste(file.dump, '-topFits.pdf', sep = ''), height = 8, width = 13)
    # lapply(runs, plot, annotation = annotation)
    # dev.off()
    # cat(yellow('\n- Top fits : '), paste(file.dump, '-topFits.pdf', sep = ''))
    #
    # pdf(paste(file.dump, '-bestFit.pdf', sep = ''), height = 8, width = 13)
    # plot(best, annotation = annotation)
    # dev.off()
    # cat(yellow('\n- Best fit : '), paste(file.dump, '-bestFit.pdf', sep = ''))
    #
    # pdf(paste(file.dump, '-boxplotFitScores.pdf', sep = ''), height = 5, width = 8)
    # .plot.fit.summary(all.tests)
    # dev.off()
    # cat(yellow('\n-  Boxplot : '), paste(file.dump, '-boxplotFitScores.pdf', sep = ''))
  }


  return(model)
}


.dbpmm.EM <- function(X, K = 3, init = 'peaks', tail = TRUE, epsilon = 1e-10,
                      maxIter = 1000, is_verbose = FALSE, fit.type = 'MM'){

  stopifnot(fit.type %in% c('MLE', 'MM'))
  stopifnot(tail | K > 0)

  ##=============================##
  # Create a BetaParetoMM object  #
  ##=============================##
  fit           = list()
  class(fit) <- "dbpmm"

  fit$X         = X
  fit$Call      = match.call()

  fit$fit.type  = fit.type # MLE or MM
  fit$Kbeta     = K
  fit$K         = K + 1 # One extra pareto Mode -- specially assigned to slot #1

  fit$N         = length(X)   # Number of samples
  fit$z_nk      = matrix(0, nrow = fit$N, ncol = fit$K)        # Responsibilities
  fit$N.k       = NULL       # Clustering assignments

  fit$pdf.w     = matrix(0, nrow = fit$N, ncol = fit$K)        # Weighted PDFs
  fit$all.NLL   = vector(mode="numeric")           # Hold NLL for all EM iterations
  fit$NLL       = 1e+40                            # Initialize Negative Log Likelihood

  fit$labels    = NULL # hard clustering assignments

  ### Compute initial conditions
  names.BetaC = paste('C', 1:K, sep = '')
  names.ParetoC = 'Tail'

  fit$init      = .initializer(X, fit$Kbeta, init)

  fit$a         = unlist(fit$init$beta['a', ])
  fit$b         = unlist(fit$init$beta['b', ])
  fit$shape     = fit$init$shape
  fit$scale     = fit$init$scale

  fit$beta      = rbind(fit$a, fit$b)
  colnames(fit$beta) = names.BetaC
  rownames(fit$beta) = c('a', 'b')

  fit$tail      = NA
  if(tail) fit$tail = list(shape = fit$shape, scale = fit$scale)

  fit$pi        = c(0, rep(1/fit$Kbeta, fit$Kbeta)) # Without tail, pi[1] = 0
  if(tail) fit$pi = rep(1/fit$K, fit$K)
  fit$pi.init = fit$pi

  fit$scores = NULL

  if (is_verbose) {
    print.dbpmm(fit)
    cat('\n')
  }

  logX = log(X)

  ##=========================================
  # Run Expectation Maximization  algorithm #
  ##=========================================
  for (i in 1:maxIter) {              # Loop until convergence
    prevNLL  <- fit$NLL               # Store NLL to check for convergence in MLE fit
    prevpi = fit$pi                   # Store pi to check for convergence in MM fit

    ##===================
    #       E-Step      #
    ##===================
    # When pi(Pareto) --> 0 the MLE fit for shape --> 0/0 = NaN and thus at some point we get NaN here

    # print(fit$shape)

    for (k in 1:fit$K) fit$pdf.w[, k] = ddbpmm(fit, data = X, components = k, log = TRUE)

    # ddbpmm(fit, data = X, components = 4, log = TRUE)

    # if(any(is.infinite(fit$pdf.w))) {
    #   cat(red('KITKAT'))
    #   break;
    # }
    # print(any(is.infinite(fit$pdf.w)))

    # Calculate probabilities using the logSumExp trick for numerical stability
    Z          = apply(fit$pdf.w, 1, .log_sum_exp)
    fit$z_nk   = fit$pdf.w - Z
    fit$z_nk   = apply(fit$z_nk, 2, exp)    # Exponentiate to get actual probabilities
    fit$NLL    = -sum(Z)  # Evaluate the NLL
    fit$all.NLL   <- c(fit$all.NLL, fit$NLL)    # Keep all NLL in a vector

    if(any(is.infinite(fit$z_nk))) stop('Error? All latent variables (z_nk) are Infinite.')


    ##===================
    #       M-Step      #
    ##===================
    N.k   <- colSums(fit$z_nk)            # Sum of responsibilities for each cluster
    fit$pi  <- N.k/fit$N                  # Update mixing proportions for each cluster

    # cat('\n')
    # print(N.k)
    # print(fit$pi)


    # PARETO: NUMERICAL VERSION NOT REQUIRED
    # Pfun = NLLParetoMix(X, z_nk, pi, scale)
    # fit = mle(Pfun, start = list(shape = as.numeric(shape)))
    # shape = coef(fit)['shape']

    # PARETO: analytical MLE
    fit$shape = as.numeric(-1 * (sum(fit$z_nk[, 1])) / (fit$z_nk[, 1] %*% (log(fit$scale) - logX)))

    # BETA: numerical MLE or analytical MM
    for (k in 2:fit$K)
    {
      if(fit.type == 'MLE') # MLE
      {
        # Compute a functional of the negative logLik, which we minimize
        MLE.fit = stats4::mle(
          minuslogl = .NLLBetaMix(k, X, fit$z_nk, fit$pi),
          start = list(a = fit$a[k-1], b = fit$b[k-1]))

        # print(MLE.fit)
        # print(as.numeric(MLE.fit['a']))
        # print(coef(MLE.fit))
        # print(coef(MLE.fit))
        # print(fit)
        # print(as.numeric(coef(MLE.fit)['a']))
        # stop('xxxx')

        fit$a[k-1] = as.numeric(stats4::coef(MLE.fit)['a'])
        fit$b[k-1] = as.numeric(stats4::coef(MLE.fit)['b'])
      }
      else # Moments Matching
      {
        mean = as.numeric((fit$z_nk[, k] %*% X)/(fit$N * fit$pi[k]))
        var = as.numeric((fit$z_nk[, k] %*% ((X - mean)**2))/(fit$N * fit$pi[k]))

        if(is.na(mean) & is.na(var)) {
          warning('Possible singularity in one Beta component a/b --> Inf.')
        }
        else {
          par = .estBetaParams(mean, var)
          fit$a[k-1] = par$a
          fit$b[k-1] = par$b
        }
      }
    }

    ## Convergency test
    if(.stoppingCriterion(i, prevNLL, fit$NLL, prevpi, fit$pi, fit.type, epsilon, is_verbose)) break;

  } #End of Expectation Maximization loop.

  if (is_verbose) cat('\n', bgGreen("STATUS"), ifelse(i == maxIter, red('NOT CONVERGED'), green('CONVERGED')))
  fit$status = (i < maxIter)

  ########### Add names to the estimated variables for clarity, populate summary matrices
  names(fit$pi) =  colnames(fit$z_nk) = colnames(fit$pdf.w) = c(names.ParetoC, names.BetaC)
  names(fit$a) = names(fit$b) = names.BetaC

  fit$beta = rbind(fit$beta, fit$a, fit$b)

  # Means and variances of estimated parameters
  means = NULL
  for(i in 1:K)
    means = cbind(means, unlist(.MeanVarBeta(fit$a[i], fit$b[i])))

  fit$beta = rbind(fit$beta, means)
  rownames(fit$beta) = c('a.init', 'b.init', 'a', 'b', 'mean', 'var')

  if(tail) {
    names(fit$tail) = c('shape.init', 'scale.init')
    fit$tail = append(fit$tail, list(shape = fit$shape, scale = fit$scale))
    fit$tail = append(fit$tail, .MeanVarPareto(fit$shape, fit$scale))
  }

  # Cluster labels of each data point. Each data point is assigned to the cluster
  # with the highest posterior responsibility. Hard assignment.
  fit$labels =  unlist(
    apply(fit$z_nk, 1,
          function(x) {
            names(fit$pi)[which(x == max(x, na.rm = TRUE))[1]]
          }))

  # Summary numbers
  fit$N.k = rep(0, fit$K)
  names(fit$N.k) = names(fit$pi)
  obFreq = table(fit$labels)

  fit$N.k[names(obFreq)] = obFreq

  ##==============================
  # Scores for model selection   #
  ##==============================
  if(tail)
    numParams = fit$K + 2 * fit$Kbeta + 1          # Total number of parameters i.e. pi (Beta + Pareto) + 2 * Kbeta (Beta) + 1 (Pareto)
  else
    numParams = (fit$K - 1) + 2 * fit$Kbeta         # Total number of parameters i.e. pi (Beta) + 2 * Kbeta (Beta)

  BIC <- 2 * fit$NLL + numParams * log(fit$N)     # BIC = -2*ln(L) + params*ln(N)
  AIC <- 2 * fit$NLL + 2 * numParams              # AIC = -2*ln(L) + 2*params


  # Integrated Complete Likelihood criterion -- uses standard entropy
  entropy <- -sum(fit$z_nk * log(fit$z_nk), na.rm = TRUE)
  ICL <- BIC + entropy

  # Integrated Complete Likelihood criterion with reduced entropy (only for
  # latent variable that involve subclones -- i.e., Betas). I think this is also a sort of
  # conditional entropy where we condition on the MAP estimate of a mutation being part
  # of a clone, rather than tail.

  # HTake the MAP estimates of z_nk, and select only entries that are assigned
  # to a Beta component (i.e. those with arg_max != Tail)
  fit$cz_nk = fit$z_nk[fit$labels != 'Tail', 2:ncol(fit$z_nk), drop = FALSE]

  # This is un-normalized -- we compute the empirical normalizing constant (C)
  C = rowSums(fit$cz_nk)
  for(i in 1:nrow(fit$cz_nk)) fit$cz_nk [i, ] = fit$cz_nk [i, ]/C[i]

  # The reduced entropy is the entropy of this distribution
  rentropy = -sum(fit$cz_nk  * log(fit$cz_nk ), na.rm = TRUE)

  # Integrated Complete Likelihood criterion with reduced entropy
  reICL <- BIC + rentropy

  fit$scores = data.frame(NLL = fit$NLL, BIC = BIC,
                          AIC = AIC, entropy = entropy, ICL = ICL,
                          reICL = reICL,
                          size = numParams)

  if (is_verbose) {
    print.dbpmm(fit)
    cat('\n')
  }

  return(fit)
}


# BetaParetoMM.mselection = function(X, K, init = 'random', tail = TRUE, epsilon = 1e-6, maxIter = 10000, is_verbose = FALSE, fit.type = 'MM', restarts = 10, parallel = FALSE, cores.ratio =.8, file = NA)
# {
#   SCORES = NULL
#   FIT = NULL
#
#   if(!is.na(file)) pdf(file, width = 8, height = 8)
#
#   # cl = NULL
#   # if(parallel) cl = setup_parallel(cores.ratio = cores.ratio)
#
#   for(k in K)
#   {
#     fit = BetaParetoMM.fit(X = X, K = k, init = init, tail = tail, epsilon = epsilon, maxIter = maxIter, fit.type = fit.type, is_verbose = is_verbose, parallel = parallel, cores.ratio = cores.ratio, restarts = restarts)
#     SCORES = append(SCORES, list(fit$scores))
#     FIT = append(FIT, list(fit))
#
#     names(SCORES)[length(SCORES)] =
#       names(FIT)[length(FIT)] = paste('K =', k)
#
#     if(!is.na(file)) plot(fit, cex = 1)
#   }
#
#   # if(parallel) stop_parallel(cl)
#
#   SCORES = Reduce(rbind, SCORES)
#   rownames(SCORES) = paste('K =', K)
#
#   cat(bgRed('\n\nBEST MODEL'), red('\n--------------------------------------------------------------\n'))
#   print(SCORES)
#   cat('\n')
#   for(i in 1:ncol(SCORES))
#     cat(bgRed(colnames(SCORES)[i]), rownames(SCORES)[which.min(SCORES[,i])], '\n')
#   cat(red('--------------------------------------------------------------\n'))
#
#   if(!is.na(file))
#   {
#     cur = par()$mfrow
#   par(mfrow = c(3,2))
#   for(i in 1:ncol(SCORES)) {
#     plot(K, SCORES[,i], type = 'l', lty = 2, xlab = 'K', ylab = '', col = 'gray')
#     points(K, SCORES[,i], pch = 19, col = 'black')
#     m = which.min(SCORES[,i])
#     points(K[m], SCORES[m,i], pch = 19, col = 'red')
#
#     title(colnames(SCORES)[i])
#   }
#   par(mfrow = cur)
#
#  dev.off()
#   }
#
#
#   return(list(SCORES = SCORES, FIT = FIT))
# }





#' Title
#'
#' @param X
#' @param output.folder
#' @param fit.tail
#' @param second.fit
#' @param do.plots
#' @param annotation
#' @param dbpmm.fit.params
#' @param DP.fit.params
#' @param vbdbmm.fit.params
#'
#' @return
#' @export
#'
#' @import crayon
#' @import parallel
#' @import doParallel
#' @import ggplot2
#' @import DPpackage
#' @import vbdbmm
#' @import BMix
#'
#' @examples will make some
dbpmm.2steps.fit = function(
  X,
  output.folder = '.',
  fit.tail = TRUE,
  second.fit = c('DP', 'vbdbmm', 'BMix'),
  do.plots = TRUE,
  annotation = 'dbpmm.tumour',
  dbpmm.fit.params = list(
    K = 1,
    samples = 1,
    init = 'peaks',
    tail = TRUE,
    epsilon = 1e-10,
    maxIter = 6000,
    is_verbose = FALSE,
    fit.type = 'MM',
    parallel = F,
    cores.ratio = .8,
    file.dump = NA,
    seed = 12345,
    top = 10,
    annotation = NULL
  ),
  DP.fit.params = list(
    alpha_0 = 1e-4,
    a1 = 1,
    b1 = 1,
    ngrid = 100,
    nburn = 5000,
    nsave = 10000,
    nskip = 3,
    ndisplay = 100,
    state = NULL,
    status = TRUE
  ),
  vbdbmm.fit.params = list(
    K = 5,
    alpha_0 = 1e-8,
    a_0 = 1,
    b_0 = 1,
    max_iter = NA,
    epsilon_conv = 1e-10,
    restarts = 10,
    parallel = TRUE,
    silent = TRUE
  ),
  bmix.fit.params = list(
    K.Binomials = 0:2,
    K.BetaBinomials = 0:2,
    epsilon = 1e-8,
    samples = 10
  )
)
{
  if(output.folder != '.') {
    current.dir = getwd()
    dir.create(output.folder)
    setwd(output.folder)
  }

  if(do.plots) pdf('plot.pdf')


  stopifnot(ncol(X) >= 3)
  stopifnot(is.data.frame(X))

  X = X[, 1:3]
  colnames(X) == c("NV", "DP", "VAF")

  pio::pioTit("Input data")
  pio::pioDisp(X)

  # X$VAF = X[, 1]/X[, 2]

  fits.table = NULL
  fits.tail.table = NULL

  # First, if required fit a tail. If you do, then remove tail mutations
  best.fit = NULL
  if (fit.tail)
  {
    cat(bgBlue(white('\n\n ===== Tail-detection via Dirichlet Beta-Pareto Mixtures ===== \n\n')))

    tail.fit = dbpmm.fit(
      X = X$VAF,
      K = dbpmm.fit.params$K,
      samples = dbpmm.fit.params$samples,
      init = dbpmm.fit.params$init,
      tail = dbpmm.fit.params$tail,
      epsilon = dbpmm.fit.params$epsilon,
      maxIter = dbpmm.fit.params$maxIter,
      is_verbose = dbpmm.fit.params$is_verbose,
      fit.type = dbpmm.fit.params$fit.type,
      parallel = dbpmm.fit.params$parallel,
      cores.ratio = dbpmm.fit.params$cores.ratio,
      file.dump = dbpmm.fit.params$file.dump,
      seed = dbpmm.fit.params$seed,
      top = dbpmm.fit.params$top,
      annotation = NA
    )

    # Best fit
    best.fit = tail.fit[[1]]

    fits.tail.table = best.fit$beta[c('a', 'b', 'mean'), , drop = FALSE]

    if(do.plots) plot(best.fit, silent = FALSE)

    save(best.fit, file = 'dbpmm.fit.RData')

    # remove tail mutations for the best fit, if the tail is used
    if(all(is.na(best.fit$tail)))
    {
      cat(red('\n\n ===== Best dbpmm fit does not have a tail ===== \n\n'))
      save(X, file = paste('Input-postMOBSTER.RData', sep = ''))
    }
    else
    {
      curX = nrow(X)

      which.tail = which(best.fit$labels != 'Tail')

      X = X[which.tail, ]
      save(X, file = paste('Input-postMOBSTER.RData', sep = ''))

      cat(green('\n\n ===== Best dbpmm fit with tail, removed',
                (curX-nrow(X)),
                'observations; reduction to',
                round(nrow(X)/curX, 2) * 100,
                '% ===== \n\n'))

      stopifnot(ncol(X) > 0)
    }
  }

  fit.vbdbmm = fit.DP = fit.bmix = NULL

  #### DP package
  if('DP' %in% second.fit)
  {
    cat(bgBlue(white('\n\n ===== Clustering via Dirichlet Process Binomial Mixtures ===== \n\n')))

    # prior -- pointwise or Bayesian
    prior = NULL
    if(length(DP.fit.params$alpha_0) == 1) prior = list(alpha = DP.fit.params$alpha_0, a1 = DP.fit.params$a1, b1 = DP.fit.params$b1)
    if(length(DP.fit.params$alpha_0) == 2) prior = list(a0 = DP.fit.params$alpha_0[1], b0 = DP.fit.params$alpha_0[2], a1 = DP.fit.params$a1, b1 = DP.fit.params$b1)

    # fittin algorithm from DPpackage
    fit.DP = DPpackage::DPbetabinom(
      y = X,
      ngrid = DP.fit.params$ngrid,
      prior = prior,
      mcmc = list(nburn = DP.fit.params$nburn,
                  nsave = DP.fit.params$nsave,
                  nskip = DP.fit.params$nskip,
                  ndisplay = DP.fit.params$ndisplay),
      state = DP.fit.params$state,
      status = DP.fit.params$status
    )

    fits.table = append(fits.table, list(.extract.DP.fit(fit.DP)))
    names(fits.table)[length(fits.table)] = 'DP'

    save(fit.DP, file = paste('DP.fit.RData'))
  }

  #### vbdbmm package
  if('vbdbmm' %in% second.fit)
  {
    cat(bgBlue(white('\n\n ===== Clustering via Variational Bayes Binomial Mixtures ===== \n\n')))

    fit.vbdbmm = vbdbmm::vb_bmm1D_fit(
      X = X,
      K = vbdbmm.fit.params$K,
      alpha_0 = vbdbmm.fit.params$alpha_0,
      a_0 = vbdbmm.fit.params$a_0,
      b_0 = vbdbmm.fit.params$b_0,
      max_iter = vbdbmm.fit.params$max_iter,
      epsilon_conv = vbdbmm.fit.params$epsilon_conv,
      restarts = vbdbmm.fit.params$restarts,
      parallel = vbdbmm.fit.params$parallel,
      silent = vbdbmm.fit.params$silent
    )

    fits.table = append(fits.table, list(.extract.vbdbmm.fit(fit.vbdbmm)))
    names(fits.table)[length(fits.table)] = 'vbdbmm'

    if(do.plots) vbdbmm::vb_bmm_summary(fit.vbdbmm)

    save(fit.vbdbmm, file = paste('vbdbmm.fit.RData'))
  }

  #### BMix package
  if('BMix' %in% second.fit)
  {
    require(BMix)

    cat(bgBlue(white('\n\n ===== Clustering via Dirichlet Process Binomial Mixtures ===== \n\n')))

    fit.bmix = BMix::bmixfit(
      X,
      K.Binomials = bmix.fit.params$K.Binomials,
      K.BetaBinomials = bmix.fit.params$K.BetaBinomials,
      epsilon = bmix.fit.params$epsilon,
      samples = bmix.fit.params$samples
    )

    fits.table = append(fits.table, list(.extract.Bmix.fit(fit.bmix)))
    names(fits.table)[length(fits.table)] = 'Bmix'

    if(do.plots) BMix:::plot.bmix(fit.bmix, data = X, coverage = mean(X[, 2]))

    save(fit.bmix, file = paste('BMix.fit.RData'))
  }

  # second.fit = c('DP', 'vbdbmm', 'BMix'),


  cat(bgBlue(white('\n\n ===== MOBSTER ===== \n\n')))

  cat(bgBlue(white('TAIL')), '\n')
  if(!all(is.null(best.fit$tail))) print(data.frame(best.fit$tail, row.names = 'TAIL'))
  else print("NONE")

  cat(bgBlue(white('\nBeta components')), '\n')
  print(fits.tail.table)

  cat(bgBlue(white('\nCLONES')), '\n')
  fits.table = Reduce(rbind, fits.table)
  rownames(fits.table) = second.fit
  print(fits.table)

  save(fits.table, file = 'fits.table.RData')

  if(output.folder != '.') setwd(current.dir)
  if(do.plots) dev.off()

  fits = NULL
  fits$MOBSTER = best.fit
  fits$vbdbmm = fit.vbdbmm
  fits$DP = fit.DP
  fits$BMix = fit.bmix

  fits$table = fits.table

  fits
}

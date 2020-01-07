# The actual fitting function
.dbpmm.EM <-
  function(X,
           K = 3,
           init = 'peaks',
           tail = TRUE,
           epsilon = 1e-10,
           maxIter = 1000,
           is_verbose = FALSE,
           fit.type = 'MM',
           trace = FALSE,
           pi_cutoff = 0.02,
           N_cutoff = 10,
           description = "My MOBSTER model.")
  {
    stopifnot(fit.type %in% c('MLE', 'MM'))
    stopifnot(tail | K > 0)
    
    .onLoad(NULL, NULL)
    
    # suppressMessages(require(tidyverse))
    # suppressMessages(require(pio))
    # suppressMessages(require(crayon))

    ##=============================##
    # Create a BetaParetoMM object  #
    ##=============================##
    fit           = list()
    class(fit) <- "dbpmm"

    fit$data      = X
    fit$Call      = match.call()
    fit$description = description

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
    fit$Clusters  = mobster:::.initializer(X$VAF, K = fit$Kbeta, tail = tail, init = init)
    fit$Clusters$fit.value = fit$Clusters$init.value

    # Extract Beta values
    fit$a = mobster:::.params_Beta(fit)$a
    fit$b = mobster:::.params_Beta(fit)$b

    names(fit$a) = names(fit$b) = names.BetaC

    # Extract Tail values
    fit$shape = fit$scale = NA
    if(fit$fit.tail)
    {
      fit$shape = mobster:::.params_Pareto(fit)$Shape
      fit$scale = mobster:::.params_Pareto(fit)$Scale
    }

    # Extract Mixin Prop
    fit$pi = mobster:::.params_Pi(fit)

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
        fit$pdf.w[, k] = mobster::ddbpmm(fit,
                                data = fit$data$VAF,
                                components = k,
                                a = fit$a,
                                b = fit$b,
                                pi = fit$pi,
                                shape = fit$shape,
                                scale = fit$scale,
                                log = TRUE)

      # Calculate probabilities using the logSumExp trick for numerical stability
      Z          = apply(fit$pdf.w, 1, mobster:::.log_sum_exp)
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
      fit = mobster:::.set_params_Pi(fit, fit$pi)

      # PARETO: NUMERICAL VERSION NOT REQUIRED
      # Pfun = NLLParetoMix(X, z_nk, pi, scale)
      # fit = mle(Pfun, start = list(shape = as.numeric(shape)))
      # shape = coef(fit)['shape']

      # PARETO: analytical MLE
      fit$shape = as.numeric(-1 * (sum(fit$z_nk[, 1])) / (fit$z_nk[, 1] %*% (log(fit$scale) - logX)))
      fit = mobster:::.set_params_Pareto(fit, fit$shape, fit$scale)

      # BETA: numerical MLE or analytical MM
      for (k in 2:fit$K)
      {
        if (fit.type == 'MLE')
          # MLE
        {
          # Compute a functional of the negative logLik, which we minimize
          MLE.fit = stats4::mle(
            minuslogl = mobster:::.NLLBetaMix(k, X$VAF, fit$z_nk, fit$pi),
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

          if (is.na(mean) | is.na(var)) {
            {
              warning('Possible singularity in one Beta component a/b --> Inf.')
            }
          }
          else {
            par = mobster:::.estBetaParams(mean, var)
            fit$a[k - 1] = par$a
            fit$b[k - 1] = par$b
          }

          names(fit$a) = names(fit$b) = names.BetaC
          fit = mobster:::.set_params_Beta(fit, fit$a, fit$b)
        }
      }

      ## Convergency test
      if (mobster:::.stoppingCriterion(i,
                             prevNLL,
                             fit$NLL,
                             prevpi,
                             fit$pi,
                             fit.type,
                             epsilon,
                             is_verbose,
                             fit$K))
        break


    } #End of Expectation Maximization loop.

    if(is_verbose) cat("EM completed!")

    # Status of convergence: TRUE/ FALSE
    fit$status = (i < maxIter)

    ########### Add names to the estimated variables for clarity, populate summary matrices
    names(fit$pi) =  colnames(fit$z_nk) = colnames(fit$pdf.w) = c(names.ParetoC, names.BetaC)
    names(fit$a) = names(fit$b) = names.BetaC

    # Update fit table
    fit = mobster:::.set_params_Beta(fit, fit$a, fit$b)
    fit = mobster:::.set_params_Pareto(fit, fit$shape, fit$scale)
    fit = mobster:::.set_params_Pi(fit, fit$pi)

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
    
    # print("APPLICO FILTRI")
    
    # Now apply the cluster-selection heuristic in function choose_clusters
    fit = choose_clusters(fit,  pi_cutoff = pi_cutoff, N_cutoff = N_cutoff)

    # print("RINOMINO")
    
    # ... and re-order the Beta cluster ID by mean ...
    fit = rename_Beta_clusters(fit)
    
    return(fit)
  }



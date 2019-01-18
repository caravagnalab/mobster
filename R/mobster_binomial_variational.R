#' Variational fit Dirichelt Mixture Model with Binomial data
#'
#' @param K
#' @param alpha_0
#' @param a_0
#' @param b_0
#' @param max_iter
#' @param epsilon_conv
#' @param restarts
#' @param is_verbose
#' @param iterative
#' @param parallel
#' @param x 
#' @param pi_cutoff 
#' @param silent 
#' @param q_init 
#' @param cores.ratio
#' @param trace
#'
#' @return
#' @export
#'
#' @examples
mobster_fit_binomial = function(x,
                                K = 2 * length(x$samples),
                                alpha_0 = 1e-6,
                                a_0 = 1,
                                b_0 = 1,
                                max_iter = 5000,
                                epsilon_conv = 1e-10,
                                pi_cutoff = 1e-2,
                                restarts = 10,
                                is_verbose = FALSE,
                                iterative = FALSE,
                                parallel = FALSE,
                                cores.ratio = .8,
                                silent = FALSE,
                                q_init = 'prior', # or kmeans, or private
                                trace = FALSE
)
{
  best = NULL
  
  pioHdr('MOBSTER Binomial mixture (variational multivariate fit)')
  
  # Get data
  x_DP = DP_table(x)
  x_NV = NV_table(x)
  
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")
  
  # Console header
  cat(bold("\n\tINPUT"))
  
  pioStr("\n  Points", paste0('N = ', nrow(x_DP)))
  pioStr("\nClusters", paste0('K = ', K), suffix = "(max)\n")
  
  pioStr("\nDirichlet", paste0('alpha = ', alpha_0), suffix = "(conc.)")
  pioStr("\n     Beta", paste0('a0 = ', a_0, '; b0 =', b_0), suffix = "(shape)\n")
  pioStr("\n     Beta (posterior)", ifelse(q_init, "With kmeans", "Prior"), suffix = "(variational distributions, `q``)")
  
  
  pioStr(
    "\n Optimize",
    paste0(
      'epsilon = ',
      epsilon_conv,
      '; steps =',
      max_iter,
      '; r = ',
      restarts
    )
  )
  
  # Extra parameters
  if (is_verbose)
    cat(green('\n ~ Verbose mode ...\n'))
  if (iterative)
    cat(green('\n ~ Press return at each step... \n'))
  if (parallel)
    cat(green('\n ~ Using parallel (output is suppressed) \n'))
  
  # Horrible scoping..
  single_model = function() {
    vb_bmm_MV(
      x_DP = x_DP,
      x_NV = x_NV,
      K = K,
      alpha_0 = alpha_0,
      a_0 = a_0,
      b_0 = b_0,
      max_iter = max_iter,
      epsilon_conv = epsilon_conv,
      is_verbose = is_verbose,
      iterative = iterative,
      silent = silent,
      trace = trace,
      q_init = q_init
    )
  }
  
  models = NULL
  best_ELBO = -Inf
  
  # Sequential implementation
  if (!parallel)
  {
    for (r in 1:restarts)
    {
      flush.console()
      
      pioStr('\n   #', r, suffix = '\n')
      
      obj = single_model() # fit
      
      if (max(obj$ELBO) > best_ELBO) {
        cat(bgGreen(black("New best ELBO.\n")))
        best_ELBO = max(obj$ELBO)
      }
      
      # Store results
      models = append(models, list(obj))
    }
  }
  else
    # Parallel implementation
  {
    # Setup clusters for parallel computing
    cl = .setup_parallel(cores.ratio = cores.ratio)
    
    # perform parallel inferences
    r = foreach(
      num = 1:restarts,
      .packages = "crayon",
      .export = ls(globalenv())
    ) %dopar%
    {
      obj = single_model()
    }
    
    models = r
    
    .stop_parallel(cl)
  }
  
  # Get best output a posteriori
  best = NULL
  best_ELBO = -Inf
  
  for (i in 1:restarts) {
    obj = models[[i]]
    
    if (max(obj$ELBO) > best_ELBO) {
      best = obj
      best_ELBO = max(obj$ELBO)
    }
  }
  
  # Print some output
  cat(bold("\n\nBEST FIT\n\n"))
  print(best)
  
  # Select clusters
  pioTit(
    paste0(
      "Selecting output clusters: requested minimum size is ",
      pi_cutoff
    ))
  
  best = choose_clusters(best, pi_cutoff)
  
  # Final output
  cat(bold(paste0("\n\n\tOUTPUT\n")))
  print(best)
  
  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")
  cat(
    bold("\n\nCOMPLETED: ") %+% cyan(round(TIME, 2), 'mins, with status', best$status, '\n')
  )
  
  x$fit.Binomial = best
  x = logOp(x, "Binomial clustering of read counts")
  
  return(x)
}


# Fit Variational Dirichelt Mixture Model with Binomial data
vb_bmm_MV <-
  function(x_DP,
           x_NV,
           K,
           alpha_0,
           a_0,
           b_0,
           max_iter,
           epsilon_conv,
           is_verbose,
           iterative,
           silent,
           trace = FALSE,
           q_init = 'prior' # or kmeans, or private
  )
  {
    # Use the log sum exp trick for having numeric stability
    log_sum_exp <- function(x) {
      # Computes log(sum(exp(x))
      offset <- max(x)
      s <- log(sum(exp(x - offset))) + offset
      i <- which(!is.finite(s))
      if (length(i) > 0) {
        s[i] <- offset
      }
      return(s)
    }
    
    # Check input
    stopifnot('id' %in% colnames(x_DP))
    stopifnot('id' %in% colnames(x_NV))
    stopifnot(ncol(x_NV) == ncol(x_DP))
    
    X = full_join(x_DP, x_NV, by = 'id')
    
    Sn = as.matrix(x_NV %>% select(-id))
    Tn = as.matrix(x_DP %>% select(-id))
    TminusSn = Tn - Sn
    
    ### Auxiliary variables
    N = nrow(X)              # Number of observations
    W = ncol(Sn)             # Number of dimensions (biopsies)
    L = -Inf                 # ELBO values
    
    cluster_names = paste0("C", 1:K)
    sample_names = x_DP$id
    dimensions_names = gsub('.NV', '', colnames(Sn))
    
    r_nk = log_r_nk = log_lambda_nk = matrix(0, nrow = N, ncol = K)
    rownames(r_nk) = rownames(log_r_nk) = rownames(log_lambda_nk) = sample_names
    colnames(r_nk) = colnames(log_r_nk) = colnames(log_lambda_nk) = cluster_names
    
    log_pi = rep(0, K)
    names(log_pi) = cluster_names
    
    ### Priors
    # Dirichlet [alpha_0 ... alpha_0]
    alpha_0 = rep(alpha_0, K)
    names(alpha_0) = cluster_names
    
    # Beta (a_0, b_0) with small randomness (~10e-3) to render different components
    a_0_scalar = a_0
    a_0 = rep(a_0, K * W) + runif(K * W) / 100
    b_0 = rep(b_0, K * W) + runif(K * W) / 100
    
    # .. we turn these into a KxW matrix
    a_0 = matrix(a_0, ncol = K)
    b_0 = matrix(b_0, ncol = K)
    
    rownames(a_0) = rownames(b_0) = dimensions_names
    colnames(a_0) = colnames(b_0) = cluster_names
    
    # Posterior approximations, 3 different ways to initialize the 
    # variational distribution q (via `q_init`)
    # - random : q_init = to the prior's value 
    # - kmeans : via kmeans clustering
    # - private: assigning clusters to private mutations (plus the prior)
    alpha = alpha_0
    
    if(q_init == 'prior')
    {
      pio::pioStr("Initializing variational distribution", paste0(q_init, " (=prior)"))
      a = a_0
      b = b_0
    }
    
    if(q_init == 'kmeans') # Special case: kmeans initialization
    {
      pio::pioStr("Initializing variational distribution", paste0(q_init, " (kmeans)"))
      
      # get parameters from kmeans
      kmeans_params = initial_condition_Binomial_kmeans(x_NV, x_DP, K, a_0_scalar)
      
      a = kmeans_params$a
      b = kmeans_params$b
    }  
    
    if(q_init == 'private') # Special case: prior plus private clusters
    {
      pio::pioStr("Initializing variational distribution", paste0(q_init, " (private clusters, plus kmeans)"))
      
      # get private clusters, 
      prv_params = initial_condition_prv_clusters(x_NV, x_DP)
      K_prv_params = ncol(prv_params$a)
      
      a = prv_params$a
      b = prv_params$b
      
      # other clusters are requested, so we take them with kmeans
      if(K > K_prv_params) 
      {
        X_not_ass = prv_params$X %>% 
          as_tibble() %>%
          filter(is.na(private.cluster)) %>% pull(id)
        
        kmeans_params = initial_condition_Binomial_kmeans(
          x_NV %>% filter(id %in% X_not_ass), 
          x_DP %>% filter(id %in% X_not_ass), 
          K - K_prv_params, 
          a_0_scalar)
        
        a = cbind(a, kmeans_params$a)
        b = cbind(b, kmeans_params$b)
      }
    }  
    
    ### normalization constants for the priors, required for the ELBO
    # ln Cn -- log of Binomial data
    # ln C(alpha_0) -- log of the Dirichlet
    # ln Beta(a0, b0) -- log of K Beta
    log_Cn = lchoose(Tn, Sn)
    log_C_alpha_0 = lgamma(sum(alpha_0)) - sum(lgamma(alpha_0))
    log_beta_ab_0 = lbeta(a, b)
    
    # Computation status to be updated in the end
    status = ''
    
    # The trace of the fit
    fit_trace = NULL
    
    # Iterate to find optimal parameters
    i = 1
    repeat {
      # Variational E-Step
      
      # Digamma values for the current Beta distributions
      # that we need in the E-Step
      dga = digamma(a)      # Psi(a)
      dgb = digamma(b)      # Psi(b)
      dgab = digamma(a + b) # Psi(a+b)
      
      # Update expectations
      log_pi = digamma(alpha) - digamma(sum(alpha)) # ln \pi -- Dirichlet entropy
      log_theta = dga - dgab                        # ln \theta
      log_1min_theta = dgb - dgab                   # ln 1-\theta
      
      # Log of lambda terms
      for (k in 1:K)
        log_lambda_nk[, k] <-
        rowSums(log_pi[k] + sapply(1:W, function(w)
          Sn[, w] * log_theta[w, k] + TminusSn[, w] * log_1min_theta[w, k]))
      # log_pi[k] + Sn * log_theta[k] + TminusSn * log_1min_theta[k]
      
      
      # Sn %*% log_theta[, k]
      #
      # matrix(1:8, ncol = 2, byrow = T) %*% c(2,1)
      
      
      
      # Calculate probabilities using the logSumExp trick for numerical stability
      Z        <- apply(log_lambda_nk, 1, log_sum_exp)
      log_r_nk <- log_lambda_nk - Z              # log of r_nk
      r_nk     <- apply(log_r_nk, 2, exp)        # r_nk
      
      # Save this entry in the trace
      if(trace) {
        fit_trace = bind_rows(fit_trace,
                              tibble(
                                cluster.Binomial = mobster:::latent_vars_hard_assignments(lv = list(`z_nk` = r_nk, `pi` = alpha_0)),
                                step = i))
      }
      
      
      
      # Variational M-Step
      
      # Auxiliary quantities
      N_k <- colSums(r_nk) + 1e-10
      s_star = sapply(1:K,
                      function(i) {
                        sapply(1:W, function(w)
                          return(Sn[, w] %*% r_nk[, i]))
                      })
      
      t_star =  sapply(1:K,
                       function(i) {
                         sapply(1:W, function(w)
                           return(Tn[, w] %*% r_nk[, i]))
                       })
      
      t_star = matrix(t_star, ncol = ncol(a))
      s_star = matrix(s_star, ncol = ncol(a))
      
      colnames(t_star) = colnames(s_star) = colnames(a)
      rownames(t_star) = rownames(s_star) = rownames(a)
      
      # Update equations
      alpha = alpha_0 + N_k       # Dirichlet
      a = a_0 + s_star            # Beta
      b = b_0 + t_star - s_star   # Beta
      
      # Compute posterior's expected values of model's parameter
      pi_k <-
        (alpha_0 + N_k) / (K * alpha_0 + N)  # mixing proportions
      theta_k = a / (a + b)                  # Binomial parameter (mean of the posterior)
      
      ###### ELBO for convergency
      
      # Updated values for the new posterior estimates
      dga = digamma(a)        # Psi(a)
      dgb = digamma(b)        # Psi(b)
      dgab = digamma(a + b)   # Psi(a+b)
      log_beta_ab = lbeta(a, b)  # Beta(a,b)
      
      # These are the quantities for the ELBO equation
      hat_digamma_ab = dga - dgb
      hat_digamma_ba = dgb - dgab
      log_rho = digamma(alpha) - digamma(sum(alpha))
      
      eta_k = lapply(1:W,
                     function(w)
                       log_beta_ab[w,] - log_beta_ab_0[w,] + (a_0[w,] - a[w,]) * dga[w,] + (b_0[w,] - b[w,]) * dgb[w,] + (a[w,] - a_0[w,] + b[w,] - b_0[w,]) * dgab[w,])
      eta_k = Reduce(rbind, eta_k)
      eta_k = matrix(eta_k, ncol = K)
      
      rownames(eta_k) = rownames(a)
      colnames(eta_k) = names(pi_k)
      
      # Dirichlet normalization constant -- for numerical stability in log format
      log_C_alpha = lgamma(sum(alpha)) - sum(lgamma(alpha))
      
      # The data likelihood
      my_dotprodfun = function(k) {
        dmdotprod = sapply(1:W, function(w)
          Sn[, w] * hat_digamma_ab[w, k] + Tn[, w] * hat_digamma_ba[w, k] - log_Cn[, w] + log_rho[k] - log_r_nk[, k])
        dmdotprod = rowSums(dmdotprod)
        
        return(r_nk[, k] %*% dmdotprod)
      }
      
      # The actual ELBO -- we round it to the 9th digit to avoid some weird numerical issues
      ELBO = log_C_alpha_0 - log_C_alpha + sum(sapply(1:K, my_dotprodfun)) + (alpha_0 - alpha) %*% (log_rho) + sum(colSums(eta_k))
      ELBO = round(ELBO, 9)
      L = c(L, ELBO)
      
      # if (is_verbose)
      # {
      #   cat('\n**** POSTERIOR APPROX. #', (i - 1), '\n')
      #
      #   df = cbind(round(a, 4),
      #              round(b, 4),
      #              round(alpha, 4),
      #              round(pi_k, 4),
      #              round(theta_k, 4))
      #
      #   colnames(df) = c('a',
      #                    'b',
      #                    'alpha',
      #                    'E[pi_k]',
      #                    'E[theta_k]')
      #   rownames(df) = paste('Mixture', 1:K)
      #
      #   df = df[order(df[, 'E[pi_k]'], decreasing = TRUE), ]
      #   print(df)
      # }
      
      if (i > 1 & !silent) {
        cat("\r")
        cat(sprintf('It#: %-8s \t ELBO: %-20s\t Diff: %-20s', i,  L[i], abs(L[i] - L[i -
                                                                                       1])))
      }
      
      # Checks -- decreasing ELBO
      if (i > 1 && L[i] < L[i - 1])
        warning("ELBO decreased by ", (L[i - 1]-L[i]) ,"at step #", i)
      
      # Checks -- convergence by epsilon_conv
      if (!is.na(epsilon_conv) && i > 1 &&
          abs(L[i] - L[i - 1]) < epsilon_conv) {
        if (!silent)
          cat(green("\nConverged with eps.", epsilon_conv, "at step", i))
        
        status = 'CONVERGED'
        flush.console()
        
        break
      }
      
      # Checks -- convergence by max_iter
      if (!is.na(max_iter) && i == max_iter) {
        if (!silent)
          cat(red("\nVB did not converged and was stopped at step", max_iter))
        
        status = 'INTERRUPTED'
        
        break
      }
      
      if (iterative)
        readline(paste(i, ':: (Hit ret)'))
      
      if (!silent)
        cat("\r")
      i = i + 1
    }
    
    cat('\n')
    
    # MAP clustering assignment, and add it to data together with the ID
    labels = tibble(cluster.Binomial = mobster:::latent_vars_hard_assignments(lv = list(`z_nk` = r_nk, `pi` = pi_k)))
    X = bind_cols(labels, X)
    
    if(trace) {
      fit_trace = bind_rows(fit_trace,
                            tibble(
                              cluster.Binomial = mobster:::latent_vars_hard_assignments(lv = list(`z_nk` = r_nk, `pi` = alpha_0)),
                              step = i))
    }
    
    
    obj <-
      structure(
        list(
          X = X,
          K = K,
          N = N,
          pi_k = pi_k,
          theta_k = theta_k,
          alpha = alpha,
          r_nk = r_nk,
          labels = labels,
          a = a,
          b = b,
          a_0 = a_0,
          b_0 = b_0,
          alpha_0 = alpha_0,
          epsilon_conv = epsilon_conv,
          status = status,
          ELBO = L,
          trace = fit_trace
        ),
        class = "vb_bmm",
        call = match.call(),
        modelname = "Variational Bayes fit for Dirichlet-Binomial Mixture-Models (VBDBMM) -- 1D_data"
      )
    
    return(obj)
  }


############## S3 functions for the object

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.vb_bmm <- function(x, ...) {
  print.vb_bmm(x, ...)
}

#' Title
#'
#' @param x
#' @param pi.cutoff
#'
#' @return
#' @export
#'
#' @examples
print.vb_bmm <- function(x, pi.cutoff = 0.01) {
  stopifnot(inherits(x, "vb_bmm"))
  
  pioHdr('MOBSTER Binomial mixture')
  
  pio::pioStr("\n Points", paste0("N = ",  nrow(x$X)))
  pio::pioStr("\nSamples", paste0("W = ",  nrow(x$theta_k)))
  pio::pioStr("\n\n Status",  paste0(x$status, ' (', length(x$ELBO), ' steps)'))
  
  pio::pioStr("\n\nBinomial parameters", '\n')
  
  pio::pioDisp(round(x$theta_k, 2))
  
  pio::pioStr("\nProportions", '', suffix = '\n')
  
  pio::pioDisp(x$pi_k)
}



########### Private function

# Subset clusters with freq below a cutoff -- maintain the object consistency
choose_clusters = function(x, pi.cutoff)
{
  table_pi = tibble(cluster = names(x$pi_k), pi = x$pi_k)
  table_pi = table_pi %>%
    mutate(accepted = pi > pi.cutoff) %>%
    arrange(desc(pi))
  
  # partitions of clusters
  detect.clones = table_pi %>% filter(pi > pi.cutoff)
  rejected.clones = table_pi %>% filter(pi <= pi.cutoff)
  
  K = nrow(detect.clones)
  K.rj = nrow(rejected.clones)
  
  # Mapping old labels to new ones
  table_pi$new.labels = paste0('C', 1:nrow(table_pi))
  
  mapping = table_pi$new.labels
  names(mapping) = table_pi$cluster
  
  # rename all entries with reference to the labels
  x$X$Binomial.cluster = mapping[x$X$cluster.Binomial]
  
  names(x$alpha) = names(x$pi_k) = names(x$alpha_0) = mapping[names(x$alpha_0)]
  
  colnames(x$a) = colnames(x$b) = colnames(x$a_0) = colnames(x$b_0) =
    colnames(x$r_nk) = colnames(x$theta_k) =
    mapping[colnames(x$theta_k)]
  
  # print some info about what has been selected
  table_pi = table_pi %>% mutate(cluster = new.labels) %>% select(-new.labels)
  
  pio::pioStr("\nNumber of input clusters", paste0('K = ', x$K))
  pio::pioStr(
    "\nFit",
    paste0('K = ', K),
    prefix = '\n',
    suffix = paste0('[pi > ', pi.cutoff, ']\n\n')
  )
  
  print(table_pi)
  
  # now drop useless clusters and re-normalize the variables
  detect.clones = table_pi %>% filter(pi > pi.cutoff) %>% pull(cluster)
  
  x$alpha = x$alpha[detect.clones]
  x$pi_k = x$pi_k[detect.clones]
  x$alpha_0 = x$alpha_0[detect.clones]
  x$alpha = x$alpha[detect.clones]
  
  x$a = x$a[, detect.clones, drop = FALSE]
  x$b = x$b[, detect.clones, drop = FALSE]
  x$a_0 = x$a_0[, detect.clones, drop = FALSE]
  x$b_0 = x$b_0[, detect.clones, drop = FALSE]
  x$r_nk = x$r_nk[, detect.clones, drop = FALSE]
  x$theta_k = x$theta_k[, detect.clones, drop = FALSE]
  
  # renormalize latent variables
  C = rowSums(x$r_nk)
  for (i in 1:nrow(x$r_nk))
    x$r_nk[i,] = x$r_nk[i,] / C[i]
  
  # renormalize mixing proportions
  C = sum(x$pi_k)
  x$pi_k = x$pi_k / C
  
  # recompute clustering assignments..
  labels = mobster:::latent_vars_hard_assignments(lv = list(`z_nk` = x$r_nk,
                                                            `pi` = x$pi_k))
  
  x$X$cluster.Binomial = labels
  
  # Update num of clusters
  x$K = K
  
  x
}

initial_condition_prv_clusters = function(x_NV, x_DP)
{
  # reconstruct thw raw VAF profile for each sample
  vaf = (x_NV %>% select(-id))/(x_DP %>% select(-id))
  colnames(vaf) = gsub('.NV', '', colnames(vaf))
  vaf[is.na(vaf)] = 0
  
  ret_vaf = vaf %>% as_tibble()
  ret_vaf$id = x_NV$id
  ret_vaf$private.cluster = NA
  
  W = ncol(vaf)
  
  # Private clusters have only one positive entry in the vaf
  prv = apply(vaf, 1, function(w) sum(w>0)) == 1
  ret_vaf$private.cluster = apply(ret_vaf, 1,
                                  function(w)
                                  {
                                    pos_entries = as.numeric(w[1:W]) > 0
                                    
                                    if (sum(pos_entries) > 1)
                                      return(NA)
                                    else
                                      return(colnames(ret_vaf)[which(pos_entries)])
                                  })
  
  prv = vaf[prv, ]
  
  
  # Make a list of params
  prv = lapply(1:ncol(prv), function(w) unlist(prv[prv[, w] > 0, w]))
  names(prv) = fitgp$samples
  
  # Where we have the 'a' and the "b", for each mean and variance value
  prv_beta_params = rbind(
    sapply(prv, mean),
    sapply(prv, var))
  
  prv_size = sapply(prv, length)
  
  # Assemble everything into a dataframe
  prv_beta_params = apply(prv_beta_params, 2, function(w) data.frame(mobster:::.estBetaParams(w[1], w[2])))
  prv_beta_params = Reduce(rbind, prv_beta_params)
  rownames(prv_beta_params) = names(prv)
  
  prv_beta_params = cbind(prv_beta_params, n = prv_size)
  
  # Now, we create the matrix of a and b parameters for the Beta where each dimension has
  # pdf concentrated on 0 (b >> a), and the private parameters per cluster
  
  a_Beta = b_Beta = matrix(NA, ncol = nrow(prv_beta_params), nrow = W)
  rownames(a_Beta) = rownames(b_Beta) = rownames(prv_beta_params)
  
  a_Beta = apply(a_Beta, c(1,2), function(w) 1+runif(1)) 
  b_Beta = apply(a_Beta, c(1,2), function(w) 30+runif(1))   
  
  
  for(i in 1:nrow(prv_beta_params)) {
    # Get a and b for a single dimension
    a_Beta[i, i] = prv_beta_params[i, 'a']
    b_Beta[i, i] = prv_beta_params[i, 'b']
  }
  
  return(list(a=a_Beta,b=b_Beta, X=ret_vaf))
}

initial_condition_Binomial_kmeans = function(x_NV, x_DP, K, a_0)
{
  # Data
  XVAF = (x_NV %>% select(-id))/(x_DP %>% select(-id)) %>% as_tibble()
  colnames(XVAF) = gsub('.NV', '', colnames(XVAF))
  XVAF[is.na(XVAF)] = 0
  
  W = ncol(XVAF)
  
  # Kmeans with K clusters
  km_cl = kmeans(XVAF, centers = K, nstart = 25)
  XVAF$cluster.kmeans = km_cl$cluster
  
  # Beta means as Gaussian means
  beta_mu = km_cl$centers[, 1:W]
  
  # Beta variances from posterior assignments
  beta_var = XVAF %>%
    reshape2::melt(id = 'cluster.kmeans') %>% as_tibble() %>%
    arrange(cluster.kmeans) %>%
    group_by(cluster.kmeans, variable) %>% summarise(value = var(value))  %>%
    spread(variable, value) %>%
    ungroup() %>%
    select(-cluster.kmeans)
  
  a_0 = rep(a_0, K * W) + runif(K * W) / 100
  a_0 = matrix(a_0, ncol = K)
  
  a = b = a_0
  
  beta_mu = t(beta_mu)
  beta_var = t(beta_var)
  
  for(i in 1:nrow(beta_mu))
  {
    for(j in 1:ncol(beta_mu))
    {
      a[i,j] = mobster:::.estBetaParams(beta_mu[i,j], beta_var[i,j])$a
      b[i,j] = mobster:::.estBetaParams(beta_mu[i,j], beta_var[i,j])$b
    }
  }
  
  extr = function(x){
    if(is.na(x)) return(runif(1) / 100)
    x
  }
  a = apply(a, c(1,2), extr)
  b = apply(b, c(1,2), extr)
  
  return(list(a=a, b=b, X = XVAF %>% as_tibble()))
}


# 
# fit_trace
# Xr = VAF_table(x)
# Xr$cluster = fit_trace %>% filter(step == 46) %>% pull(cluster.Binomial)
# 
# ggplot(Xr, aes(x=Set7_55.VAF, y= Set7_57.VAF, color = factor(cluster)))+
#   geom_point(alpha = .3) + 
#   facet_wrap(~factor(cluster)) 
# # geom_point(data = bp, color = 'black') +
# # guides(color = FALSE)

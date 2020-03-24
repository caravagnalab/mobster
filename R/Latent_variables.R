# Return latent variables (re-computed), plus other variables.
latent_vars = function(x, clusters_tibble = NULL) {

  if(all(!is.null(clusters_tibble)))
    x$Clusters = clusters_tibble

  # Extract Beta values
  B = mobster:::.params_Beta(x)
  a = B$a
  b = B$b

  # print(B)
  # print(a)
  # print(b)
  # print(x$Clusters)
  
  # if(is.null(a)) save(x, file = '~/G.RData')
  # if(is.null(a)) cat(" -> NULL ", x$K)
    
  
  names(a) = names(b) = B$cluster

  # Extract Tail values
  shape = ifelse(x$fit.tail, mobster:::.params_Pareto(x)$Shape, NA)
  scale = ifelse(x$fit.tail, mobster:::.params_Pareto(x)$Scale, NA)

  # Extract Mixing Proportions
  pi = mobster:::.params_Pi(x)

  # Reconstruct latent variables, responsibilities and weighted densities
  N = nrow(x$data)
  z_nk  = matrix(0, nrow = N, ncol = x$K)
  colnames(z_nk) = names(pi)
  
  pdf.w = z_nk

  # Get log density per component
  for (k in 1:x$K)
    pdf.w[, k] = ddbpmm(x,
                        data = x$data$VAF,
                        components = paste(k),
                        a = a,
                        b = b,
                        pi = pi,
                        shape = shape,
                        scale = scale,
                        log = TRUE)

  # Calculate probabilities using the logSumExp trick for numerical stability
  Z = apply(pdf.w, 1, .log_sum_exp)
  z_nk   = pdf.w - Z
  z_nk   = apply(z_nk, 2, exp)

  NLL    = -sum(Z)  # Evaluate the NLL

  return(list(NLL = NLL, z_nk = z_nk, pdf.w = pdf.w, a = a, b = b, pi = pi, shape = shape, scale = scale))
}

# Compute hard clustering assignments
latent_vars_hard_assignments = function(lv) {
  # Cluster labels of each data point. Each data point is assigned to the cluster
  # with the highest posterior responsibility. Hard assignment.
  unlist(apply(lv$z_nk, 1,
               function(x) {
                 names(lv$pi)[which(x == max(x, na.rm = TRUE))[1]]
               }))
}

# Compute scores for model selection 
latent_vars_scores = function(lv, K, tail, cluster)
{
  ##==============================
  # Scores for model selection   #
  ##==============================
  
  # K = Kbeta + 1
  pi_components = (K - 1) + as.numeric(tail)
  Beta_components = (K-1) * 2
  Pareto_components = 2 * as.numeric(tail)
  
  numParams = pi_components + Beta_components + Pareto_components
  
  # Kbeta = K + as.numeric(tail)
  Kbeta = K + as.numeric(tail)
  
  N = nrow(lv$z_nk)
  
  # if (tail)
  #   numParams = K + 2 * Kbeta + 1          # Total number of parameters i.e. pi (Beta + Pareto) + 2 * Kbeta (Beta) + 1 (Pareto)
  # else
  #   numParams = (K - 1) + 2 * Kbeta         # Total number of parameters i.e. pi (Beta) + 2 * Kbeta (Beta)
  
  BIC <-
    2 * lv$NLL + numParams * log(N)     # BIC = -2*ln(L) + params*ln(N)
  AIC <-
    2 * lv$NLL + 2 * numParams              # AIC = -2*ln(L) + 2*params
  
  # Integrated Complete Likelihood criterion -- uses standard entropy
  entropy <- -sum(lv$z_nk * log(lv$z_nk), na.rm = TRUE)
  ICL <- BIC + entropy
  
  # Integrated Complete Likelihood criterion with reduced entropy (only for
  # latent variable that involve subclones -- i.e., Betas). I think this is also a sort of
  # conditional entropy where we condition on the MAP estimate of a mutation being part
  # of a clone, rather than tail.
  
  # HTake the MAP estimates of z_nk, and select only entries that are assigned
  # to a Beta component (i.e. those with arg_max != Tail)
  cz_nk = lv$z_nk[cluster != 'Tail', 2:ncol(lv$z_nk), drop = FALSE]
  
  # print("BEFORE")
  # print(cz_nk %>% head)
  
  # This is un-normalized -- we compute the empirical normalizing constant (C)
  C = rowSums(cz_nk)
  # for (i in 1:nrow(cz_nk))
  #   cz_nk [i, ] = cz_nk [i,  ] / C[i]
  # 
  
  cz_nk = cz_nk/C
  # print("A1")
  
  # The reduced entropy is the entropy of this distribution
  rentropy = -sum(cz_nk  * log(cz_nk), na.rm = TRUE)
  # print("A2")
  
  # Integrated Complete Likelihood criterion with reduced entropy
  reICL <- BIC + rentropy
  
  scores = data.frame(
    NLL = lv$NLL,
    BIC = BIC,
    AIC = AIC,
    entropy = entropy,
    ICL = ICL,
    reduced.entropy = rentropy,
    reICL = reICL,
    size = numParams
  )
  
  scores
}


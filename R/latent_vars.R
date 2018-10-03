latent_vars = function(x, clusters_tibble = NULL) {
  
  if(all(!is.null(clusters_tibble))) 
    x$Clusters = clusters_tibble
  
  # Extract Beta values
  B = mobster:::.params_Beta(x)
  a = B$a
  b = B$b
  
  names(a) = names(b) = rownames(B)
  
  # Extract Tail values
  shape = ifelse(x$fit.tail, mobster:::.params_Pareto(x)$Shape, NA)
  scale = ifelse(x$fit.tail, mobster:::.params_Pareto(x)$Scale, NA)
  
  # Extract Mixin Prop
  pi = mobster:::.params_Pi(x)  
  
  # Reconstruct latent variables, responsibilities and weighted densities
  N = nrow(x$data)
  z_nk  = matrix(0, nrow = N, ncol = x$K)       
  pdf.w = z_nk
  
  # Get density per component
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
  Z = apply(pdf.w, 1, mobster:::.log_sum_exp)
  z_nk   = pdf.w - Z
  z_nk   = apply(z_nk, 2, exp)    
  
  NLL    = -sum(Z)  # Evaluate the NLL
  
  return(list(NLL = NLL, z_nk = z_nk, pdf.w = pdf.w, a = a, b = b, pi = pi, shape = shape, scale = scale))
}

# lv = latent_vars(x, clusters_tibble)
# 
latent_vars_hard_assignments = function(lv) {
  
  # Cluster labels of each data point. Each data point is assigned to the cluster
  # with the highest posterior responsibility. Hard assignment.
  unlist(apply(lv$z_nk, 1,
               function(x) {
                 names(lv$pi)[which(x == max(x, na.rm = TRUE))[1]]
               }))
  
}

latent_vars_scores = function(lv, K, tail, cluster)
{
  ##==============================
  # Scores for model selection   #
  ##==============================
  
  Kbeta = K + as.numeric(tail)
  N = nrow(lv$z_nk)

  if (tail)
    numParams = K + 2 * Kbeta + 1          # Total number of parameters i.e. pi (Beta + Pareto) + 2 * Kbeta (Beta) + 1 (Pareto)
  else
    numParams = (K - 1) + 2 * Kbeta         # Total number of parameters i.e. pi (Beta) + 2 * Kbeta (Beta)
  
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
  
  # This is un-normalized -- we compute the empirical normalizing constant (C)
  C = rowSums(cz_nk)
  for (i in 1:nrow(cz_nk))
    cz_nk [i, ] = cz_nk [i, ] / C[i]
  
  # The reduced entropy is the entropy of this distribution
  rentropy = -sum(cz_nk  * log(cz_nk), na.rm = TRUE)
  
  # Integrated Complete Likelihood criterion with reduced entropy
  reICL <- BIC + rentropy
  
  scores = data.frame(
    NLL = lv$NLL,
    BIC = BIC,
    AIC = AIC,
    entropy = entropy,
    ICL = ICL,
    reICL = reICL,
    size = numParams
  )
  
  scores
}
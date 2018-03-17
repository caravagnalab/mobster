find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}


mmdensity = function(pi, p, K, trials)
{
  # Trials, 1 to ...
  x = 1:trials
  # x = round(seq(1, trials, trials/40))
  
  # Simulate a Binomial sampling with "trials", and report the probability of having "x" successes
  # according to the mixture
  density = NULL
  mixture = 0.0
  for (i in 1:K)
  {
    # Calculate the estimated density (rescaled for the proprortions)
    mixture.points = dbinom(x, size = trials, prob = p[i]) * pi[i]
    density = append(density, list(mixture.points))
    mixture <- mixture + mixture.points
  }
  
  return(list(
    x = x,
    mixture = mixture,
    density = density
  ))
}

dbmixture = function(x, pi, p, K, Log = TRUE)
{
  density = 0
  for (i in 1:K)
  {
    # Calculate the estimated density (rescaled for the proprortions)
    # density = density + sum(dbinom(x$s_n, size = x$t_n, prob = p[i], log = TRUE) + log(pi[i]))
    density = density + sum(dbinom(x$s_n, size = x$t_n, prob = p[i]) * pi[i])
  }
  
  if(Log) density = log(density)
  
  # print(density)
  
  return(density)
}


kldivergence = function(p, q, trials, lowerbound.precision = 1e-10)
{
  p_density = mmdensity(p$pi, p$p, p$K, trials)$mixture
  q_density = mmdensity(q$pi, q$p, q$K, trials)$mixture
  
  
  nozeros = which(p_density > lowerbound.precision)
  p_density = p_density[nozeros]
  q_density = q_density[nozeros]
  
  nozeros = which(q_density > lowerbound.precision)
  p_density = p_density[nozeros]
  q_density = q_density[nozeros]
  
  # print(nozeros)
  # print(p_density)
  # print(q_density)

  if(length(p_density) == 0) return(Inf)
  
  return(as.numeric(p_density %*% log(p_density / q_density)))
}

clusters_precision = function(true_model, fit, pi.cutoff = 0.01, p.collapse = 0.02, p.range = 0.02)
{
  # Throw away mixtures at proportions below cutoff
  model_compress = function(m, c)
  {
    n = NULL
    active.modes = which(m$pi_k > c)
    n$pi_k = m$pi_k[active.modes]
    n$theta_k = m$theta_k[active.modes]
    n$K = length(m$theta_k) 
    
    return(n)
  }
  
  # merge mixtures components within c-breaks
  model_merge = function(m, c)
  {
    ival = seq(min(m$theta_k) - c, max(m$theta_k) + c, c)
    mapping = floor(m$theta_k / c)
    
    df = data.frame(m$theta_k, m$pi_k, mapping)
    df = split(df, f = mapping)
    
    collapsed.th = collapsed.pi.k = c()
    for(t in 1:length(df))
    {
      group = df[[t]]
      
      collapsed.th = c(collapsed.th, mean(group$m.theta_k))
      collapsed.pi.k = c(collapsed.pi.k, sum(group$m.pi_k))
    }
    
    collapsed.K = length(collapsed.th)
    
    n = NULL
    n$pi_k = collapsed.pi.k
    n$theta_k = collapsed.th
    n$K = collapsed.K
    
    return(n)
  }
  
  cat('** Comparator:',  pi.cutoff,  'min.pi --', p.collapse, 'binning --', p.range, 'tolerance.\n')
  
  # Compress models
  fit = model_compress(fit, pi.cutoff)
  true_model = model_compress(true_model, pi.cutoff)
  
  # Merge models (which renders easy to identify mixtures)
  fit = model_merge(fit, p.collapse)
  true_model = model_merge(true_model, p.collapse)
  
  # get dataframes
  tm = data.frame(theta_k = true_model$theta_k, pi_k = true_model$pi_k)
  fm = data.frame(theta_k = fit$theta_k, pi_k = fit$pi_k)
  
  
  tm = tm[order(tm$theta_k, decreasing = TRUE), , drop = F]
  fm = fm[order(fm$theta_k, decreasing = TRUE), , drop = F]
  
  # rounding to the 2nd digit renders comparisons less strict
  tm = data.frame(apply(tm, 2, round, digits = 5))
  fm = data.frame(apply(fm, 2, round, digits = 5))
  
  # orrible R's autoconversions
  if(fit$K == 1) fm = t(fm)
  if(true_model$K == 1) tm = t(tm)
  
  v.tm = tm[, 1]
  v.fm = fm[, 1]
  
  cat('@TM ', toStringCl(true_model), '\n')
  cat('@FIT', toStringCl(fit), '\n')
  
  # check if it's in a range with error tolerance r
  isinRange = function(x,y,r){return(any(abs(x-y) <= r))}
  
  which.tp = sapply(v.fm, isinRange, y = v.tm, r = p.range)
  which.fn = !sapply(v.tm, isinRange, y = v.fm, r = p.range)
  
  print(v.tm)
  print(v.fm)
  
  
  tp = v.fm[which.tp]
  fp = v.fm[!which.tp]
  fn = v.tm[which.fn]
      
  ntp = length(which(which.tp))
  nfp = length(v.fm) - ntp
  nfn = length(which(which.fn))
  
  cat('\nMixture parameters\n')
  cat('TP = ', ntp, ' : ', tp, '\n')
  cat('FP = ', nfp, ' : ', fp, '\n')
  cat('FN = ', nfn, ' : ', fn, '\n')
  
  
  return(list(
    ntp = ntp,   # num of TPs 
    nfp = nfp,   # num of FPs 
    nfn = nfn,   # num of FNs
    tp = tp,     # TPs
    fp = fp,     # FPs 
    fn = fn,     # FNs
    K_true = true_model$K,      # K (true model)
    K_fit = fit$K               # K (fit)
    ))
}


# analyze_experiment = function(
#   dir.output = '.',
#   dir.data = '.',
#   file.output = 'vb_bmm1D_fit', 
#   descr = '', 
#   NUM_TRIALS_PER_CONFIGURATION, 
#   COVERAGE_SPAN, 
#   NUM_OF_CLONES_SPAN,
#   pi.cutoff = 0.01, 
#   p.collapse = 0.02, 
#   p.range = 0.02
#   )
# {
#   COVERAGE_SPAN = as.character(COVERAGE_SPAN)
#   NUM_OF_CLONES_SPAN = as.character(NUM_OF_CLONES_SPAN)
#   
#   numtests =  length(COVERAGE_SPAN) * length(NUM_OF_CLONES_SPAN) * NUM_TRIALS_PER_CONFIGURATION
#   cat('** Analyze experiment: processing', numtests,' files in 2 second ...\n')
#   Sys.sleep(2)
#   
#   df_tp = df_fp = df_fn = df_k = NULL
#     
#   for (C in COVERAGE_SPAN)
#   {
#     for (K in NUM_OF_CLONES_SPAN)
#     {
#       for (i in 1:NUM_TRIALS_PER_CONFIGURATION)
#       {
#         fit.fname = paste(dir.output, '/', file.output, '.', descr,'.EXP-C', C, '-K', K, '-TRIAL',i, '.RData', sep = '')
#         cat('\n * ', fit.fname, '\n')
#         load(fit.fname)
# 
#         tm.fname = paste(dir.data, '/', 'EXP-C', C, '-K', K, '-TRIAL',i, '.RData', sep = '')
#         load(tm.fname)
#         
#         stats = clusters_precision(data, fit, pi.cutoff, p.collapse, p.range)
#         
#         df_tp = rbind(df_tp, data.frame(coverage = C, K = K, tp = stats$ntp))
#         df_fp = rbind(df_fp, data.frame(coverage = C, K = K, fp = stats$nfp))
#         df_fn = rbind(df_fn, data.frame(coverage = C, K = K, fn = stats$nfn))
#         df_k = rbind(df_k, data.frame(coverage = C, K = K, K_ratio = stats$K_true/stats$K_fit))
#       }
#     }
#   }
#   
#   return(list(df_tp = df_tp, df_fp = df_fp,  df_fn = df_fn,  df_k = df_k))
# }


analyze_experiment = function(
  dir.output = '.',
  dir.data = '.',
  file.output = 'vb_bmm1D_fit', 
  descr = '', 
  pi.cutoff = 0.01, 
  p.collapse = 0.02, 
  p.range = 0.02
)
{
  data.files = list.files(path = dir.data, pattern = '*.RData')
  fit.files = list.files(path = dir.output, pattern = '*.RData')
  
  if(!all(endsWith(fit.files, data.files))) warning('Files in data and fit folders do not match.')
  
  cat('** Analyze experiment: processing', length(fit.files),' files in 2 second ...\n')
  Sys.sleep(2)
  
  df = NULL
  
  missing.files = NULL
  
  for (f in 1:length(data.files))
  {
    cat('\n * ', data.files[f], '\n')
    load(paste(dir.data, '/', data.files[f], sep = ''))
    
    fitfile = paste(dir.output, '/', file.output, '.', descr, '.', data.files[f], sep = '') 
    if(!file.exists(fitfile)){
      missing.files = c(missing.files, fitfile)
      next;
    }
    load(fitfile)
        
    stats = clusters_precision(data, fit, pi.cutoff, p.collapse, p.range)
    
    df = rbind(
      df, 
      data.frame(
        coverage = data$coverage, 
        K = data$K, 
        purity = data$purity,
        tp = stats$ntp,
        fp = stats$nfp,
        fn = stats$nfn,
        K_ratio = stats$K_true/stats$K_fit
      ))
  }
  
  df = df[order(df$coverage), ]
  df = df[order(df$K), ]
  df = df[order(df$purity), ]
  
  # df$K = as.character(df$K)    
  # df$coverage = as.character(df$coverage)    
  # df$purity = as.character(df$purity)    
  df$K = as.factor(df$K)    
  df$coverage = as.factor(df$coverage)    
  df$purity = as.factor(df$purity)    
  
  rownames(df) = NULL
  
  cat('\n\n** Sample size: ', nrow(df), '\n')
  print(table(df[, c('coverage', 'K')]))
  
   
  return(df)
}


analyze_experiment_PopGen = function(
  dir.output = '.',
  dir.data = '.',
  file.output = 'vb_bmm1D_fit', 
  descr = '', 
  pi.cutoff = 0.01, 
  p.collapse = 0.02, 
  p.range = 0.02
)
{
  data.files = list.files(path = dir.data, pattern = '*.RData')
  fit.files = list.files(path = dir.output, pattern = '*.RData')
  
  if(!all(endsWith(fit.files, data.files))) warning('Files in data and fit folders do not match.')
  
  cat('** Analyze experiment: processing', length(fit.files),' files in 2 second ...\n')
  Sys.sleep(2)
  
  df = NULL
  
  missing.files = NULL
  
  for (f in 1:length(data.files))
  {
    cat('\n * ', data.files[f], '\n')
    load(paste(dir.data, '/', data.files[f], sep = ''))
    
    fitfile = paste(dir.output, '/', file.output, '.', descr, '.', data.files[f], sep = '') 
    if(!file.exists(fitfile)){
      missing.files = c(missing.files, fitfile)
      next;
    }
    load(fitfile)
    
    stats = clusters_precision(data, fit, pi.cutoff, p.collapse, p.range)
    
    # print(stats)
    
    ## TODO -- nella genetazione poopGen mettere purity in data
    
    if(is.null(data$purity))
    {
      print(data.files[f])
      foo = substr(data.files[f], 1, nchar(data.files[f]) - 6)
      foo = strsplit(foo, '_')[[1]][4]
      foo = substr(foo, 2, nchar(foo))
      data$purity = as.numeric(foo)
    }

    df = rbind(
      df, 
      data.frame(
        coverage = data$coverage, 
        K = data$K, 
        purity = data$purity,
        tp = stats$ntp,
        fp = stats$nfp,
        fn = stats$nfn,
        K_ratio = stats$K_true/stats$K_fit
      ))
  }
  
  df = df[order(df$coverage), ]
  df = df[order(df$K), ]
  df = df[order(df$purity), ]
  
  # df$K = as.character(df$K)    
  # df$coverage = as.character(df$coverage)    
  # df$purity = as.character(df$purity)    
  df$K = as.factor(df$K)    
  df$coverage = as.factor(df$coverage)    
  df$purity = as.factor(df$purity)    
  
  rownames(df) = NULL
  
  cat('\n\n** Sample size: ', nrow(df), '\n')
  print(table(df[, c('coverage', 'K')]))
  
  
  return(df)
}

analyze_alpha_posterior = function(
  dir.output = '.',
  dir.data = '.',
  file.output = 'vb_bmm1D_fit', 
  descr = '', 
  NUM_TRIALS_PER_CONFIGURATION, 
  COVERAGE_SPAN, 
  NUM_OF_CLONES_SPAN
)
{
  COVERAGE_SPAN = as.character(COVERAGE_SPAN)
  NUM_OF_CLONES_SPAN = as.character(NUM_OF_CLONES_SPAN)
  
  numtests =  length(COVERAGE_SPAN) * length(NUM_OF_CLONES_SPAN) * NUM_TRIALS_PER_CONFIGURATION
  cat('** Analyze experiment: processing', numtests,' files in 2 second ...\n')
  Sys.sleep(2)
  
  df_alpha  = NULL
  
  for (C in COVERAGE_SPAN)
  {
    for (K in NUM_OF_CLONES_SPAN)
    {
      for (i in 1:NUM_TRIALS_PER_CONFIGURATION)
      {
        fit.fname = paste(dir.output, '/', file.output, '.', descr,'.EXP-C', C, '-K', K, '-TRIAL',i, '.RData', sep = '')
        cat('\n * ', fit.fname, '\n')
        load(fit.fname)
        
        stats = fit$alpha
        
        df_alpha = rbind(df_alpha, data.frame(coverage = C, K = K, alpha = fit$alpha))
      }
    }
  }
  
  return(df_alpha)
}

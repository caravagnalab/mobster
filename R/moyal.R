dmoyal <-  function(x, loc, scale, log = FALSE) {
  
  standardized_value = (x - loc) / scale
  std_lk = standard_moyal(standardized_value)
  
  lk = std_lk - log(scale)
  
  if(!log) lk = exp(lk)
  
  return(lk)
  
}

dtruncmoyal <- function(x, loc, scale, lower, upper, log = FALSE) {
  
  moyal_lk = dmoyal(x, loc, scale, FALSE)
  moyal_lk[x < lower] = 0
  moyal_lk[x > upper] = 0
  
  norm = integrate(dmoyal, lower, upper,loc = loc, scale = scale)
  truncmoyal_lk = moyal_lk / norm$value
  
  if(log) truncmoyal_lk = log(truncmoyal_lk)
  
  return(truncmoyal_lk)
}

standard_moyal <-  function(x) {
  
  norm_const = -0.5 * log(2 * pi)
  exponent = -0.5 * (x + exp(-x))
  
  return(exponent + norm_const)
}


rtruncmoyal <- function(n, loc, scale, lower, upper) {

  smpl = -1
  while(any(smpl > upper | smpl < lower)){
    rand = runif(n = n)
    smpl = moyalppf(rand,loc, scale)
  }
  return(smpl)
}

moyalppf <- function(value, loc, scale){
  inv_err = -qnorm(0,1,p = value * 0.5) * sqrt(1/2)
  return(-log(2 * inv_err ** 2) * scale + loc)
}
  
 

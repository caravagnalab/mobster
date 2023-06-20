

data {

    int<lower=0> m;
    real<lower=0> Ntot;
    real<lower=0,upper = 1> ccf;
    real<lower=0> length_dipl_genome;
    real<lower=0> mu;
    real<lower=0> u;
    real<lower=0> sigma;
}

parameters {
 
    real<lower=0> t_mrca;
    real<lower=0,upper=t_mrca> t_driver;
    real<lower=0> s;
    
 }
 
model{

// priors
  
  t_mrca ~ uniform(0,10^12/log(2));
  t_driver ~ uniform(0,t_mrca);
  s ~ lognormal(u,sigma);
  
//  Likelihood mutations

  m ~ poisson(2*mu*log(2)*length_dipl_genome*((t_driver) + (1+s)*(t_mrca-t_driver)));

// Likelihood branching process

 target += exponential_lpdf(Ntot*ccf|exp(-log(2)*(1+s)*(log(Ntot*(1-ccf))/log(2) - t_driver)));
  
}




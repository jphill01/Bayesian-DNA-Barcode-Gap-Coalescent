// A Bayesian implementation of the DNA barcode gap coalescent
// Computes proportional probability of overlap/separation between intraspecfic and interspexific genetic distance distrubitions for a species of interest

data {
  int<lower = 0> N; // number of genetic distances
  vector<lower = 0, upper = 1>[N] intra; // intraspecific (within-species) genetic distances
  vector<lower = 0, upper = 1>[N] inter; // interspecific (among-species) genetic distances
}

transformed data {
  int<lower = 0, upper = N> y_lwr, y_upr;
  real<lower = 0, upper = 1> min_inter, max_intra;
  
  y_lwr = 0;
  y_upr = 0;
  min_inter = min(inter);
  max_intra = max(intra);
  
  for (n in 1:N) {
    y_lwr += intra[n] >= min_inter;
    y_upr += inter[n] <= max_intra;
  }
}

parameters {
  real<lower = 0, upper = 1> p_lwr; 
  real<lower = 0, upper = 1> p_upr; 
}

model {
  y_lwr ~ binomial(N, p_lwr);
  y_upr ~ binomial(N, p_upr);
  
  // priors uaing Laplace's rule of succession - updating from a beta(1, 1) = U(0,1)
  p_lwr ~ beta(y_lwr + 1, N - y_lwr + 1);
  p_upr ~ beta(y_upr + 1, N - y_upr + 1);
}


// A Bayesian model of the DNA barcode gap coalescent

// When p_lwr is close to 0, it suggests that the probability of intraspecific distances being larger than interspecific distances is low on averqge, 
// while the probability of interspecific distances being larger than intraspecific distances is high on average; that is, there is evidence for a DNA barcode gap.

// When p_upr is close to 1, it indicates that the probability of intraspecific distances being larger than interspecific distances is high on average, 
// while the probability of interspecific distances being larger than intraspecific distances is low on average; there is no evidence for a DNA barcode gap.

// If min_inter is relatively large and max_intra is relatively small, then p_lwr represents the extent to which intraspecific distances 
// tend to be larger than interspecific distances at and beyond the minimum interspecific distance and at and below the maximum intraspecific distance.

// if max_intra is relatively large and min_inter is relatively small, p_upr represents the extent to which interspecific distances 
// tend to be larger than intraspecific distances below the maximum intraspecific distance and beyond the minimum interspecific distance.

data {
  int<lower = 0> N; // number of genetic distances
  vector<lower = 0, upper = 1>[N] intra; // intraspecific (within-species) genetic distances
  vector<lower = 0, upper = 1>[N] inter; // interspecific (among-species) genetic distances
}

transformed data {
  int<lower = 0, upper = N> y_lwr; // count of intraspecific distances equalling or exceeding min_inter
  int<lower = 0, upper = N> y_upr; // count of interspecific distances less than or equal to max_intra
  real<lower = 0, upper = 1> min_inter = min(inter); // minimum interspecific distance
  real<lower = 0, upper = 1> max_intra = max(intra); // maximum intraspecific distance

  y_lwr = 0;
  y_upr = 0;

  for (n in 1:N) {
    y_lwr += (intra[n] >= min_inter); // count intraspecific distances equalling or exceeding min_inter
    y_upr += (inter[n] <= max_intra); // count interspecific distances less than or equal to max_intra
  }
 
}

parameters {
  real<lower = 0, upper = 1> p_lwr; // parameter representing the proprtional overlap/separation between intraspecific and interspecific distances
  real<lower = 0, upper = 1> p_upr; // parameter representing the proportional overlap/separation between interspecific and intraspecific distances
}

model {
  // beta(1, 1) = U(0, 1) is conjugate for binomial(n, p)
  y_lwr ~ binomial(N, p_lwr); // likelihood for intraspecific distances equalling or exceeding min_inter
  y_upr ~ binomial(N, p_upr); // likelihood for interspecific distances equalling or falling below max_intra
}

generated quantities {
  real log10_p_lwr; // log10 of p_lwr
  real log10_p_upr; // log10 of p_upr
  int<lower = 0, upper = N> y_lwr_rep[N]; // replicates of counts for y_lwr
  int<lower = 0, upper = N> y_upr_rep[N]; // replicates of counts for y_upr
  int<lower = 0, upper = 1> mean_y_lwr; // indicator variable for mean_y_lwr
  int<lower = 0, upper = 1> mean_y_upr; // indicator variable for mean_y_upr
  
  // Generate replicates of counts for y_lwr and y_upr
  for (n in 1:N) {
    y_lwr_rep[n] = binomial_rng(N, p_lwr);
    y_upr_rep[n] = binomial_rng(N, p_upr);
  }

  // Compute mean indicator variables
  mean_y_lwr = mean(y_lwr_rep) > y_lwr;
  mean_y_upr = mean(y_upr_rep) > y_upr;

  // Compute log10 of p_lwr and p_upr
  log10_p_lwr = log10(p_lwr); 
  log10_p_upr = log10(p_upr);
}


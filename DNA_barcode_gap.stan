////// A Bayesian model of the DNA barcode gap coalescent //////

// check these //

// When p_lwr is close to 0, it suggests that the probability of intraspecific distances being larger than interspecific distances is low on averqge,
// while the probability of interspecific distances being larger than intraspecific distances is high on average; that is, there is evidence for a DNA barcode gap

// When p_lwr is close to 1, it suggests that the probability of intraspecific distances being larger than interspecific distances is high on averqge,
// while the probability of interspecific distances being larger than intraspecific distances is low on average; that is, there is no evidence for a DNA barcode gap.\

// When p_upr is close to 0, it suggests that the probability of interspecific distances being larger than intraspecific distances is high on average,
// while the probability of intraspecific distances being larger than interspecific distances is low on average; that is, there is evidence for a DNA barcode gap

// When p_upr is close to 1, it suggests that the probability of interspecific distances being larger than intraspecific distances is low on average,
// while the probability of intraspecific distances being larger than interspecific distances is high on average; that is, there is no evidence for a DNA barcode gap


// When p_lwr_prime is close to 0, it suggests that the probability of intraspecific distances being larger than combined interspecific distances for a target species and its nearest neighbour species is low on averqge,
// while the probability of combined interspecific distances for a target species and its nearest neighbour species being larger than intraspecific distances is high on average; that is, there is evidence for a DNA barcode gap

// When p_lwr_prime is close to 1, it suggests that the probability of intraspecific distances being larger than combined interspecific distances for a target species and its nearest neighbour species is high on averqge,
// while the probability of combined interspecific distances for a target species and its nearest neighbour species being larger than intraspecific distances is low on average; that is, there is no evidence for a DNA barcode gap

// When p_upr_prime is close to 0, it suggests that the probability of combined interspecific distances for a target species and its nearest neighbour species being larger than intraspecific distances is high on averqge,
// while the probability of intraspecific distances being larger than combined interspecific distances for a target species and its nearest neighbour species is low on average; that is, there is evidence for a DNA barcode gap

// When p_upr_prime is close to 1, it indicates that the probability of combined interspecific distances for a target species and its nearest neighbour species being larger than intraspecific distances is low on average,
// while the probability of intraspecific distances being larger than combined interspecific distances for a target species and its nearest neighbour species is high on average; that is, there is no evidence for a DNA barcode gap


// If max_intra is relatively large and min_inter is relatively small, p_lwr represents the extent to which intraspecific distances tend to be larger than interspecific distances at and beyond min_inter and at and below max_intra
// If min_inter is relatively large and max_intra is relatively small, p_upr represents the extent to which interspecific distances tend to be larger than intraspecific distances at and below max_intra and at and beyond min_inter

// If max_intra is relatively small and min_comb is relatively large, p_lwr_prime represents the extent to which intraspecific distances tend to be larger than combined interspecific distances for a target species and its nearest neighbour species at and beyond min_comb and at and below max_intra
// If min_comb is relatively small and max_inter is relatively large, p_upr_prime represents the extent to which combined interspecific distances for a target species and its nearest neighbour species tend to be larger than intraspecific distances at and below max_intra and at and beyond min_comb


data {
  int<lower = 1> K; // number of species in genus
  int<lower = 1> N[K]; // number of intraspecific (within-species) genetic distances for each species
  int<lower = 1> M; // number of interspecific (among-species) genetic distances for all species
  int<lower = 1> C[K]; // number of combined interspecific distances for a target species and its nearest neighbour species
  vector<lower = 0, upper = 1>[sum(N)] intra; // intraspecific genetic distances for each species
  vector<lower = 0, upper = 1>[M] inter; // interspecific genetic distances for all species
  vector<lower = 0, upper = 1>[sum(C)] comb; // interspecific genetic distances for a target species and its nearest neighbour species
}


transformed data {
  int start_n[K + 1];
  int start_c[K + 1];
  
  start_n[1] = 1;
  start_c[1] = 1;
  
  for (k in 2:(K + 1)) {
    start_n[k] = start_n[k - 1] + N[k - 1];
    start_c[k] = start_c[k - 1] + C[k - 1];
  }
  
  real<lower = 0, upper = 1> min_inter = min(inter); // minimum interspecific genetic distance for all species
  vector<lower = 0, upper = 1>[K] max_intra; // maximum intraspecific genetic distance for each species
  vector<lower = 0, upper = 1>[K] min_comb; // minimum combined interspecific genetic distance for a target species and its nearest neighbour species
  
  for (k in 1:K) {
    max_intra[k] = max(segment(intra, start_n[k], N[k]));
    min_comb[k] = min(segment(comb, start_c[k], C[k]));
  }
  
  int<lower = 0, upper = max(N)> y_lwr[K] = rep_array(0, K); // count of intraspecific genetic distances for each species equalling or exceeding min_inter for all species
  int<lower = 0, upper = max(N)> y_lwr_prime[K] = rep_array(0, K); // count of intraspecific genetic distances for each species equalling or exceeding the minimum combined interspecific distance for a target species and its nearest neighbour species
  int<lower = 0, upper = M> y_upr[K] = rep_array(0, K); // count of interspecific genetic distances for all species less than or equal to max_intra for each species
  int<lower = 0, upper = max(C)> y_upr_prime[K] = rep_array(0, K); // count of combined interspecific genetic distances for a target species and its nearest neighbour species less than or equal to max_intra for each species

  for (k in 1:K) {
    for (n in 1:N[k]) {
      y_lwr[k] += (intra[start_n[k] + n - 1] >= min_inter);
      y_lwr_prime[k] += (intra[start_n[k] + n - 1] >= min_comb[k]);
    }
    
    for (m in 1:M) {
      y_upr[k] += (inter[m] <= max_intra[k]);
    }
    
    for (c in 1:C[k]) {
      y_upr_prime[k] += (comb[start_c[k] + c - 1] <= max_intra[k]);
    }
  }
  
}


parameters {
  vector<lower = 0, upper = 1>[K] p_lwr; // parameter representing the proportional overlap/separation between intraspecific genetic distances for each species and interspecific distances for all species
  vector<lower = 0, upper = 1>[K] p_upr; // parameter representing the proportional overlap/separation between interspecific genetic distances for all species and intraspecific distances for each species

  vector<lower = 0, upper = 1>[K] p_lwr_prime; // parameter representing the proportional overlap/separation between intraspecific genetic distances for each species and combined interspecific distances for a target species and its nearest neighbour species
  vector<lower = 0, upper = 1>[K] p_upr_prime; // parameter representing the proportional overlap/separation between intraspecific and intraspecific genetic distances for a target species and its nearest neighbour species

}


model {
    // Priors //
    
    // p_lwr ~ uniform(0, 1);
    // p_upr ~ uniform(0, 1);
    // p_lwr_prime ~ uniform(0, 1);
    // p_upr_prime ~ uniform(0, 1);

    // equivalent to above

    p_lwr ~ beta(1, 1);
    p_upr ~ beta(1, 1);
    p_lwr_prime ~ beta(1, 1);
    p_upr_prime ~ beta(1, 1);

    // places greater density at extemes - may cause divergent transitions etc.

    // p_lwr ~ beta(0.5, 0.5);
    // p_upr ~ beta(0.5, 0.5);
    // p_lwr_prime ~ beta(0.5, 0.5);
    // p_upr_prime ~ beta(0.5, 0.5);
    
    // Likelihood
    y_lwr ~ binomial(N, p_lwr); // likelihood for intraspecific genetic distances equalling or exceeding min_inter
    y_upr ~ binomial(M, p_upr); // likelihood for interspecific genetic distances equalling or falling below max_intra

    y_lwr_prime ~ binomial(N, p_lwr_prime); // likelihood for intraspecific genetic distances equalling or exceeding min_comb
    y_upr_prime ~ binomial(C, p_upr_prime); // likelihood for combined interspecific genetic distances for a target species and its nearest neighbour species equalling or falling below max_intra
    
}


generated quantities {
  
  // Posterior Predictive Checks

  int<lower = 0> ppc_y_lwr[K];
  int<lower = 0> ppc_y_upr[K];
  int<lower = 0> ppc_y_lwr_prime[K];
  int<lower = 0> ppc_y_upr_prime[K];

  ppc_y_lwr = binomial_rng(N, p_lwr);
  ppc_y_upr = binomial_rng(M, p_upr);
  ppc_y_lwr_prime = binomial_rng(N, p_lwr_prime);
  ppc_y_upr_prime = binomial_rng(C, p_upr_prime);
  

}

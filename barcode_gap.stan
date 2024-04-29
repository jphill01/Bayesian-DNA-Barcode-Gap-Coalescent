// A Bayesian model of the DNA barcode gap coalescent

// When p_lwr is close to 0, it suggests that the probability of intraspecific distances being larger than interspecific distances is low on averqge,
// while the probability of interspecific distances being larger than intraspecific distances is high on average; that is, there is evidence for a DNA barcode gap.

// When p_upr is close to 1, it indicates that the probability of intraspecific distances being larger than interspecific distances is high on average,
// while the probability of interspecific distances being larger than intraspecific distances is low on average; there is no evidence for a DNA barcode gap.

// If min_inter is relatively large and max_intra is relatively small, then p_lwr represents the extent to which intraspecific distances
// tend to be larger than interspecific distances at and beyond the minimum interspecific distance and at and below the maximum intraspecific distance.

// if max_intra is relatively large and min_inter is relatively small, p_upr represents the extent to which interspecific distances
// tend to be larger than intraspecific distances below the maximum intraspecific distance and beyond the minimum interspecific distance.

// Add info about p_lwr_prime, p_upr_prime, comb, and min(comb)


data {
  int<lower = 1> K; // number of species in genus
  int<lower = 1> N[K]; // number of intraspecific (within-species) genetic distances for each species
  int<lower = 1> M; // number of interspecific (among-species) genetic distances for all species
  int<lower=1> C[K]; // number of combined interspecfic distances for a target species and its nearest neighbour species
  array[K] row_vector<lower = 0, upper = 1>[max(N)] intra; // intraspecific genetic distances for each species
  vector<lower = 0, upper = 1>[M] inter; // interspecific genetic distances for all species
  array[K] row_vector<lower = 0, upper = 1>[max(C)] comb; // interspecific genetic distances for a target species and its nearest neighbour species
}

transformed data {
  int<lower = 0, upper = max(N)> y_lwr[K]; // count of intraspecific distances for each species equalling or exceeding min_inter for all species
  int<lower = 0, upper = M> y_upr[K]; // count of interspecific distances for all species less than or equal to max_intra for each species

  int<lower = 0, upper = max(N)> y_lwr_prime[K]; // count of intraspecific distances for each species equalling or exceeding the minimum combined interspecfic distance for a target species and its nearest neighbour
  int<lower = 0, upper = max(C)> y_upr_prime[K]; // count of interspecific distances for a target species and its nearest neighbour species less than or equal to max_intra for each species

  real<lower=0, upper=1> min_inter; // minimum interspecific distance for all species
  real<lower=0, upper=1> max_intra[K]; // maximum intraspecific distance for each species
  real<lower=0, upper=1> min_comb; // minimum combined distance for all species

  min_inter = min(inter);

  for (k in 1:K) {
    max_intra[k] = max(intra[k, 1:N[k]]);
    min_comb = min(comb[k, 1:C[k]]);
       
    y_lwr[k] = 0;
    y_upr[k] = 0;

    y_lwr_prime[k] = 0;
    y_upr_prime[k] = 0;

    for (n in 1:N[k]) {
      y_lwr[k] += (intra[k, n] >= min_inter); // count intraspecific distances for each species equalling or exceeding min_inter for all species
      y_lwr_prime[k] += (intra[k, n] >= min_comb); // count intraspecific distances for each species equalling or exceeding min_comb for all species
    }
    
    for (m in 1:M) {
      y_upr[k] += (inter[k] <= max_intra[k]); // count interspecific distances for all species less than or equal to max_intra for each species
    }
    
    for (c in 1:C[k]) {
      y_upr_prime[k] += (comb[k, c] <= max_intra[k]); // count combined interspecfic distances for a target species and its nearest neighbour species less than or equal to max_intra for each species
     }
     
  }

}

parameters {
  real<lower = 0, upper = 1> p_lwr[K]; // parameter representing the proportional overlap/separation between intraspecific distances for each species and interspecific distances for all species
  real<lower = 0, upper = 1> p_upr[K]; // parameter representing the proportional overlap/separation between interspecific distances for all species and intraspecific distances
  
  real<lower = 0, upper = 1> p_lwr_prime[K]; // parameter representing the proportional overlap/separation between intraspecific distances for each species and combined interspecfic distances for a target species and its nearest neighbour species
  real<lower = 0, upper = 1> p_upr_prime[K]; // parameter representing the proportional overlap/separation between intraspecific and intraspecific distances for a target species and is nearest neighbour species

}

model {
  // beta(1, 1) = U(0, 1) prior is conjugate for binomial(n, p), so posterior is beta(y_lwr + 1, N - y_lwr + 1) and beta(y_upr + 1, N - y_upr + 1) for y_lwr and y_upr, respectively
  
  for (k in 1:K) {
    // likelihood
    y_lwr[k] ~ binomial(N[k], p_lwr[k]); // likelihood for intraspecific distances equalling or exceeding min_inter
    y_upr[k] ~ binomial(M, p_upr[k]); // likelihood for interspecific distances equalling or falling below max_intra
    
    y_lwr_prime[k] ~ binomial(N[k], p_lwr_prime[k]); // likelihood for intraspecific distances equalling or exceeding min_comb
    y_upr_prime[k] ~ binomial(C[k], p_upr_prime[k]); // likelihood for combined distances equalling or falling below max_intra
  }
}


 ////// Fix this //////

generated quantities {
  real log10_p_lwr[K]; // log10 of p_lwr
  real log10_p_upr[K]; // log10 of p_upr
  real log10_p_lwr_prime[K]; // log10 of p_lwr
  real log10_p_upr_prime[K]; // log10 of p_upr
//   int<lower = 0, upper = N> y_lwr_rep[N]; // replicates of counts for y_lwr
//   int<lower = 0, upper = N> y_upr_rep[N]; // replicates of counts for y_upr
//   int<lower = 0, upper = 1> mean_y_lwr; // indicator variable for mean_y_lwr
//   int<lower = 0, upper = 1> mean_y_upr; // indicator variable for mean_y_upr
//
//   // Generate replicates of counts for y_lwr and y_upr
//   for (n in 1:N) {
//     y_lwr_rep[n] = binomial_rng(N, p_lwr);
//     y_upr_rep[n] = binomial_rng(N, p_upr);
//   }
//
//   // Compute mean indicator variables
//   mean_y_lwr = mean(y_lwr_rep) > y_lwr;
//   mean_y_upr = mean(y_upr_rep) > y_upr;
//
  // Compute log10 of p_lwr and p_upr
  log10_p_lwr[K] = log10(p_lwr[K]);
  log10_p_upr[K] = log10(p_upr[K]);
  log10_p_lwr_prime[K] = log10(p_lwr_prime[K]);
  log10_p_upr_prime[K] = log10(p_upr_prime[K]);
}





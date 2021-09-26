
data {
  int<lower=0> n_i;
  int<lower=0> n_k;
  int y[n_i,n_k];
}


parameters {
  vector[n_i] alphas;
  vector[n_k] betas;
  real<lower=0, upper=1> inv_omegas[n_k];
  real<lower=0> sigma_alpha;
  real mu_beta;
  real<lower=0> sigma_beta;
}

model {
  alphas ~ normal(0, sigma_alpha);
  betas ~ normal(mu_beta, sigma_beta);
  
  for(k in 1:n_k) {
    real omega_k_m1;
    omega_k_m1 = inv(inv(inv_omegas[k]) - 1) ;
    for (i in 1:n_i) {
      real xi_i_k;
      xi_i_k = omega_k_m1 * exp(alphas[i] + betas[k]) ;
      y[i,k] ~ neg_binomial(xi_i_k, omega_k_m1);
    }
  }
}

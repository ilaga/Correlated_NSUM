
data {
  int<lower=0> n_i;
  int<lower=0> n_k;
  int<lower=0> x_size;
  matrix[n_i, x_size] x_cov;
  matrix[n_i, n_k] z_cov;
  vector[n_i] age;
  vector[n_i] age_sq;
  int y[n_i,n_k];
}

parameters {
  vector[n_i] logdi;
  real<lower=0> sigma_di;
  matrix[n_i,n_k] eps;
  vector[x_size] beta_par;
  row_vector[n_k] age_par;
  row_vector[n_k] age_sq_par;
  vector<lower=0>[n_k] tau_N;
  cholesky_factor_corr[n_k] L_Omega;
  vector[n_k] prevalence;
  vector[n_k] alpha_par;
  real mu_prev;
  real<lower=0> sigma_prev;
}

transformed parameters {
  vector[n_k] mu;
  vector[n_k] tau;
  matrix[n_i,n_k] bias;
  matrix[n_i,n_k] prev_mean;
  
  prev_mean = exp(rep_matrix(prevalence, n_i)' + age * age_par + age_sq * age_sq_par + z_cov .* rep_matrix(alpha_par, n_i)' + rep_matrix(sigma_di * logdi + x_cov * beta_par, n_k));
  
  mu = log(1.0 ./ sqrt(1 + square(tau_N)));
  tau = sqrt(log(1 + square(tau_N)));
  bias = exp(rep_matrix(mu, n_i)' + (diag_pre_multiply(tau, L_Omega) * eps')');
}

model {
  logdi ~ normal(0, 1);
  sigma_di ~ cauchy(0, 2.5);
  tau_N ~ cauchy(0, 2.5); // Half-cauchy suggested in stan-users-guide
  to_vector(eps) ~ std_normal();
 
  prevalence ~ normal(mu_prev, sigma_prev);
  
  for(k in 1:n_k){
    y[,k] ~ poisson(prev_mean[,k] .* bias[,k]);
  }
}

generated quantities{
  matrix[n_k, n_k] Corr;
  Corr = L_Omega * L_Omega';
}



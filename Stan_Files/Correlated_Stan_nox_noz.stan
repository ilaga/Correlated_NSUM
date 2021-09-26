


data {
  int<lower=0> n_i;
  int<lower=0> n_k;
  int y[n_i,n_k];
}


parameters {
  vector[n_i] logdi;
  real<lower=0> sigma_di;
  matrix[n_i,n_k] eps;
  vector<lower=0>[n_k] tau_N;
  cholesky_factor_corr[n_k] L_Omega;
  vector[n_k] prevalence;
  real mu_prev;
  real<lower=0> sigma_prev;
}

transformed parameters {
  vector[n_k] mu;
  vector[n_k] tau;
  matrix[n_i,n_k] bias;
  vector[n_i] degree;
  
  degree = exp(sigma_di * logdi);
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
    y[,k] ~ poisson(degree * exp(prevalence[k]) .* bias[,k]);
  }
}

generated quantities{
  matrix[n_k, n_k] Corr;
  Corr = L_Omega * L_Omega';
}




data {
  int<lower=0> n_i;
  int<lower=0> n_k;
  array[n_i,n_k] int y;
}

parameters {
  vector[n_i] delta;
  real<lower=0> sigma_delta;
  matrix[n_i,n_k] eps;
  vector<lower=0>[n_k] tau_N;
  cholesky_factor_corr[n_k] L_Omega;
  vector[n_k] rho;
  real mu_rho;
  real<lower=0> sigma_rho;
}

transformed parameters {
  vector[n_k] mu;
  vector[n_k] tau;
  matrix[n_i,n_k] bias;
  matrix[n_i,n_k] prev_mean;

  prev_mean = exp(rep_matrix(rho, n_i)' + rep_matrix(sigma_delta * delta, n_k));

  mu = log(1.0 ./ sqrt(1.0 + square(tau_N)));
  tau = sqrt(log(1.0 + square(tau_N)));
  bias = exp(rep_matrix(mu, n_i)' + (diag_pre_multiply(tau, L_Omega) * eps')');
}

model {
  delta ~ normal(0, 1);
  sigma_delta ~ cauchy(0, 2.5);
  tau_N ~ cauchy(0, 2.5); // Half-cauchy suggested in stan-users-guide
  to_vector(eps) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  rho ~ normal(mu_rho, sigma_rho);

  for(k in 1:n_k){
    y[,k] ~ poisson(prev_mean[,k] .* bias[,k]);
  }
}

generated quantities{
  matrix[n_k, n_k] Corr;
  Corr = L_Omega * L_Omega';
}



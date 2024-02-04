
data {
  int<lower=0> n_i;
  int<lower=0> n_k;
  int<lower=0> z_global_size;
  int<lower=0> z_subpop_size;
  matrix[n_i, n_k] x;
  matrix[n_i, z_global_size] z_global;
  matrix[n_i, z_subpop_size] z_subpop;
  array[n_i,n_k] int y;
}

parameters {
  vector[n_i] delta;
  real<lower=0> sigma_delta;
  matrix[n_i,n_k] eps;
  vector[n_k] alpha;
  vector[z_global_size] beta_global;
  matrix[z_subpop_size, n_k] beta_subpop;
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
  matrix[n_i,n_k] prev_mean = exp(rep_matrix(rho, n_i)' + x .* rep_matrix(alpha, n_i)' + z_subpop * beta_subpop + rep_matrix(sigma_delta * delta + z_global * beta_global, n_k));

  mu = log(1.0 ./ sqrt(1.0 + square(tau_N)));
  tau = sqrt(log(1.0 + square(tau_N)));
  bias = exp(rep_matrix(mu, n_i)' + (diag_pre_multiply(tau, L_Omega) * eps')');
}

model {
  delta ~ std_normal();
  sigma_delta ~ cauchy(0, 2.5); // Half-cauchy suggested in stan-users-guide
  tau_N ~ cauchy(0, 2.5); // Half-cauchy suggested in stan-users-guide
  to_vector(eps) ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
  alpha ~ normal(0, 10);
  to_vector(beta_subpop) ~ normal(0, 10);
  beta_global ~ normal(0, 10);
  sigma_rho ~ cauchy(0, 2.5); // Half-cauchy suggested in stan-users-guide
  mu_rho ~ normal(0, 10);
  rho ~ normal(mu_rho, sigma_rho);

  for(k in 1:n_k){
    y[,k] ~ poisson(prev_mean[,k] .* bias[,k]);
  }
}

generated quantities{
  matrix[n_k, n_k] Corr;
  Corr = L_Omega * L_Omega';
}



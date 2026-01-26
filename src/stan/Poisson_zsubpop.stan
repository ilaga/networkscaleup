data {
  int<lower=0> n_i;
  int<lower=0> n_k;
  int<lower=0> z_global_size;
  matrix[n_i, z_global_size] z_global;
  array[n_i,n_k] int y;
}
parameters {
  vector[n_i] alphas;
  vector[n_k] betas;
  vector[z_global_size] beta_global;
  real<lower=0> sigma_alpha;
  real mu_beta;
  real<lower=0> sigma_beta;
}
transformed parameters {
  matrix<lower=0>[n_i,n_k] mu;
  for(i in 1:n_i){
	for(k in 1:n_k){
		mu[i,k] = exp(alphas[i] + betas[k] + z_global[i,] * beta_global);
	}
  }
}
model {
  alphas ~ normal(0, sigma_alpha);
  betas ~ normal(mu_beta, sigma_beta);
  beta_global ~ normal(0, 10);
  for(k in 1:n_k) {
    for (i in 1:n_i) {
      y[i,k] ~ poisson(mu[i,k]);
    }
  }
} 
generated quantities {
  matrix[n_i,n_k] log_lik;
  for(k in 1:n_k) {
    for (i in 1:n_i) {
      log_lik[i,k] = poisson_lpmf(y[i,k] | mu[i,k]);
    }
  }
} 

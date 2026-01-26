data {
  int<lower=0> n_i;
  int<lower=0> n_k;
  int<lower=0> z_subpop_size;
  matrix[n_i, z_subpop_size] z_subpop;
  array[n_i,n_k] int y;
}
parameters {
  vector[n_i] alphas;
  vector[n_k] betas;
  matrix[z_subpop_size, n_k] beta_subpop;
  vector<lower=0, upper=1>[n_k] inv_omegas;
  real<lower=0> sigma_alpha;
  real mu_beta;
  real<lower=0> sigma_beta;
}
transformed parameters {
  vector<lower=0>[n_k] omegas = 1.0 ./ inv_omegas;
  matrix<lower=0>[n_i,n_k] par1;
  vector<lower=0>[n_k] par2;
  for(i in 1:n_i){
    for(k in 1:n_k){
    	par1[i,k] = exp(alphas[i] + betas[k] + z_subpop[i,] * beta_subpop[,k]) / (omegas[k] - 1.0);
    }
  }
  for(k in 1:n_k){
	  par2[k] = 1.0 / (omegas[k] - 1.0);
  }
}
model {
  alphas ~ normal(0, sigma_alpha);
  betas ~ normal(mu_beta, sigma_beta);
  to_vector(beta_subpop) ~ normal(0, 10);
  for(k in 1:n_k) {
    for (i in 1:n_i) {
      y[i,k] ~ neg_binomial(par1[i,k], par2[k]);
    }
  }
} 
generated quantities {
  matrix[n_i,n_k] log_lik;
  for(k in 1:n_k) {
    for (i in 1:n_i) {
      log_lik[i,k] = neg_binomial_lpmf(y[i,k] | par1[i,k], par2[k]);
    }
  }
} 

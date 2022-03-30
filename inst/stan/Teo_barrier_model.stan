
data {
  int<lower=0> N;
  int<lower=0> n_i;
  int<lower=0> n_k;
  int<lower=0> n_known;
  int<lower=0> n_unknown;
  int<lower=0> z_size;
  matrix[n_i, n_unknown] x;
  matrix[n_i, z_size] z;
  vector<lower=0, upper=N>[n_known] S_K;
  int y[n_i,n_k];
}

parameters {
  real<lower=0,upper=10> lambda;
  vector<lower=0>[n_i] alpha;
  real<lower=0,upper=10> tau;
  vector[n_unknown] beta;
  matrix[z_size, n_k] gamma;
  real<lower=0> sigma;
  real<lower=0> sigmag;
  vector<lower=0,upper=N>[n_unknown] S_H;
}

model {
  alpha ~ lognormal(0, tau);
  beta ~ normal(0, 1 / sigma);
  for(k in 1:n_k){
    gamma[,k] ~ normal(0, 1 / sigmag);
  }
  sigma ~ gamma(1, 0.1);
  sigmag ~ gamma(1, 0.1);
  
  for(k in 1:n_known){
	y[,k] ~ poisson(lambda * alpha .* exp(z * gamma[,k]) * S_K[k]);
  }
  
  for(k in 1:n_unknown){
	y[,n_known + k] ~ poisson(lambda * alpha .* exp(x[,k] * beta[k]) .* exp(z * gamma[,n_known + k]) * S_H[k]);
  }
}



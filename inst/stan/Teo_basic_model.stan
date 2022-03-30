
data {
  int<lower=0> N;
  int<lower=0> n_i;
  int<lower=0> n_k;
  int<lower=0> n_known;
  int<lower=0> n_unknown;
  vector<lower=0, upper=N>[n_known] S_K;
  int y[n_i,n_k];
}

parameters {
  real<lower=0,upper=10> lambda;
  vector<lower=0>[n_i] alpha;
  real<lower=0,upper=10> tau;
  vector<lower=0,upper=N>[n_unknown] S_H;
}

model {
  alpha ~ lognormal(0, tau);
  
  for(k in 1:n_known){
	y[,k] ~ poisson(lambda * alpha * S_K[k]);
  }
  
  for(k in 1:n_unknown){
	y[,n_known + k] ~ poisson(lambda * alpha * S_H[k]);
  }
}



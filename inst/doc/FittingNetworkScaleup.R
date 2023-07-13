## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(networkscaleup)

## ----simulation---------------------------------------------------------------
set.seed(1998)
N_i = 50
N_k = 5
N = 1e5
sizes = rbinom(N_k, N, prob = runif(N_k, min = 0.01, max = 0.15))
degrees = round(exp(rnorm(N_i, mean = 5, sd = 1)))

ard = matrix(NA, nrow = N_i, ncol = N_k)
for(k in 1:N_k){
  ard[,k] = rbinom(N_i, degrees, sizes[k] / N)
}

## Create some artificial covariates for use later
x = matrix(sample(1:5, size = N_i * N_k, replace = T),
           nrow = N_i,
           ncol = N_k)
z = cbind(rbinom(N_i, 1, 0.3), rnorm(N_i), rnorm(N_i), rnorm(N_i))

## ----prep---------------------------------------------------------------------
x = scale(x)
z = scale(z)
z_subpop = z[,1:2]
z_global = z[,3:4]


## ----pimle--------------------------------------------------------------------
pimle.est = killworth(ard,
  known_sizes = sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  N = N, model = "PIMLE")


plot(degrees ~ pimle.est$degrees, xlab = "Estimated PIMLE degrees", ylab = "True Degrees")
abline(0, 1, col = "red")

round(data.frame(true = sizes[c(3, 5)],
                 pimle = pimle.est$sizes))

## ----mle----------------------------------------------------------------------
mle.est = killworth(ard,
  known_sizes = sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  N = N, model = "MLE")


plot(degrees ~ mle.est$degrees, xlab = "Estimated MLE degrees", ylab = "True Degrees")
abline(0, 1, col = "red")

round(data.frame(true = sizes[c(3, 5)],
                 pimle = mle.est$sizes))

## ----overdisp-----------------------------------------------------------------
overdisp_gibbs_metrop_est = overdispersed(
  ard,
  known_sizes = sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  G1_ind = 1,
  G2_ind = 2,
  B2_ind = 4,
  N = N,
  warmup = 500,
  iter = 1000,
  verbose = TRUE,
  init = "MLE"
)

overdisp_stan = overdispersedStan(
  ard,
  known_sizes = sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  G1_ind = 1,
  G2_ind = 2,
  B2_ind = 4,
  N = N,
  chains = 2,
  cores = 2,
  warmup = 250,
  iter = 500,
)


round(data.frame(true = sizes,
                 gibbs_est = colMeans(overdisp_gibbs_metrop_est$sizes),
                 stan_est = colMeans(overdisp_stan$sizes)))

plot(degrees ~ colMeans(overdisp_stan$degrees), xlab = "Overdispersed Degree Estimates", ylab = "True Degrees")
abline(0, 1, col = "red")

## ----correlated---------------------------------------------------------------
correlated_cov_stan = correlatedStan(
  ard,
  known_sizes = sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  model = "correlated",
  scaling = "weighted",
  x = x,
  z_subpop = z_subpop,
  z_global = z_global,
  N = N,
  chains = 2,
  cores = 2,
  warmup = 250,
  iter = 500,
)

correlated_nocov_stan = correlatedStan(
  ard,
  known_sizes = sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  model = "correlated",
  scaling = "all",
  N = N,
  chains = 2,
  cores = 2,
  warmup = 250,
  iter = 500,
)


uncorrelated_cov_stan = correlatedStan(
  ard,
  known_sizes = sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  model = "uncorrelated",
  scaling = "all",
  x = x,
  z_subpop = z_subpop,
  z_global = z_global,
  N = N,
  chains = 2,
  cores = 2,
  warmup = 250,
  iter = 500,
)

uncorrelated_x_stan = correlatedStan(
  ard,
  known_sizes = sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  model = "uncorrelated",
  scaling = "all",
  x = x,
  N = N,
  chains = 2,
  cores = 2,
  warmup = 250,
  iter = 500,
)



round(data.frame(true = sizes,
                 corr_cov_est = colMeans(correlated_cov_stan$sizes),
                 corr_nocov_est = colMeans(correlated_nocov_stan$sizes),
                 uncorr_cov_est = colMeans(uncorrelated_cov_stan$sizes),
                 uncorr_x_est = colMeans(uncorrelated_x_stan$sizes)))

plot(degrees ~ colMeans(correlated_cov_stan$degrees), xlab = "Correlated Covariate Degree Estimates", ylab = "True Degrees")
abline(0, 1, col = "red")

## Examine parameter estimates
colMeans(correlated_cov_stan$alpha)
colMeans(correlated_cov_stan$beta_global)
colMeans(correlated_cov_stan$beta_subpop)


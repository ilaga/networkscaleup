# Fit basic Poisson and Negative Binomial models using glmmTMB

Fit basic Poisson and Negative Binomial models using glmmTMB

## Usage

``` r
fit_mle(
  ard,
  x_cov_global = NULL,
  x_cov_local = NULL,
  family = c("poisson", "nbinomial")
)
```

## Arguments

- ard:

  n_i by n_k ARD matrix

- x_cov_global:

  n_i by p_global covariate matrix of global covariates

- x_cov_local:

  n_i by p_local covariate matrix of local covariates

- family:

  distribution to fit, either "poisson" or "nbinomial"

## Value

list containing fitted model and extracted parameters

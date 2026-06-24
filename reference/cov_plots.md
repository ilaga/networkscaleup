# Covariance plots

Plots of the estimated covariance structure from a given fitted model

## Usage

``` r
cov_plots(
  ard,
  model_fit,
  x_cov,
  resid_type = c("rqr", "pearson_residuals"),
  method = "lm",
  se = F
)
```

## Arguments

- ard:

  ard matrix

- model_fit:

  a fitted object from \[fit_mle()\] or \[fit_map()\]

- x_cov:

  covariate matrix

- resid_type:

  the type of residuals to use

- method:

  the method to use

- se:

  whether to compute standard errors of estimates

## Value

a list of ggplots, corresponding to covariance structure

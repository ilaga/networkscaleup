# Plot residuals against fitted values

Plot residuals against fitted values

## Usage

``` r
plot_fitted(ard, model_fit = NULL, resid = c("rqr", "pearson", "surrogate"))
```

## Arguments

- ard:

  ARD matrix (may be needed)

- model_fit:

  fitted model

- resid:

  the type of residuals to be used

## Value

a ggplot showing fitted values against residuals

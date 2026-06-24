# Tracy-Widom test for residual group correlation

Tracy-Widom test for residual group correlation

## Usage

``` r
tw_group_corr_test(model_fit, correction = c("none", "half"), plot = TRUE)
```

## Arguments

- model_fit:

  fitted model object

- correction:

  correction constant, either "none", "half"

- plot:

  a logical, whether to return a ggplot density plot of TW with observed
  statistic

## Value

a list containing test statistic, p-value, and diagnostic plots

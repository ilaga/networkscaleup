# Hanging Rootogram for Fitted ARD Model

Hanging Rootogram for Fitted ARD Model

## Usage

``` r
hang_rootogram_ard(ard, model_fit, width = 0.9, x_max = NULL, by_group = FALSE)
```

## Arguments

- ard:

  ard matrix

- model_fit:

  fitted model object

- width:

  width of bars

- x_max:

  the maximum x value to display

- by_group:

  logical; if TRUE, create separate rootograms for each column (group)

## Value

a ggplot of the hanging rootogram (single plot if by_group=FALSE,
combined plot if by_group=TRUE)

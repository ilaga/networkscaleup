# Scale raw log degree and log prevalence estimates

This function scales estimates from either the overdispersed model or
from the correlated models. Several scaling options are available.

## Usage

``` r
scaling(
  log_degrees,
  log_prevalences,
  scaling = c("all", "overdispersed", "weighted", "weighted_sq"),
  known_sizes = NULL,
  known_ind = NULL,
  Correlation = NULL,
  G1_ind = NULL,
  G2_ind = NULL,
  B2_ind = NULL,
  N = NULL
)
```

## Arguments

- log_degrees:

  The matrix of estimated raw log degrees from either the overdispersed
  or correlated models.

- log_prevalences:

  The matrix of estimates raw log prevalences from either the
  overdispersed or correlated models.

- scaling:

  An character vector providing the name of scaling procedure should be
  performed in order to transform estimates to degrees and subpopulation
  sizes. Scaling options are \`overdispersed\`, \`all\` (the default),
  \`weighted\`, or \`weighted_sq\` (\`weighted\` and \`weighted_sq\` are
  only available if \`Correlation\` is provided. Further details are
  provided in the Details section.

- known_sizes:

  The known subpopulation sizes corresponding to a subset of the columns
  of `ard`.

- known_ind:

  The indices that correspond to the columns of `ard` with known_sizes.
  By default, the function assumes the first `n_known` columns, where
  `n_known` corresponds to the number of `known_sizes`.

- Correlation:

  The estimated correlation matrix used to calculate scaling weights.
  Required if \`scaling = weighted\` or \`scaling = weighted_sq\`.

- G1_ind:

  If \`scaling = overdispersed\`, a vector of indices corresponding to
  the subpopulations that belong to the primary scaling groups, i.e. the
  collection of rare girls' names in Zheng, Salganik, and Gelman (2006).
  By default, all known_sizes are used. If G2_ind and B2_ind are not
  provided, \`C = C_1\`, so only G1_ind are used. If G1_ind is not
  provided, no scaling is performed.

- G2_ind:

  If \`scaling = overdispersed\`, a vector of indices corresponding to
  the subpopulations that belong to the first secondary scaling groups,
  i.e. the collection of somewhat popular girls' names.

- B2_ind:

  If \`scaling = overdispersed\`, a vector of indices corresponding to
  the subpopulations that belong to the second secondary scaling groups,
  i.e. the collection of somewhat popular boys' names.

- N:

  The known total population size.

## Value

The named list containing the scaled log degree, degree, log prevalence,
and size estimates

## Details

The \`scaling\` options are described below:

- overdispersed:

  The scaling procedure outlined in Zheng et al. (2006) is performed. In
  this case, at least \`Pg1_ind\` must be provided. See
  [overdispersedStan](https://ilaga.github.io/networkscaleup/reference/overdispersedStan.md)
  for more details.

- all:

  All subpopulations with known sizes are used to scale the parameters,
  using a modified scaling procedure that standardizes the sizes so each
  population is weighted equally. Additional details are provided in
  Laga et al. (2021).

- weighted:

  All subpopulations with known sizes are weighted according their
  correlation with the unknown subpopulation size. Additional details
  are provided in Laga et al. (2021)

- weighted_sq:

  Same as \`weighted\`, except the weights are squared, providing more
  relative weight to subpopulations with higher correlation.

## References

Zheng, T., Salganik, M. J., and Gelman, A. (2006). How many people do
you know in prison, *Journal of the American Statistical Association*,
**101:474**, 409–423

Laga, I., Bao, L., and Niu, X (2021). A Correlated Network Scaleup
Model: Finding the Connection Between Subpopulations

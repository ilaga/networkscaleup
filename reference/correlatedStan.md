# Fit ARD using the uncorrelated or correlated model in Stan This function fits the ARD using either the uncorrelated or correlated model in Laga et al. (2021) in Stan. The population size estimates and degrees are scaled using a post-hoc procedure.

Fit ARD using the uncorrelated or correlated model in Stan This function
fits the ARD using either the uncorrelated or correlated model in Laga
et al. (2021) in Stan. The population size estimates and degrees are
scaled using a post-hoc procedure.

## Usage

``` r
correlatedStan(
  ard,
  known_sizes = NULL,
  known_ind = NULL,
  N = NULL,
  model = c("correlated", "uncorrelated"),
  scaling = c("all", "overdispersed", "weighted", "weighted_sq"),
  x = NULL,
  z_global = NULL,
  z_subpop = NULL,
  G1_ind = NULL,
  G2_ind = NULL,
  B2_ind = NULL,
  chains = 3,
  cores = 1,
  warmup = 1000,
  iter = 1500,
  thin = 1,
  return_fit = FALSE,
  ...
)
```

## Arguments

- ard:

  The \`n_i x n_k\` matrix of non-negative ARD integer responses, where
  the \`(i,k)th\` element corresponds to the number of people that
  respondent \`i\` knows in subpopulation \`k\`.

- known_sizes:

  The known subpopulation sizes corresponding to a subset of the columns
  of `ard`.

- known_ind:

  The indices that correspond to the columns of `ard` with known_sizes.
  By default, the function assumes the first `n_known` columns, where
  `n_known` corresponds to the number of `known_sizes`.

- N:

  The known total population size.

- model:

  A character vector denoting which of the two models should be fit,
  either 'uncorrelated' or 'correlated'. More details of these models
  are provided below. The function decides which covariate model is
  needed based on the covariates provided below.

- scaling:

  An optional character vector providing the name of scaling procedure
  should be performed in order to transform estimates to degrees and
  subpopulation sizes. If \`NULL\`, the parameters will be returned
  unscaled. Alternatively, scaling may be performed independently using
  the
  [scaling](https://ilaga.github.io/networkscaleup/reference/scaling.md)
  function. Scaling options are \`NULL\`, \`overdispersed\`, \`all\`,
  \`weighted\`, or \`weighted_sq\` (\`weighted\` and \`weighted_sq\` are
  only available if \`model = "correlated"\`. Further details are
  provided in the Details section.

- x:

  A matrix with dimensions \`n_i x n_unknown\`, where \`n_unknown\`
  refers to the number of unknown subpopulation sizes. In the language
  of Teo et al. (2019), these represent the individual's perception of
  each hidden population.

- z_global:

  A matrix with dimensions \`n_i x p_global\`, where \`p_global\` is the
  number of demographic covariates used. This matrix represents the
  demographic information about the respondents in order to capture the
  barrier effects.

- z_subpop:

  A matrix with dimensions \`n_i x p_subpop\`, where \`p_subpop\` is the
  number of demographic covariates used. This matrix represents the
  demographic information about the respondents in order to capture the
  barrier effects.

- G1_ind:

  A vector of indices denoting the columns of \`ard\` that correspond to
  the primary scaling groups, i.e. the collection of rare girls' names
  in Zheng, Salganik, and Gelman (2006). By default, all known_sizes are
  used. If G2_ind and B2_ind are not provided, \`C = C_1\`, so only
  G1_ind are used. If G1_ind is not provided, no scaling is performed.

- G2_ind:

  A vector of indices denoting the columns of \`ard\` that correspond to
  the subpopulations that belong to the first secondary scaling groups,
  i.e. the collection of somewhat popular girls' names.

- B2_ind:

  A vector of indices denoting the columns of \`ard\` that correspond to
  the subpopulations that belong to the second secondary scaling groups,
  i.e. the collection of somewhat popular boys' names.

- chains:

  A positive integer specifying the number of Markov chains.

- cores:

  A positive integer specifying the number of cores to use to run the
  Markov chains in parallel.

- warmup:

  A positive integer specifying the total number of samples for each
  chain (including warmup). Matches the usage in
  [stan](https://mc-stan.org/rstan/reference/stan.html).

- iter:

  A positive integer specifying the number of warmup samples for each
  chain. Matches the usage in
  [stan](https://mc-stan.org/rstan/reference/stan.html).

- thin:

  A positive integer specifying the interval for saving posterior
  samples. Default value is 1 (i.e. no thinning).

- return_fit:

  A logical indicating whether the fitted \`stanfit\` object should be
  return. Defaults to \`FALSE\`.

- ...:

  Additional arguments to be passed to
  [stan](https://mc-stan.org/rstan/reference/stan.html).

## Value

Either the full fitted Stan model if `return_fit = TRUE`, else a named
list with the estimated parameters extracted using
[extract](https://mc-stan.org/rstan/reference/stanfit-method-extract.html)
(the default). The estimated parameters are named as follows (if
estimated in the corresponding model), with additional descriptions as
needed:

- sigma_delta:

  Standard deviation of delta

- rho:

  Log prevalence, if scaled, else raw rho parameters

- mu_rho:

  Mean of rho

- sigma_rho:

  Standard deviation of rho

- alpha:

  Slope parameters corresponding to z

- beta_global:

  Slope parameters corresponding to x_global

- beta_subpop:

  Slope parameters corresponding to x_subpop

- tau_N:

  Standard deviation of random effects b

- Corr:

  Correlation matrix, if \`Correlation = TRUE\`

If scaled, the following additional parameters are included:

- degree:

  Scaled degrees

- log_prevalences:

  Scaled log prevalences

- sizes:

  Subpopulation size estimates

## Details

This function currently fits a variety of models proposed in Laga et al.
(2022+). The user may provide any combination of \`x\`, \`z_global\`,
and \`z_subpop\`. Additionally, the user may choose to fit a
uncorrelated version of the model, where the correlation matrix is equal
to the identity matrix.

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
  Laga et al. (2022+).

- weighted:

  All subpopulations with known sizes are weighted according their
  correlation with the unknown subpopulation size. Additional details
  are provided in Laga et al. (2022+)

- weighted_sq:

  Same as \`weighted\`, except the weights are squared, providing more
  relative weight to subpopulations with higher correlation.

## References

Laga, I., Bao, L., and Niu, X (2021). A Correlated Network Scaleup
Model: Finding the Connection Between Subpopulations

## Examples

``` r
if (FALSE) { # \dontrun{
data(example_data)

x <- example_data$x
z_global <- example_data$z[, 1:2]
z_subpop <- example_data$z[, 3:4]

basic_corr_est <- correlatedStan(example_data$ard,
  known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  N = example_data$N,
  model = "correlated",
  scaling = "weighted",
  chains = 1,
  cores = 1,
  warmup = 50,
  iter = 100
)

cov_uncorr_est <- correlatedStan(example_data$ard,
  known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  N = example_data$N,
  model = "uncorrelated",
  scaling = "all",
  x = x,
  z_global = z_global,
  z_subpop = z_subpop,
  chains = 1,
  cores = 1,
  warmup = 50,
  iter = 100
)

cov_corr_est <- correlatedStan(example_data$ard,
  known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
  known_ind = c(1, 2, 4),
  N = example_data$N,
  model = "correlated",
  scaling = "all",
  x = x,
  z_subpop = z_subpop,
  chains = 1,
  cores = 1,
  warmup = 50,
  iter = 100
)

# Compare size estimates
round(data.frame(
  true = example_data$subpop_sizes,
  corr_basic = colMeans(basic_corr_est$sizes),
  uncorr_x_zsubpop_zglobal = colMeans(cov_uncorr_est$sizes),
  corr_x_zsubpop = colMeans(cov_corr_est$sizes)
))

# Look at z slope parameters
colMeans(cov_uncorr_est$beta_global)
colMeans(cov_corr_est$beta_subpop)
colMeans(cov_uncorr_est$beta_subpop)

# Look at x slope parameters
colMeans(cov_uncorr_est$alpha)
colMeans(cov_corr_est$alpha)
} # }
```

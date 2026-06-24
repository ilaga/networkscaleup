# Fit ARD using the Overdispersed model in Stan

This function fits the ARD using the Overdispersed model in Stan. The
population size estimates and degrees are scaled using a post-hoc
procedure. For the Gibbs-Metropolis algorithm implementation, see
[overdispersed](https://ilaga.github.io/networkscaleup/reference/overdispersed.md).

## Usage

``` r
overdispersedStan(
  ard,
  known_sizes = NULL,
  known_ind = NULL,
  G1_ind = NULL,
  G2_ind = NULL,
  B2_ind = NULL,
  N = NULL,
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

- N:

  The known total population size.

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

  A logical indicating whether the fitted Stan model should be returned
  instead of the rstan::extracted and scaled parameters. This is FALSE
  by default.

- ...:

  Additional arguments to be passed to
  [stan](https://mc-stan.org/rstan/reference/stan.html).

## Value

Either the full fitted Stan model if `return_fit = TRUE`, else a named
list with the estimated parameters extracted using
[extract](https://mc-stan.org/rstan/reference/stanfit-method-extract.html)
(the default). The estimated parameters are named as follows, with
additional descriptions as needed:

- betas:

  Log prevalence, if \`scaling = TRUE\`, else raw beta parameters

- inv_omegas:

  Inverse of overdispersion parameters

- sigma_alpha:

  Standard deviation of alphas

- mu_beta:

  Mean of betas

- sigma_beta:

  Standard deviation of betas

- omegas:

  Overdispersion parameters

If \`scaling = TRUE\`, the following additional parameters are included:

- degrees:

  Degree estimates

- sizes:

  Subpopulation size estimates

## Details

This function fits the overdispersed NSUM model using the
Gibbs-Metropolis algorithm provided in Zheng et al. (2006).

## References

Zheng, T., Salganik, M. J., and Gelman, A. (2006). How many people do
you know in prison, *Journal of the American Statistical Association*,
**101:474**, 409–423

## Examples

``` r
# Analyze an example ard data set using Zheng et al. (2006) models
# Note that in practice, both warmup and iter should be much higher
if (FALSE) { # \dontrun{
data(example_data)

ard <- example_data$ard
subpop_sizes <- example_data$subpop_sizes
known_ind <- c(1, 2, 4)
N <- example_data$N

overdisp.est <- overdispersedStan(ard,
  known_sizes = subpop_sizes[known_ind],
  known_ind = known_ind,
  G1_ind = 1,
  G2_ind = 2,
  B2_ind = 4,
  N = N,
  chains = 1,
  cores = 1,
  warmup = 250,
  iter = 500
)

# Compare size estimates
round(data.frame(
  true = subpop_sizes,
  basic = colMeans(overdisp.est$sizes)
))

# Compare degree estimates
plot(example_data$degrees, colMeans(overdisp.est$degrees))

# Look at overdispersion parameter
colMeans(overdisp.est$omegas)
} # }
```

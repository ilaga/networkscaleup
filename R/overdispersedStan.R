#' Fit ARD using the Overdispersed model in Stan
#'
#' This function fits the ARD using the Overdispersed model in Stan. The
#' population size estimates and degrees are scaled using a post-hoc procedure.
#' For the Gibbs-Metropolis algorithm implementation, see
#' \link[networkscaleup]{overdispersed}.
#'
#' @param ard The `n_i x n_k` matrix of non-negative ARD integer responses,
#'   where the `(i,k)th` element corresponds to the number of people that
#'   respondent `i` knows in subpopulation `k`.
#' @param known_sizes The known subpopulation sizes corresponding to a subset of
#'   the columns of \code{ard}.
#' @param known_ind The indices that correspond to the columns of \code{ard}
#'   with known_sizes. By default, the function assumes the first \code{n_known}
#'   columns, where \code{n_known} corresponds to the number of
#'   \code{known_sizes}.
#' @param G1_ind A vector of indices denoting the columns of `ard` that
#'   correspond to the primary scaling groups, i.e. the collection of rare
#'   girls' names in Zheng, Salganik, and Gelman (2006). By default, all
#'   known_sizes are used. If G2_ind and B2_ind are not provided, `C = C_1`, so
#'   only G1_ind are used. If G1_ind is not provided, no scaling is performed.
#' @param G2_ind A vector of indices denoting the columns of `ard` that
#'   correspond to the subpopulations that belong to the first secondary scaling
#'   groups, i.e. the collection of somewhat popular girls' names.
#' @param B2_ind A vector of indices denoting the columns of `ard` that
#'   correspond to the subpopulations that belong to the second secondary
#'   scaling groups, i.e. the collection of somewhat popular boys' names.
#' @param N The known total population size.
#' @param chains A positive integer specifying the number of Markov chains.
#' @param cores A positive integer specifying the number of cores to use to run
#'   the Markov chains in parallel.
#' @param warmup A positive integer specifying the total number of samples for
#'   each chain (including warmup). Matches the usage in \link[rstan]{stan}.
#' @param iter A positive integer specifying the number of warmup samples for
#'   each chain. Matches the usage in \link[rstan]{stan}.
#' @param thin A positive integer specifying the interval for saving posterior
#'   samples. Default value is 1 (i.e. no thinning).
#' @param return_fit A logical indicating whether the fitted Stan model should
#'   be returned instead of the rstan::extracted and scaled parameters. This is
#'   FALSE by default.
#' @param ... Additional arguments to be passed to \link[rstan]{stan}.
#'
#' @details This function fits the overdispersed NSUM model using the
#'   Gibbs-Metropolis algorithm provided in Zheng et al. (2006).
#'
#' @return Either the full fitted Stan model if \code{return_fit = TRUE}, else a
#'   named list with the estimated parameters extracted using
#'   \link[rstan]{extract} (the default). The estimated parameters are named as
#'   follows, with additional descriptions as needed:
#'
#'   \describe{\item{alphas}{Log degree, if `scaling = TRUE`, else raw alpha parameters}
#'   \item{betas}{Log prevalence, if `scaling = TRUE`, else raw beta parameters}
#'   \item{inv_omegas}{Inverse of overdispersion parameters}
#'   \item{sigma_alpha}{Standard deviation of alphas}
#'   \item{mu_beta}{Mean of betas}
#'   \item{sigma_beta}{Standard deviation of betas}
#'   \item{omegas}{Overdispersion parameters}}
#'
#'   If `scaling = TRUE`, the following additional parameters are included:
#'   \describe{\item{mu_alpha}{Mean of log degrees}
#'   \item{degrees}{Degree estimates}
#'   \item{sizes}{Subpopulation size estimates}}
#' @references Zheng, T., Salganik, M. J., and Gelman, A. (2006). How many
#'   people do you know in prison, \emph{Journal of the American Statistical
#'   Association}, \bold{101:474}, 409--423
#' @export
#'
#' @examples
#' # Analyze an example ard data set using Zheng et al. (2006) models
#' # Note that in practice, both warmup and iter should be much higher
#' \dontrun{
#' data(example_data)
#'
#' ard = example_data$ard
#' subpop_sizes = example_data$subpop_sizes
#' known_ind = c(1, 2, 4)
#' N = example_data$N
#'
#' overdisp.est = overdispersedStan(ard,
#' known_sizes = subpop_sizes[known_ind],
#' known_ind = known_ind,
#' G1_ind = 1,
#' G2_ind = 2,
#' B2_ind = 4,
#' N = N,
#' chains = 1,
#' cores = 1,
#' warmup = 250,
#' iter = 500)
#'
#' # Compare size estimates
#' round(data.frame(true = subpop_sizes,
#' basic = colMeans(overdisp.est$sizes)))
#'
#' # Compare degree estimates
#' plot(example_data$degrees, colMeans(overdisp.est$degrees))
#'
#' # Look at overdispersion parameter
#' colMeans(overdisp.est$omegas)
#' }
overdispersedStan <-
  function(ard,
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
           ...) {
    N_i = nrow(ard)
    N_k = ncol(ard)



    known_prevalences = known_sizes / N
    prevalences_vec = rep(NA, N_k)
    prevalences_vec[known_ind] = known_prevalences
    if (!is.null(G1_ind)) {
      Pg1 = sum(prevalences_vec[G1_ind])
    }
    if (!is.null(G2_ind)) {
      Pg2 = sum(prevalences_vec[G2_ind])
    }
    if (!is.null(B2_ind)) {
      Pb2 = sum(prevalences_vec[B2_ind])
    }

    stan_data = list(n_i = N_i,
                     n_k = N_k,
                     y = ard)


    ## Fit model
    overdispersed_fit = rstan::sampling(
      object = stanmodels$Overdispersed_Stan,
      data = stan_data,
      chains = chains,
      cores = cores,
      iter = iter,
      warmup = warmup,
      refresh = 100,
      thin = thin,
      ...
    )

    ## Extract variables


    draws = rstan::extract(overdispersed_fit)


    betas = draws$betas
    mu_beta = draws$mu_beta
    alphas = draws$alphas

    mu_alpha = rep(NA, nrow(betas))

    ## Perform scaling
    if (is.null(G1_ind)) {
      ## Perform no scaling

    } else if (is.null(G2_ind) |
               is.null(B2_ind)) {
      ## Perform scaling with only main

      for (ind in 1:nrow(betas)) {
        C1 = log(sum(exp(betas[ind, G1_ind]) / Pg1))
        C = C1

        alphas[ind, ] = alphas[ind, ] + C
        mu_alpha[ind] = C
        betas[ind, ] = betas[ind, ] - C
        mu_beta[ind] = mu_beta[ind] - C
      }

      draws$betas = betas
      draws$alphas = alphas
      draws$mu_beta = mu_beta
      draws$mu_alpha = mu_alpha
      draws$degrees = exp(draws$alphas)
      draws$sizes = exp(draws$betas) * N

    } else{
      ## Perform scaling with secondary groups
      for (ind in 1:nrow(betas)) {
        C1 = log(sum(exp(betas[ind, G1_ind]) / Pg1))
        C2 = log(sum(exp(betas[ind, B2_ind]) / Pb2)) - log(sum(exp(betas[ind, G2_ind]) / Pg2))
        C = C1 + 1 / 2 * C2

        alphas[ind, ] = alphas[ind, ] + C
        mu_alpha[ind] = C
        betas[ind, ] = betas[ind, ] - C
        mu_beta[ind] = mu_beta[ind] - C
      }

      draws$betas = betas
      draws$alphas = alphas
      draws$mu_beta = mu_beta
      draws$mu_alpha = mu_alpha
      draws$degrees = exp(draws$alphas)
      draws$sizes = exp(draws$betas) * N
    }

    ## Return values
    if (return_fit) {
      return(overdispersed_fit)
    } else{
      return(draws)
    }

  }

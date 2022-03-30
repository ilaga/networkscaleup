#' Fit Overdispersed model to ARD (Gibbs-Metropolis)
#'
#' This function fits the ARD using the Overdispersed model using the original
#' Gibbs-Metropolis Algorithm provided in Zheng, Salganik, and Gelman (2006).
#' The population size estimates and degrees are scaled using a post-hoc
#' procedure. For the Stan implementation, see
#' \link[networkscaleup]{overdispersedStan}.
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
#' @param G1_ind A vector of indices corresponding to the subpopulations that
#'   belong to the primary scaling groups, i.e. the collection of rare girls'
#'   names in Zheng, Salganik, and Gelman (2006). By default, all known_sizes
#'   are used. If G2_ind and B2_ind are not provided, `C = C_1`, so only G1_ind
#'   are used. If G1_ind is not provided, no scaling is performed.
#' @param G2_ind A vector of indices corresponding to the subpopulations that
#'   belong to the first secondary scaling groups, i.e. the collection of
#'   somewhat popular girls' names.
#' @param B2_ind A vector of indices corresponding to the subpopulations that
#'   belong to the second secondary scaling groups, i.e. the collection of
#'   somewhat popular boys' names.
#' @param N The known total population size.
#' @param warmup A positive integer specifying the number of warmup samples.
#' @param iter A positive integer specifying the total number of samples
#'   (including warmup).
#' @param refresh An integer specifying how often the progress of the sampling
#'   should be reported. By default, resorts to every 10%. Suppressed if
#'   \code{verbose = FALSE}.
#' @param verbose Logical value, specifying whether sampling progress should be
#'   reported.
#' @param thin A positive integer specifying the interval for saving posterior
#'   samples. Default value is 1 (i.e. no thinning).
#' @param alpha_tune A positive numeric indicating the standard deviation used
#'   as the jumping scale in the Metropolis step for alpha. Defaults to 0.4,
#'   which has worked well for other ARD datasets.
#' @param beta_tune A positive numeric indicating the standard deviation used as
#'   the jumping scale in the Metropolis step for beta Defaults to 0.2, which
#'   has worked well for other ARD datasets.
#' @param omega_tune A positive numeric indicating the standard deviation used
#'   as the jumping scale in the Metropolis step for omega Defaults to 0.2,
#'   which has worked well for other ARD datasets.
#' @param init A named list with names corresponding to the first-level model
#'   parameters, name 'alpha', 'beta', and 'omega'. By default the 'alpha' and
#'   'beta' parameters are initialized at the values corresponding to the
#'   Killworth MLE estimates (for the missing 'beta'), with all 'omega' set to
#'   20. Alternatively, \code{init = 'random'} simulates 'alpha' and 'beta' from
#'   a normal random variable with mean 0 and standard deviation 1.
#'
#' @details This function fits the overdispersed NSUM model using the
#'   Metropolis-Gibbs sampler provided in Zheng et al. (2006).
#'
#' @return A named list with the estimated posterior samples. The estimated
#' parameters are named as follows, with additional descriptions as needed:
#'
#'   \describe{\item{alphas}{Log degree, if scaled, else raw alpha parameters}
#'   \item{betas}{Log prevalence, if scaled, else raw beta parameters}
#'   \item{inv_omegas}{Inverse of overdispersion parameters}
#'   \item{sigma_alpha}{Standard deviation of alphas}
#'   \item{mu_beta}{Mean of betas}
#'   \item{sigma_beta}{Standard deviation of betas}
#'   \item{omegas}{Overdispersion parameters}}
#'
#'   If scaled, the following additional parameters are included:
#'   \describe{\item{mu_alpha}{Mean of log degrees}
#'   \item{degrees}{Degree estimates}
#'   \item{sizes}{Subpopulation size estimates}}
#' @references Zheng, T., Salganik, M. J., and Gelman, A. (2006). How many
#'   people do you know in prison, \emph{Journal of the American Statistical
#'   Association}, \bold{101:474}, 409--423
#' @import rstan
#' @importFrom LaplacesDemon rinvchisq
#' @export
#'
#' @examples
#' # Analyze an example ard data set using Zheng et al. (2006) models
#' # Note that in practice, both warmup and iter should be much higher
#' data(example_data)
#'
#' overdisp.est = overdispersed(example_data$ard,
#' known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
#' known_ind = c(1, 2, 4),
#' G1_ind = 1,
#' G2_ind = 2,
#' B2_ind = 4,
#' N = example_data$N,
#' warmup = 250,
#' iter = 500)
#'
#' # Compare size estimates
#' data.frame(true = example_data$subpop_sizes,
#' basic = colMeans(overdisp.est$betas))
#'
#' # Compare degree estimates
#' plot(example_data$degrees, colMeans(overdisp.est$alphas))
#'
#' # Look at overdispersion parameter
#' colMeans(overdisp.est$omegas)
overdispersed <-
  function(ard,
           known_sizes = NULL,
           known_ind = NULL,
           G1_ind = NULL,
           G2_ind = NULL,
           B2_ind = NULL,
           N = NULL,
           warmup = 1000,
           iter = 1500,
           refresh = NULL,
           thin = 1,
           verbose = FALSE,
           alpha_tune = 0.4,
           beta_tune = 0.2,
           omega_tune = 0.2,
           init = "MLE") {
    ## Extract dimensions
    N_i = nrow(ard) ## Number of respondents
    N_k = ncol(ard) ## Number of subpopulations

    if (is.null(refresh)) {
      ## By default, refresh every 10% of iterations
      refresh = round(iter / 10)
    }

    known_prevalences = known_sizes / N
    if (!is.null(G1_ind)) {
      Pg1 = sum(known_prevalences[G1_ind])
    }
    if (!is.null(G2_ind)) {
      Pg2 = sum(known_prevalences[G2_ind])
    }
    if (is.null(B2_ind)) {
      Pb2 = sum(known_prevalences[B2_ind])
    }


    ## Initialize parameters
    alphas = matrix(NA, nrow = iter, ncol = N_i)
    betas = matrix(NA, nrow = iter, ncol = N_k)
    omegas = matrix(NA, nrow = iter, ncol = N_k)
    mu_alpha = mu_beta = sigma_sq_alpha = sigma_sq_beta = rep(NA, iter)
    C1 = C2 = C = NA

    if (class(init) == "list") {
      if ("alpha" %in% names(init)) {
        alphas[1,] = init$alpha
      } else{
        alphas[1,] = stats::rnorm(N_i)
      }

      if ("beta" %in% names(init)) {
        betas[1,] = init$beta
      } else{
        betas[1,] = stats::rnorm(N_k)
      }

      if ("omega" %in% names(init)) {
        omegas[1,] = init$omega
      } else{
        omegas[1,] = 20
      }

    } else if (init == "random") {
      alphas[1,] = stats::rnorm(N_i)
      betas[1,] = stats::rnorm(N_k)
      omegas[1,] = 20
    } else{
      ## Based on MLE
      killworth_init = networkscaleup::killworth(ard, known_sizes, known_ind, N, model = "MLE")
      alphas[1,] = log(killworth_init$degrees)
      alphas[1, which(is.infinite(alphas[1,]))] = -10
      beta_vec = rep(NA, N_k)
      beta_vec[known_ind] = known_sizes
      beta_vec[-known_ind] = killworth_init$sizes
      beta_vec = log(beta_vec / N)
      betas[1,] = beta_vec
      omegas[1,] = 20
    }

    mu_alpha[1] = mean(alphas[1,])

    sigma_alpha_hat = mean((alphas[1,] - mu_alpha[1]) ^ 2)
    sigma_sq_alpha[1] = LaplacesDemon::rinvchisq(1, N_i - 1, sigma_alpha_hat)

    mu_beta[1] = mean(betas[1,])

    sigma_beta_hat = mean((betas[1,] - mu_beta[1]) ^ 2)
    sigma_sq_beta[1] = LaplacesDemon::rinvchisq(1, N_k - 1, sigma_beta_hat)



    for (ind in 2:iter) {
      ## Step 1
      for (i in 1:N_i) {
        alpha_prop = alphas[ind - 1, i] + stats::rnorm(1, 0, alpha_tune)
        zeta_prop = exp(alpha_prop + betas[ind - 1,]) / (omegas[ind - 1,] - 1)
        zeta_old = exp(alphas[ind - 1, i] + betas[ind - 1,]) / (omegas[ind - 1,] - 1)
        sum1 = sum(lgamma(ard[i,] + zeta_prop) - lgamma(zeta_prop) - zeta_prop * log(omegas[ind - 1,])) +
          stats::dnorm(alpha_prop, mu_alpha[ind - 1], sqrt(sigma_sq_alpha[ind - 1]), log = T)
        sum2 = sum(lgamma(ard[i,] + zeta_old) - lgamma(zeta_old) - zeta_old * log(omegas[ind - 1,])) +
          stats::dnorm(alphas[ind - 1, i], mu_alpha[ind - 1], sqrt(sigma_sq_alpha[ind - 1]), log = T)
        prob.acc = exp(sum1 - sum2)

        if (prob.acc > stats::runif(1)) {
          alphas[ind, i] = alpha_prop
        } else{
          alphas[ind, i] = alphas[ind - 1, i]
        }
      }

      ## Step 2
      for (k in 1:N_k) {
        beta_prop = betas[ind - 1, k] + stats::rnorm(1, 0, beta_tune)
        zeta_prop = exp(alphas[ind, ] + beta_prop) / (omegas[ind - 1, k] - 1)
        zeta_old = exp(alphas[ind, ] + betas[ind - 1, k]) / (omegas[ind - 1, k] - 1)
        sum1 = sum(lgamma(ard[, k] + zeta_prop) - lgamma(zeta_prop) - zeta_prop * log(omegas[ind - 1, k])) +
          stats::dnorm(beta_prop, mu_beta[ind - 1], sqrt(sigma_sq_beta[ind - 1]), log = T)
        sum2 = sum(lgamma(ard[, k] + zeta_old) - lgamma(zeta_old) - zeta_old * log(omegas[ind - 1, k])) +
          stats::dnorm(betas[ind - 1, k], mu_beta[ind - 1], sqrt(sigma_sq_beta[ind - 1]), log = T)
        prob.acc = exp(sum1 - sum2)

        if (prob.acc > stats::runif(1)) {
          betas[ind, k] = beta_prop
        } else{
          betas[ind, k] = betas[ind - 1, k]
        }
      }

      ## Step 3
      mu_alpha_hat = mean(alphas[ind,])
      mu_alpha[ind] = stats::rnorm(1, mu_alpha_hat, sqrt(sigma_sq_alpha[ind - 1] / 2))

      ## Step 4
      sigma_alpha_hat = mean((alphas[ind,] - mu_alpha[ind]) ^ 2)
      sigma_sq_alpha[ind] = LaplacesDemon::rinvchisq(1, N_i - 1, sigma_alpha_hat)

      ## Step 5
      mu_beta_hat = mean(betas[ind,])
      mu_beta[ind] = stats::rnorm(1, mu_beta_hat, sqrt(sigma_sq_beta[ind - 1] / 2))

      ## Step 6
      sigma_beta_hat = mean((betas[ind,] - mu_beta[ind]) ^ 2)
      sigma_sq_beta[ind] = LaplacesDemon::rinvchisq(1, N_k - 1, sigma_beta_hat)


      ## Step 7
      for (k in 1:N_k) {
        omega_prop = omegas[ind - 1, k] + stats::rnorm(1, 0, omega_tune)
        if (omega_prop > 1) {
          zeta_prop = exp(alphas[ind, ] + betas[ind, k]) / (omega_prop - 1)
          zeta_old = exp(alphas[ind, ] + betas[ind, k]) / (omegas[ind - 1, k] - 1)
          sum1 = sum(
            lgamma(ard[, k] + zeta_prop) - lgamma(zeta_prop) - zeta_prop * log(omega_prop) +
              ard[, k] * log((omega_prop - 1) / omega_prop)
          )
          sum2 = sum(
            lgamma(ard[, k] + zeta_old) - lgamma(zeta_old) - zeta_old * log(omegas[ind - 1, k]) +
              ard[, k] * log((omegas[ind - 1, k] - 1) / omegas[ind - 1, k])
          )
          prob.acc = exp(sum1 - sum2)

          if (prob.acc > stats::runif(1)) {
            omegas[ind, k] = omega_prop
          } else{
            omegas[ind, k] = omegas[ind - 1, k]
          }
        } else{
          omegas[ind, k] = omegas[ind - 1, k]
        }
      }

      ## Step 8
      if (is.null(G1_ind)) {
        ## Perform no scaling
        C = 0
      } else if (is.null(G2_ind) |
                 is.null(B2_ind)) {
        ## Perform scaling with only main
        C1 = log(sum(exp(betas[ind, G1_ind]) / Pg1))
        C = C1
      } else{
        ## Perform scaling with secondary groups
        C1 = log(sum(exp(betas[ind, G1_ind]) / Pg1))
        C2 = log(sum(exp(betas[ind, B2_ind]) / Pb2)) - log(sum(exp(betas[ind, G2_ind]) / Pg2))
        C = C1 + 1 / 2 * C2
      }

      alphas[ind,] = alphas[ind,] + C
      mu_alpha[ind] = mu_alpha[ind] + C
      betas[ind,] = betas[ind,] - C
      mu_beta[ind] = mu_beta[ind] - C


      if (verbose) {
        if (ind %% refresh == 0) {
          ## Match update format of Stan
          cat("Iteration: ",
              ind,
              " / ",
              iter,
              " [",
              round(ind / iter * 100),
              "%]",
              sep = "")
        }
      }
    }

    ## Burn-in and thin
    alphas = alphas[-c(1:warmup),]
    betas = betas[-c(1:warmup),]
    omegas = omegas[-c(1:warmup),]
    mu_alpha = mu_alpha[-c(1:warmup)]
    mu_beta = mu_beta[-c(1:warmup)]
    sigma_sq_alpha = sigma_sq_alpha[-c(1:warmup)]
    sigma_sq_beta = sigma_sq_beta[-c(1:warmup)]

    thin.ind = seq(1, nrow(alphas), by = thin)
    alphas = alphas[thin.ind,]
    betas = betas[thin.ind,]
    omegas = omegas[thin.ind,]
    mu_alpha = mu_alpha[thin.ind]
    mu_beta = mu_beta[thin.ind]
    sigma_sq_alpha = sigma_sq_alpha[thin.ind]
    sigma_sq_beta = sigma_sq_beta[thin.ind]


    return_list = list(
      alpha = alphas,
      degrees = exp(alphas),
      beta = betas,
      sizes = exp(betas) * N,
      omega = omegas,
      mu_alpha = mu_alpha,
      mu_beta = mu_beta,
      sigma_sq_alpha = sigma_sq_alpha,
      sigma_sq_beta = sigma_sq_beta
    )

    return(return_list)

  }

#' Generate simulated ARD
#'
#' @param n_i number of respondents (rows)
#' @param n_k number of groups (columns)
#' @param N total population size
#' @param p number of collected covariates
#' @param p_global_nonzero number of non-zero global covariates
#' @param p_local_nonzero number of non-zero local covariates
#' @param group_corr group correlation
#' @param degree_corr degree correlation
#' @param family sampling distribution
#' @param omega_range minimum and maximum omega for negative binomial overdispersion
#' @param alpha_mean mean of alphas
#' @param alpha_sd variance of alphas
#' @param eta correlation hyperparameter for LKJ prior
#' @param seed random seed
#'
#' @returns simulated ARD along with all true parameters
#' @export
#'
#' @examples make_ard(N = 10000, family = "poisson")
make_ard <- function(n_i = 500,
                     n_k = 20,
                     N = 1000000,
                     p = 0,
                     p_global_nonzero = 0,
                     p_local_nonzero = 0,
                     group_corr = FALSE,
                     degree_corr = FALSE,
                     family = c("poisson", "nbinomial"),
                     omega_range = c(1, 5),
                     alpha_mean = 5,
                     alpha_sd = 0.15,
                     eta = 3,
                     seed = NULL) {
  family <- match.arg(family)
  if (!is.null(seed)) {
    # Only set seed if provided
    set.seed(seed)
  }
  if (p_global_nonzero + p_local_nonzero > p) {
    stop("p_global_nonzero + p_local_nonzero cannot be greater than p")
  }
  if (group_corr & degree_corr) {
    stop("Cannot have group_corr and degree_corr together")
  }
  nk_prev <- stats::runif(n_k, 0.01, 0.15)
  nk_size <- round(nk_prev * N)
  omega <- stats::runif(n_k, omega_range[1], omega_range[2]) ## Unused if Poisson
  betas <- log(nk_size / N)
  alphas <- stats::rnorm(n_i, alpha_mean, alpha_sd)

  ## Correlated parameters
  tau.N <- stats::runif(n_k, min = 0.5, max = 1.5)


  if (group_corr) {
    mu <- log(1 / sqrt(1 + tau.N^2))
    tau <- sqrt(log(1 + tau.N^2))
    Omega <- trialr::rlkjcorr(1, n_k, eta = eta)
    L.omega <- t(chol(Omega))
    eps <- matrix(stats::rnorm(n_i * n_k), nrow = n_i, ncol = n_k)
    bias <- matrix(NA, nrow = n_i, ncol = n_k + 1)
    for (i in 1:n_i) {
      bias[i, -1] <- mu + diag(tau) %*% L.omega %*% eps[i, ]
    }
    bias[, 1] <- alphas
  } else if (degree_corr) {
    mu <- c(alpha_mean, log(1 / sqrt(1 + tau.N^2)))
    tau <- c(alpha_sd, sqrt(log(1 + tau.N^2)))
    Omega <- trialr::rlkjcorr(1, n_k + 1, eta = eta)
    L.omega <- t(chol(Omega))
    eps <- matrix(stats::rnorm(n_i * (n_k + 1)), nrow = n_i, ncol = n_k + 1)
    bias <- matrix(NA, nrow = n_i, ncol = n_k + 1)
    for (i in 1:n_i) {
      bias[i, ] <- mu + diag(tau) %*% L.omega %*% eps[i, ]
    }
    alphas <- bias[, 1]
  } else {
    bias <- matrix(0, nrow = n_i, ncol = n_k + 1)
    bias[, 1] <- alphas
  }

  ## Handle covariates
  if (p_global_nonzero > 0) {
    p_global_nonzero_ind <- sample(p, p_global_nonzero)
    p_local_nonzero_ind <- sample(c(1:p)[-p_global_nonzero_ind], p_local_nonzero)
    x_cov <- matrix(stats::runif(n_i * p, -1, 1), nrow = n_i, ncol = p)
    x_beta_global <- matrix(0, nrow = p, ncol = 1)
    x_beta_global[p_global_nonzero_ind, ] <- stats::runif(p_global_nonzero, -2, 2)
    x_beta_local <- matrix(0, nrow = p, ncol = n_k)
    x_beta_local[p_local_nonzero_ind, ] <- stats::runif(p_local_nonzero * n_k, -2, 2)
  } else {
    p_local_nonzero_ind <- sample(c(1:p), p_local_nonzero)
    x_cov <- matrix(stats::runif(n_i * p, -1, 1), nrow = n_i, ncol = p)
    x_beta_global <- matrix(0, nrow = p, ncol = 1)
    x_beta_local <- matrix(0, nrow = p, ncol = n_k)
    x_beta_local[p_local_nonzero_ind, ] <- stats::runif(p_local_nonzero * n_k, -2, 2)
  }

  # Center to have mean 0
  row_means <- rowMeans(x_beta_local)
  x_beta_local <- x_beta_local - row_means

  ## Simulate ARD
  ard <- matrix(NA, nrow = n_i, ncol = n_k)
  if (family == "poisson") {
    for (k in 1:n_k) {
      ard[, k] <- stats::rpois(
        n_i,
        lambda = exp(
          bias[, 1] + betas[k] +
            x_cov %*% x_beta_global +
            x_cov %*% x_beta_local[, k] +
            bias[, k + 1]
        )
      )
    }
  } else {
    for (k in 1:n_k) {
      ard[, k] <- stats::rnbinom(
        n_i,
        size = exp(
          bias[, 1] + betas[k] +
            x_cov %*% x_beta_global +
            x_cov %*% x_beta_local[, k] +
            bias[, k + 1]
        ) / (omega[k] - 1),
        prob = 1 / omega[k]
      )
    }
  }

  # The return gives all the information need to recreate this data in and out of function
  return(
    list(
      ard = ard,
      prev = nk_prev,
      size = nk_size,
      omega = omega,
      # TODO: unused if Poisson, should throw inside if statement
      eta = eta,
      # TODO: unused if uncorrelated, should throw inside if statement
      alphas = alphas,
      betas = betas,
      x_cov = x_cov,
      x_beta_global = x_beta_global,
      x_beta_local = x_beta_local,
      n_i = n_i,
      n_k = n_k,
      bias_mat = bias,
      seed = seed
    )
  )
}

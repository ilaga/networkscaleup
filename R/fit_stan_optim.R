library(cmdstanr)



#' Find MAP estimates of basic Poisson and Negative Binomial models using Stan optimization
#'
#' @param ard n_i by n_k ARD matrix
#' @param x_cov_global n_i by p_global covariate matrix of global covariates 
#' @param x_cov_local n_i by p_local covariate matrix of lobal covariates
#' @param family Distribution to fit, either "poisson" or "nbinomial"
#'
#' @return Stan fit
#' @export
#'
#' @examples
fit_stan_optim <- function(ard,
                           x_cov_global = NULL,
                           x_cov_local = NULL,
                           family = c("poisson", "nbinomial"),
                           ...) {
  ## Grab family
  family <- match.arg(family, c("poisson", "nbinomial"))

  n_local <- ncol(x_cov_local)
  n_global <- ncol(x_cov_global)

  if (is.null(x_cov_global) & is.null(x_cov_local)) {
    # No cov
    stan_data <- list(
      y = ard,
      n_i = nrow(ard),
      n_k = ncol(ard)
    )
    if (family == "poisson") {
      mod <- cmdstan_model("./Stan_Files/Poisson.stan")
    } else {
      mod <- cmdstan_model("./Stan_Files/Overdispersed.stan")
    }
  } else if (is.null(x_cov_global)) {
    # Only subpop cov
    stan_data <- list(
      y = ard,
      n_i = nrow(ard),
      n_k = ncol(ard),
      z_subpop_size = n_local,
      z_subpop = x_cov_local
    )
    if (family == "poisson") {
      mod <- cmdstan_model("./Stan_Files/Poisson_zsubpop.stan")
    } else {
      mod <- cmdstan_model("./Stan_Files/Overdispersed_zsubpop.stan")
    }
  } else if (is.null(x_cov_local)) {
    # Only global cov
    stan_data <- list(
      y = ard,
      n_i = nrow(ard),
      n_k = ncol(ard),
      z_global_size = n_global,
      z_global = x_cov_global
    )
    if (family == "poisson") {
      mod <- cmdstan_model("./Stan_Files/Poisson_zglobal.stan")
    } else {
      mod <- cmdstan_model("./Stan_Files/Overdispersed_zglobal.stan")
    }
  } else {
    # Both types of covariates
    stan_data <- list(
      y = ard,
      n_i = nrow(ard),
      n_k = ncol(ard),
      z_subpop_size = n_local,
      z_subpop = x_cov_local,
      z_global_size = n_global,
      z_global = x_cov_global
    )
    if (family == "poisson") {
      mod <- cmdstan_model("./Stan_Files/Poisson_zglobal_zsubpop.stan")
    } else {
      mod <- cmdstan_model("./Stan_Files/Overdispersed_zglobal_zsubpop.stan")
    }
  }





  fit <- mod$optimize(data = stan_data, ...)


  ## Add residuals
  if (family == "poisson") {
    pois_lambda_est <- fit$summary(variables = "mu")$estimate
    alphas <- fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(data.frame(value = c(ard)), family = "poisson", fit = pois_lambda_est)
    resids <- matrix(resid_vec, nrow = nrow(ard), ncol = ncol(ard))
  } else if (family == "nbinomial") {
    nb_prob_est <- fit$summary(variables = "inv_omegas")$estimate
    nb_size_est <- fit$summary(variables = "par1")$estimate
    alphas <- fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(
      data.frame(value = c(ard)),
      family = "nbinomial",
      size = nb_size_est,
      prob = rep(nb_prob_est, each = n_i)
    )
    resids <- matrix(resid_vec, nrow = nrow(ard), ncol = ncol(ard))
  }

  return_obj <- list(
    fit = fit,
    family = family,
    residuals = resids
  )

  return(return_obj)
}

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
fit_stan_optim <- function(ard,
                           x_cov_global = NULL,
                           x_cov_local = NULL,
                           family = c("poisson", "nbinomial"),
                           ...) {
  ## Grab family
  family <- match.arg(family, c("poisson", "nbinomial"))

  n_local <- ncol(x_cov_local)
  n_global <- ncol(x_cov_global)
  
  n_i = nrow(ard)
  n_k = ncol(ard)

  if (is.null(x_cov_global) & is.null(x_cov_local)) {
    # No cov
    stan_data <- list(
      y = ard,
      n_i = n_i,
      n_k = n_k
    )
    if (family == "poisson") {
      mod <- cmdstanr::cmdstan_model("./Stan_Files/Poisson.stan")
    } else {
      mod <- cmdstanr::cmdstan_model("./Stan_Files/Overdispersed.stan")
    }
  } else if (is.null(x_cov_global)) {
    # Only subpop cov
    stan_data <- list(
      y = ard,
      n_i = n_i,
      n_k = n_k,
      z_subpop_size = n_local,
      z_subpop = x_cov_local
    )
    if (family == "poisson") {
      mod <- cmdstanr::cmdstan_model("./Stan_Files/Poisson_zsubpop.stan")
    } else {
      mod <- cmdstanr::cmdstan_model("./Stan_Files/Overdispersed_zsubpop.stan")
    }
  } else if (is.null(x_cov_local)) {
    # Only global cov
    stan_data <- list(
      y = ard,
      n_i = n_i,
      n_k = n_k,
      z_global_size = n_global,
      z_global = x_cov_global
    )
    if (family == "poisson") {
      mod <- cmdstanr::cmdstan_model("./Stan_Files/Poisson_zglobal.stan")
    } else {
      mod <- cmdstanr::cmdstan_model("./Stan_Files/Overdispersed_zglobal.stan")
    }
  } else {
    # Both types of covariates
    stan_data <- list(
      y = ard,
      n_i = n_i,
      n_k = n_k,
      z_subpop_size = n_local,
      z_subpop = x_cov_local,
      z_global_size = n_global,
      z_global = x_cov_global
    )
    if (family == "poisson") {
      mod <- cmdstanr::cmdstan_model("./Stan_Files/Poisson_zglobal_zsubpop.stan")
    } else {
      mod <- cmdstanr::cmdstan_model("./Stan_Files/Overdispersed_zglobal_zsubpop.stan")
    }
  }
  fit <- mod$optimize(data = stan_data, ...)
  ## Add residuals
  if (family == "poisson") {
    pois_lambda_est <- fit$summary(variables = "mu")$estimate
    alphas <- fit$summary(variables = "alphas")$estimate
    # Pearson residuals
    pearson_vec <- construct_pearson(
      y = ard,
      family = "poisson",
      fit = pois_lambda_est
    )
    pearson_resids <- matrix(pearson_vec, nrow = n_i, ncol = n_k)
    # Randomized quantile residuals
    rqr_vec <- construct_rqr(
      y = ard,
      family = "poisson",
      model_fit = pois_lambda_est
    )
    rqr_resids <- matrix(rqr_vec, nrow = n_i, ncol = n_k)
    ## Return both sets of residuals
    return_obj <- list(
      fit = fit,
      family = family,
      n_i = n_i,
      n_k = n_k,
      alphas = alphas,
      pearson_residuals = pearson_resids,
      rqr = rqr_resids,
      x_cov_local = x_cov_local,
      x_cov_global = x_cov_global,
      mu = pois_lambda_est
    )
    
  } else if (family == "nbinomial") {
    nb_prob_est <- fit$summary(variables = "inv_omegas")$estimate
    nb_size_est <- fit$summary(variables = "par1")$estimate
    alphas <- fit$summary(variables = "alphas")$estimate
    # Pearson residuals
    pearson_vec <- construct_pearson(
      data.frame(value = c(ard)),
      family = "nbinomial",
      model_fit = fit)
    )
    pearson_resids <- matrix(pearson_vec, nrow = n_i, ncol = n_k)
    # Randomized quantile residuals
    rqr_vec <- construct_rqr(
      data.frame(value = c(ard)),
      family = "nbinomial",
      model_fit = fit)
    )
    rqr_resids <- matrix(rqr_vec, nrow = n_i, ncol = n_k)
    ## Return both sets of residuals
    return_obj <- list(
      fit = fit,
      family = family,
      n_i = n_i,
      n_k = n_k,
      alphas = alphas,
      pearson_residuals = pearson_resids,
      rqr = rqr_resids,
      x_cov_local = x_cov_local,
      x_cov_global = x_cov_global,
      size = nb_size_est,
      prob = nb_prob_est
    )
  }
  return_obj
}

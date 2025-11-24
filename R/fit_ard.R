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
                           family = c("poisson", "nbinomial")) {
  ## Grab family
  family <- match.arg(family, c("poisson", "nbinomial"))
  n_local <- ncol(x_cov_local)
  n_global <- ncol(x_cov_global)
  n_i <- nrow(ard)
  n_k <- ncol(ard)

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
  fit <- mod$optimize(data = stan_data)
  ## Add residuals
  if (family == "poisson") {
    fit_list <- list(fit = fit, 
                     mu = fit$summary(variables = "mu")$estimate,
                     family = "poisson")
    # pois_lambda_est <- fit$summary(variables = "mu")$estimate
    # fit$mu <- fit$summary(variables = "mu")$estimate
    alphas <- fit$summary(variables = "alphas")$estimate
    betas <- fit$summary(variables = "betas")$estimate
    # Pearson residuals
    pearson_vec <- construct_pearson(
      y = ard,
      model_fit = fit_list
    )
    pearson_resids <- matrix(pearson_vec, nrow = n_i, ncol = n_k)
    # Randomized quantile residuals
    rqr_vec <- construct_rqr(
      y = ard,
      model_fit = fit_list
    )
    rqr_resids <- matrix(rqr_vec, nrow = n_i, ncol = n_k)
    ## Return both sets of residuals
    return_obj <- list(
      fit = fit,
      family = family,
      n_i = n_i,
      n_k = n_k,
      alphas = alphas,
      betas = betas,
      pearson_residuals = pearson_resids,
      rqr = rqr_resids,
      x_cov_local = x_cov_local,
      x_cov_global = x_cov_global,
      mu = fit_list$mu
    )
    
  } else if (family == "nbinomial") {
    # nb_prob_est <- fit$summary(variables = "inv_omegas")$estimate
    # nb_size_est <- fit$summary(variables = "par1")$estimate
    par2_est <- fit$summary(variables = "par2")$estimate
    nb_prob_est <- par2_est / (par2_est + 1)
    fit_list <- list(fit = fit, 
                     size = fit$summary(variables = "par1")$estimate,
                     prob = nb_prob_est, family = "nbinomial")
    alphas <- fit$summary(variables = "alphas")$estimate
    betas <- fit$summary(variables = "betas")$estimate
    # Pearson residuals
    pearson_vec <- construct_pearson(
      y = ard,
      model_fit = fit_list)
    pearson_resids <- matrix(pearson_vec, nrow = n_i, ncol = n_k)
    # Randomized quantile residuals
    rqr_vec <- construct_rqr(
      y = ard,
      model_fit = fit_list)
    rqr_resids <- matrix(rqr_vec, nrow = n_i, ncol = n_k)
    ## Return both sets of residuals
    return_obj <- list(
      fit = fit,
      family = family,
      n_i = n_i,
      n_k = n_k,
      alphas = alphas,
      betas = betas,
      pearson_residuals = pearson_resids,
      rqr = rqr_resids,
      x_cov_local = x_cov_local,
      x_cov_global = x_cov_global,
      size = fit_list$size,
      prob = fit_list$prob
    )
  }
  return_obj
}





#' Fit basic Poisson and Negative Binomial models using glmmTMB
#'
#' @param ard n_i by n_k ARD matrix
#' @param x_cov_global n_i by p_global covariate matrix of global covariates
#' @param x_cov_local n_i by p_local covariate matrix of local covariates
#' @param family Distribution to fit, either "poisson" or "nbinomial"
#'
#' @return List containing fitted model and extracted parameters
#' @export
#' @import glmmTMB
fit_mle <- function(ard,
                    x_cov_global = NULL,
                    x_cov_local = NULL,
                    family = c("poisson", "nbinomial")) {
  ## Grab family
  family <- match.arg(family, c("poisson", "nbinomial"))
  n_i <- nrow(ard)
  n_k <- ncol(ard)

  # Reshape to long format
  y_long <- data.frame(
    i = rep(1:n_i, n_k),
    k = rep(1:n_k, each = n_i),
    count = as.vector(ard)
  )
  
  y_long$k = factor(y_long$k)
  y_long$i = factor(y_long$i)

  # Add covariates if provided
  if (!is.null(x_cov_global)) {
    # Global covariates: same value for all k's for each i
    for (j in 1:ncol(x_cov_global)) {
      y_long[[paste0("x_global_", j)]] <- x_cov_global[y_long$i, j]
    }
  }

  if (!is.null(x_cov_local)) {
    # Local covariates: vary by i (same for all k's? or by both i and k?)
    # Assuming they vary by i only
    for (j in 1:ncol(x_cov_local)) {
      y_long[[paste0("x_local_", j)]] <- x_cov_local[y_long$i, j]
    }
  }

  # Build formula
  fixed_terms <- c()
  random_terms <- "(1 | i) + (1 | k)"
  
  if (!is.null(x_cov_global)) {
    global_vars <- paste0("x_global_", 1:ncol(x_cov_global))
    fixed_terms <- c(fixed_terms, global_vars)
  }
  
  if (!is.null(x_cov_local)) {
    local_vars <- paste0("x_local_", 1:ncol(x_cov_local))
    fixed_terms <- c(fixed_terms, local_vars)
    
    # Add random slopes for local covariates with respect to k
    random_slopes <- paste0("(", local_vars, " | k)", collapse = " + ")
    random_terms <- paste(random_terms, "+", random_slopes)
  }
  
  # Construct full formula
  if (length(fixed_terms) == 0) {
    formula_str <- paste("count ~ 1 +", random_terms)
  } else {
    formula_str <- paste("count ~", paste(fixed_terms, collapse = " + "), "+", random_terms)
  }
  formula_full <- stats::as.formula(formula_str)
  # Fit model
  if (family == "poisson") {
    fit <- glmmTMB::glmmTMB(
      formula = formula_full,
      data = y_long,
      family = stats::poisson()
    )
  } else {
    fit <- glmmTMB::glmmTMB(
      formula = formula_full,
      data = y_long,
      family = glmmTMB::nbinom1,
      dispformula = ~ factor(k) - 1 # Group-specific dispersion
    )
  }

  # Extract parameters
  alphas <- glmmTMB::ranef(fit)$cond$i[[1]]
  betas <- glmmTMB::ranef(fit)$cond$k[[1]]

  # Get predictions
  y_long$mu_pred <- stats::predict(fit, type = "response")

  # Create model_fit object for residual functions
  if (family == "poisson") {
    mu_mat <- matrix(y_long$mu_pred, nrow = n_i, ncol = n_k)
    fit_list <- list(fit = fit, mu = mu_mat, family = "poisson")
    # Pearson residuals
    pearson_vec <- construct_pearson(
      y = ard,
      model_fit = fit_list
    )
    pearson_resids <- matrix(pearson_vec, nrow = n_i, ncol = n_k)

    # Randomized quantile residuals
    rqr_vec <- construct_rqr(
      y = ard,
      model_fit = fit_list
    )
    rqr_resids <- matrix(rqr_vec, nrow = n_i, ncol = n_k)

    return_obj <- list(
      fit = fit,
      family = family,
      n_i = n_i,
      n_k = n_k,
      alphas = alphas,
      betas = betas,
      pearson_residuals = pearson_resids,
      rqr = rqr_resids,
      x_cov_local = x_cov_local,
      x_cov_global = x_cov_global,
      mu = mu_mat
    )
  } else if (family == "nbinomial") {
    # Extract dispersion parameters
    phi_est <- exp(glmmTMB::fixef(fit)$disp) # phi from V = mu(1 + phi)
    omega_est <- phi_est + 1 # Convert to omega

    # Add group index
    y_long$k_idx <- as.numeric(factor(y_long$k))

    # Convert to rnbinom (size, prob) parameterization
    y_long$omega <- omega_est[y_long$k_idx]
    y_long$size <- y_long$mu_pred / (y_long$omega - 1)
    y_long$prob <- 1 / y_long$omega

    # Create matrices
    size_mat <- matrix(y_long$size, nrow = n_i, ncol = n_k)
    prob_mat <- matrix(y_long$prob, nrow = n_i, ncol = n_k)

    fit_list <- list(
      fit = fit,
      size = y_long$size,
      prob = prob_mat[1,],
      family = "nbinomial"
    )

    # Pearson residuals
    pearson_vec <- construct_pearson(
      y = ard,
      model_fit = fit_list
    )
    pearson_resids <- matrix(pearson_vec, nrow = n_i, ncol = n_k)

    # Randomized quantile residuals
    rqr_vec <- construct_rqr(
      y = ard,
      model_fit = fit_list
    )
    rqr_resids <- matrix(rqr_vec, nrow = n_i, ncol = n_k)
    return_obj <- list(
      fit = fit,
      family = family,
      n_i = n_i,
      n_k = n_k,
      alphas = alphas,
      betas = betas,
      pearson_residuals = pearson_resids,
      rqr = rqr_resids,
      x_cov_local = x_cov_local,
      x_cov_global = x_cov_global,
      size = y_long$size,
      prob = prob_mat[1,],
      omega = omega_est,
      phi = phi_est
    )
  }

  return(return_obj)
}

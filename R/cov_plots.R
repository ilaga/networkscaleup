library(ggplot2)

#' Title
#'
#' @param ard
#' @param x_cov
#'
#' @return
#' @export
#'
#' @examples
cov_plots <- function(model_fit,
                      ard,
                      x_cov,
                      cov_names = NULL,
                      family = c("poisson", "nbinomial"),
                      se = F) {
  family <- match.arg(family, c("poisson", "nbinomial"))
  
  ## Obtain residuals
  if (family == "poisson") {
    pois_lambda_est <- model_fit$summary(variables = "mu")$estimate
    alpha_est <- model_fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(data.frame(value = c(ard)), family = "poisson", fit = pois_lambda_est)
  } else if (family == "nbinomial") {
    nb_prob_est <- model_fit$summary(variables = "inv_omegas")$estimate
    nb_size_est <- model_fit$summary(variables = "par1")$estimate
    alpha_est <- model_fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(
      data.frame(value = c(ard)),
      family = "negbin",
      size = nb_size_est,
      prob = rep(nb_prob_est, each = n_i)
    )
  } else {
    stop("Invalid family argument. Must be one of poisson or nbinomial.")
  }
  
  resid_mat = matrix(resid_vec, nrow = nrow(ard), ncol = ncol(ard))
  
  ## Convert x_cov to data.frame, if not already
  if (!inherits(ard, "data.frame")) {
    ard <- data.frame(ard)
  }
  
  ## Convert x_cov to data.frame, if not already
  if (!inherits(x_cov, "data.frame")) {
    x_cov <- data.frame(x_cov)
  }
  
  if (is.null(cov_names)) {
    cov_names <- paste0("cov.", names(x_cov))
  }
  
  ard_x <- data.frame(resid = resid_mat, cov = x_cov, alpha = alpha_est)
  ard_long <- ard_x %>%
    pivot_longer(cols = starts_with("resid."),
                 names_to = "Group",
                 values_to = "resid")
  # ard_long$resid <- resid_vec
  
  
  ard_longer <- ard_long %>%
    pivot_longer(cols = starts_with("cov."),
                 names_to = "CovariateName",
                 values_to = "CovValue")
  
  ## Produce plot 1, averaged over groups
  gg1 <- ggplot(ard_longer, aes(
    x = CovValue,
    y = resid,
    col = Group,
    group = Group
  )) +
    geom_smooth(method = "lm", se = se) +
    # geom_smooth(se = se) +
    facet_wrap( ~ CovariateName, scales = "free") +
    theme_minimal() +
    labs(x = "Covariate Value", y = "Residual", title = "Residuals vs Covariates") +
    scale_color_discrete(guide = "none")
  
  gg2 <- ggplot(ard_longer, aes(x = CovValue, y = alpha, col = CovariateName)) +
    # geom_point(alpha = 0.6) +
    geom_smooth(method = "lm") +
    # facet_wrap(~ CovariateName, scales = "free_x") +
    theme_minimal() +
    labs(x = "Covariate Value", y = "Estimated Respondent Effect", title = "Respondent Effect vs Covariate")
  
  return(list(Group_plot = gg1, Respondent_plot = gg2))
}

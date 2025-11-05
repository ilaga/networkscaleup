#' Compute Surrogate Residuals for ARD Models
#'
#' @param y the ARD matrix
#' @param model_fit if poisson distribution, rate matrix
#' @param family the distribution to fit, which will choose prev pars to be
#' specified
#'
#' @returns a vector of residuals (column by column)
#' @export
get_surrogate <- function(y, model_fit = NULL,
                          family = c("binomial", "negbin", "poisson")) {
  
  n_samp <- nrow(y)
  
  family <- match.arg(family, c("poisson", "nbinomial", "binomial"))
  if (family == "poisson") {
    pois_lambda_est <- model_fit$mu
    mu_vec <- as.numeric(pois_lambda_est)
  } else if (family == "nbinomial") {
    nb_prob_est <- model_fit$prob
    nb_size_est <- model_fit$size
    size_vec <- as.numeric(nb_size_est)
    prob_vec <- as.numeric(nb_prob_est)
    prob_vec <- rep(prob_vec, each = n_samp)
  } else {
    stop("Invalid family argument. Must be one of poisson or nbinomial.",
         call. = FALSE)
  }
  ## TO DO: Add Binomial correctly here
  
  resid_len <- if (!is.null(mu_vec)) {
    length(mu_vec)
  } else if (!is.null(size_vec)) {
    length(size_vec)
  }
  resid <- rep(NA, resid_len)
  if (family == "binomial") {
    for (i in 1:length(resid)) {
      y_sim <- stats::rbinom(1, size = size_vec[i], 
                             prob = p)
      F_val <- stats::pbinom(y_sim, size = size_vec[i],
                             prob = p)
      # Inverse standard normal transformation
      resid[i] <- stats::qnorm(F_val)
    }
  } else if (family == "nbinomial") {
    for (i in 1:length(resid)) {
      y_sim <- stats::rnbinom(n = 1, size = size_vec[i],
                              prob = prob_vec[i])
      F_val <- stats::pnbinom(y_sim, size = size_vec[i], 
                              prob = prob_vec[i])
      resid[i] <- stats::qnorm(F_val)
    }
  } else if (family == "poisson") {
    for (i in 1:length(resid)) {
      y_sim <- stats::rpois(n = 1, lambda = mu_vec[i])
      F_val <- stats::ppois(y_sim, lambda = mu_vec[i])
      resid[i] <- stats::qnorm(F_val)
    }
  } else {
    stop("Invalid distribution")
  }
  resid
}

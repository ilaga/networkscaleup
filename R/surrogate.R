#' Compute Surrogate Residuals for ARD Models
#'
#' @param p if binomial distribution, success probability (single value)
#' @param fit if poisson distribution, rate matrix
#' @param size if negative binomial, size matrix
#' @param prob if negative binomial, size matrix
#' @param dist the distribution to fit, which will choose prev pars to be
#' specified
#'
#' @returns a vector of residuals (column by column)
#' @export
get_surrogate <- function(p = NULL, fit = NULL, size = NULL, prob = NULL,
                          dist = c("binomial", "negbin", "poisson")) {
  
  ## TO DO - update checks 
  if (dist == "poisson" & !is.matrix(fit)) {
    stop("Supplied Poisson fit must be a matrix",
         call. = FALSE)
  }
  if (dist == "negbin" & (!is.matrix(size) | !is.matrix(prob)) ) {
    stop("Supplied Negative Binomial fit must be specified with matrices",
         call. = FALSE)
  }
  
  ## transform matrices to vector
  mu_vec <- as.numeric(fit)
  size_vec <- as.numeric(size)
  prob_vec <- as.numeric(prob)
  
  resid_len <- if (!is.null(fit)) {
    length(fit)
  } else if (!is.null(size)) {
    length(size)
  }
  
  resid <- rep(NA, resid_len)
  
  if (dist == "binomial") {
    for (i in 1:length(resid)) {
      y_sim <- stats::rbinom(1, size = size_vec[i], 
                             prob = p)
      F_val <- stats::pbinom(y_sim, size = size_vec[i],
                             prob = p)
      # Inverse standard normal transformation
      resid[i] <- stats::qnorm(F_val)
    }
  } else if (dist == "negbin") {
    for (i in 1:length(resid)) {
      y_sim <- stats::rnbinom(n = 1, size = size_vec[i],
                              prob = prob_vec[i])
      F_val <- stats::pnbinom(y_sim, size = size_vec[i], 
                              prob = prob_vec[i])
      resid[i] <- stats::qnorm(F_val)
    }
  } else if (dist == "poisson") {
    for (i in 1:length(resid)) {
      y_sim <- stats::rpois(n = 1, lambda = mu_vec[i])
      F_val <- stats::ppois(y_sim, lambda = mu_vec[i])
      resid[i] <- stats::qnorm(F_val)
    }
  } else {
    stop("Invalid distribution")
  }
  
  return(resid)
}

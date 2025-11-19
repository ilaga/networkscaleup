#' Compute Randomized Quantile Residuals for ARD Models
#'
#' @param y ard matrix
#' @param model_fit list of details of fitted model
#' specified
#'
#' @returns a vector of residuals (column by column)
#' @export
get_rqr <- function(y, model_fit = NULL) {
  
  n_i <- nrow(y)
  
  family <- model_fit$family
  family <- match.arg(family, c("poisson", "nbinomial"))
  if (family == "poisson") {
    pois_lambda_est <- model_fit$mu
    fit_vec <- as.numeric(pois_lambda_est)
  } else if (family == "nbinomial") {
    nb_prob_est <- model_fit$prob
    nb_size_est <- model_fit$size
    size_vec <- as.numeric(nb_size_est)
    prob_vec <- as.numeric(nb_prob_est)
    prob_vec <- rep(prob_vec, each = n_i)
  } else {
    stop("Invalid family argument. Must be one of poisson or nbinomial.",
         call. = FALSE)
  }
  ## transform matrices to vector
  y_vec <- as.numeric(y)
  rqr <- rep(NA, length(y_vec))
  
  if (family == "binomial") {
    stop("Not implemented yet, can't specify p")
    for (i in 1:length(y_vec)) {
      # Get CDF at y[i] and at y[i] - 1
      F_lower <- NA
      if (y_vec[i] == 0) {
        F_lower <- 0
      } else {
        F_lower <- stats::pbinom(y_vec[i] - 1, size = size_vec[i], prob = 0.5)
      }
      F_upper <- stats::pbinom(y_vec[i], size = size_vec[i], prob = 0.5)
      
      # Sample a uniform value between F_lower and F_upper
      u <- stats::runif(1, min = F_lower, max = F_upper)
      
      # Inverse standard normal transformation
      rqr[i] <- stats::qnorm(u)
    }
  } else if (family == "nbinomial") {
    for (i in 1:length(y_vec)) {
      # Get CDF at y[i] and at y[i] - 1
      F_lower <- NA
      if (y_vec[i] == 0) {
        F_lower <- 0
      } else {
        F_lower <- stats::pnbinom(y_vec[i] - 1, size = size_vec[i],
                                  prob = prob_vec[i])
      }
      F_upper <- stats::pnbinom(y_vec[i], size = size_vec[i],
                                prob = prob_vec[i])
      
      # Sample a uniform value between F_lower and F_upper
      u <- stats::runif(1, min = F_lower, max = F_upper)
      
      # Inverse standard normal transformation
      rqr[i] <- stats::qnorm(u)
    }
  } else if (family == "poisson") {
    for (i in 1:length(y_vec)) {
      # Get CDF at y[i] and at y[i] - 1
      # F_lower <- NA
      # if (y_vec[i] == 0) {
      #   F_lower <- 0
      # } else {
      #   F_lower <- stats::ppois(y_vec[i] - 1, lambda = mu_vec[i])
      # }
      # F_upper <- stats::ppois(y_vec[i], lambda = mu_vec[i])
      # 
      # # Sample a uniform value between F_lower and F_upper
      # u <- stats::runif(1, min = F_lower, max = F_upper)
      # 
      # # Inverse standard normal transformation
      # rqr[i] <- stats::qnorm(u)
      rqr[i] <- rqr_pois_logs(y_vec[i], fit_vec[i])
    }
  } else {
    stop("Invalid distribution")
  }
  rqr
}




#' log computed uniform quantile
#' 
#' This function is not vectorized
#'
#' @param logFl log of lower value 
#' @param logFu log of upper value
#'
#' @returns log value of uniform between Flower and Fupper
log_mix_uniform <- function(logFl, logFu) {
  u <- stats::runif(1)
  if(is.infinite(logFl)) {
    logu <- logFu + log(u)
  } else{
    a <- logFu - logFl
    logu <- logFl + base::log1p( u * (base::exp(a) - 1) )
  }
  logu
}


#' compute numerically stable Poisson rqr
#'
#' @param y observed value
#' @param mu mean value of poisson
#' @param eps precision parameter
#'
#' @returns appropriate randomized quantile residual
rqr_pois_logs <- function(y, mu, eps = 1e-12) {
  logFu <- stats::ppois(y,     mu, log.p = TRUE)
  logFl <- stats::ppois(y - 1, mu, log.p = TRUE)
  logu  <- log_mix_uniform(logFl, logFu)
  # Clip in probability space *after* exponentiating
  u     <- base::exp(logu)
  u     <- base::pmin(base::pmax(u, eps), 1 - eps)
  stats::qnorm(u)
}
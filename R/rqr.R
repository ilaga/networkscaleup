#' Compute Randomized Quantile Residuals for ARD Models
#'
#' @param y ard matrix
#' @param p if binomial distribution, success probability (single value)
#' @param fit if poisson distribution, rate matrix
#' @param size if negative binomial, size matrix
#' @param prob if negative binomial, size matrix
#' @param dist the distribution to fit, which will choose prev pars to be
#' specified
#'
#' @returns a vector of residuals (column by column)
#' @export
get_rqr <- function(y, p = NULL, fit = NULL, size = NULL, prob = NULL,
                    dist = c("binomial", "negbin", "poisson")) {
  
  
  if (!is.matrix(y)) {
    stop("ARD must be a matrix", call. = FALSE)
  }
  if (dist == "poisson" & !is.matrix(fit)) {
    stop("Supplied Poisson fit must be a matrix",
         call. = FALSE)
  }
  if(dist == "poisson" & (!identical(dim(y), dim(fit))) ){
    stop("Parameters don't match ARD matrix", call. = FALSE)
  }
  if (dist == "negbin" & (!is.matrix(size) | !is.matrix(prob)) ) {
    stop("Supplied Negative Binomial fit must be specified with matrices",
         call. = FALSE)
  }
  if(dist == "negbin" & (!identical(dim(y), dim(size))) ){
    stop("Parameters don't match ARD matrix", call. = FALSE)
  }
  
  
  
  
  ## transform matrices to vector
  y_vec <- as.numeric(y)
  mu_vec <- as.numeric(fit)
  size_vec <- as.numeric(size)
  prob_vec <- as.numeric(prob)
  
  rqr <- rep(NA, length(y_vec))
  
  if (dist == "binomial") {
    for (i in 1:length(y_vec)) {
      # Get CDF at y[i] and at y[i] - 1
      F_lower <- NA
      if (y_vec[i] == 0) {
        F_lower <- 0
      } else {
        F_lower <- stats::pbinom(y_vec[i] - 1, size = size_vec[i], prob = p)
      }
      F_upper <- stats::pbinom(y_vec[i], size = size_vec[i], prob = p)
      
      # Sample a uniform value between F_lower and F_upper
      u <- stats::runif(1, min = F_lower, max = F_upper)
      
      # Inverse standard normal transformation
      rqr[i] <- stats::qnorm(u)
    }
  } else if (dist == "negbin") {
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
  } else if (dist == "poisson") {
    for (i in 1:length(y_vec)) {
      # Get CDF at y[i] and at y[i] - 1
      F_lower <- NA
      if (y_vec[i] == 0) {
        F_lower <- 0
      } else {
        F_lower <- stats::ppois(y_vec[i] - 1, lambda = mu_vec[i])
      }
      F_upper <- stats::ppois(y_vec[i], lambda = mu_vec[i])
      
      # Sample a uniform value between F_lower and F_upper
      u <- stats::runif(1, min = F_lower, max = F_upper)
      
      # Inverse standard normal transformation
      rqr[i] <- stats::qnorm(u)
    }
  } else {
    stop("Invalid distribution")
  }
  
  return(rqr)
}

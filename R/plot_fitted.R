#' Plot residuals against fitted values
#'
#' @param y ARD matrix (may be needed)
#' @param p if binomial distribution, success probability (single value)
#' @param fit if poisson distribution, rate matrix
#' @param size if negative binomial, size matrix
#' @param prob if negative binomial, prob matrix
#' @param dist the distribution to be fit
#' @param resid the type of residuals to be used
#'
#' @returns a ggplot showing fitted values against residuals
#' @export
#'
plot_fitted <- function(y, p = NULL, fit = NULL, size = NULL, prob = NULL,
                        dist = c("binomial", "negbin", "poisson"),
                        resid = c("rqr", "pearson", "surrogate")) {
  
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
  
  if(resid == "rqr") {
    resids <- get_rqr(y, p = p, fit = fit, size = size, 
                      prob = prob, dist = dist)
    plot_label <- "Randomized Quantile Residuals"
  } else if(resid == "pearson") {
    resids <- construct_pearson(y, dist = dist, fit = fit,
                                size = size, prob = prob)
    plot_label <- "Pearson Residuals"
  } else if(resid == "surrogate") {
    resids <- get_surrogate(p = p, fit = fit, size = size, 
                            prob = prob, dist = dist)
    plot_label <- "Surrogate Residuals"
  } else{
    stop("Invalid residuals specified")
  }
  # then construct data for plot
  if(dist == "poisson") {
    plot_data <- base::data.frame(fit = as.numeric(fit), resid = resids)
  } else if(dist == "negbin") {
    plot_data <- base::data.frame(fit = as.numeric(size) * 
                                    as.numeric(1 - prob)/(as.numeric(prob)),
                                  resid = resids)
  } else if(dist == "binomial") {
    plot_data <- base::data.frame(fit = as.numeric(size) * p,
                                  resid = resids)
  }
  ## then construct the plot
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$fit,
                                 y = .data$resid)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Fitted Values",
                  y = plot_label) +
    ggplot2::theme_bw()
}
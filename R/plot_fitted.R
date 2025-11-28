#' Plot residuals against fitted values
#'
#' @param y ARD matrix (may be needed)
#' @param model_fit fitted model
#' @param resid the type of residuals to be used
#'
#' @returns a ggplot showing fitted values against residuals
#' @export
#'
plot_fitted <- function(y, model_fit = NULL,
                        resid = c("rqr", "pearson", "surrogate")) {
  
  n_samp <- nrow(y)
  family <- model_fit$family
  family <- match.arg(family, c("poisson", "nbinomial", "binomial"))
  if (family == "poisson") {
    pois_lambda_est <- model_fit$mu
    fit_vec <- as.numeric(pois_lambda_est)
  } else if (family == "nbinomial") {
    nb_prob_est <- model_fit$prob
    nb_size_est <- model_fit$size
    size_vec <- as.numeric(nb_size_est)
    prob_vec <- as.numeric(nb_prob_est)
    full_prob <- rep(prob_vec, each = n_samp)
  } else {
    stop("Invalid family argument. Must be one of poisson or nbinomial.",
         call. = FALSE)
  }
  
  if(resid == "rqr") {
    resids <- construct_rqr(y, model_fit = model_fit)
    plot_label <- "Randomized Quantile Residuals"
  } else if(resid == "pearson") {
    resids <- construct_pearson(y, model_fit = model_fit)
    plot_label <- "Pearson Residuals"
  } else if(resid == "surrogate") {
    resids <- get_surrogate(y, model_fit = model_fit)
    plot_label <- "Surrogate Residuals"
  } else{
    stop("Invalid residuals specified")
  }
  # then construct data for plot
  if(family == "poisson") {
    plot_data <- base::data.frame(fit = as.numeric(fit_vec), resid = resids)
  } else if(family == "negbin") {
    plot_data <- base::data.frame(fit = as.numeric(size_vec) * 
                                    as.numeric(1 - full_prob)/(as.numeric(full_prob)),
                                  resid = resids)
  } else if(family == "binomial") {
    ## TO DO, not completed
    stop("Not completed, do not use")
    size <- 100
    p <- 0.5
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
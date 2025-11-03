#' Hanging Rootogram for Fitted ARD Model
#'
#' @param y ard matrix 
#' @param fit matrix of (estimated) means of each entry if poisson ARD model
#' @param family whether ard model poisson or negative binomial
#' @param size if negative binomial, matrix of size estimate
#' @param prob if negative binomial, matrix of prob estimate
#' @param width width of bars
#' @param x_max the maximum x value to display
#'
#' @return a ggplot of the hanging rootogram
#' @export
#' @importFrom rlang .data
hang_rootogram_ard <- function(y,
                               model_fit = NULL,         # fitted stan model
                               family = "poisson",   # default familyribution to fit
                               size = NULL,        # NB dispersion if needed
                               prob = NULL,        # NB size if needed 
                               width  = 0.9,       # bar width (0–1)
                               x_max = NULL){      
  
  family <- match.arg(family, c("poisson", "nbinomial"))
  if (family == "poisson") {
    pois_lambda_est <- model_fit$summary(variables = "mu")$estimate
    fit_vec <- as.numeric(pois_lambda_est)
  } else if (family == "nbinomial") {
    nb_prob_est <- model_fit$summary(variables = "inv_omegas")$estimate
    nb_size_est <- model_fit$summary(variables = "par1")$estimate
    size_vec <- as.numeric(nb_size_est)
    prob_vec <- as.numeric(nb_prob_est)
  } else {
    stop("Invalid family argument. Must be one of poisson or nbinomial.",
         call. = FALSE)
  }
  
  ## transform matrix to vector
  y_vec <- as.numeric(y)
  
  if(is.null(x_max)){
    ## 1. support (integer counts)
    k <- 0:max(y_vec, floor(max(y_vec) * 1.25))   
  }
  else{
    # k <- 0:max(y_vec, floor(x_max))
    k <- 0:floor(x_max)
  }
  # generous upper limit
  
  ## 2. observed bin counts
  obs_counts <- as.numeric(table(factor(y_vec, levels = k)))
  
  ## 3. expected bin counts
  exp_counts <- vapply(
    k,
    function(j) {
      if (family == "poisson")
        sum(stats::dpois(j, lambda = fit_vec))
      else if (family == "nbinomial") {
        if (is.null(size_vec) | is.null(prob_vec))
          stop("Please supply 'size' and 'prob' for the negative-binomial model")
        sum(stats::dnbinom(j, size = size_vec, prob = prob_vec))
      }
    },
    numeric(1)
  )
  
  ## 4. square-root transform and bar coordinates
  obs_root <- sqrt(obs_counts)
  exp_root <- sqrt(exp_counts)
  baseline <- exp_root                 # hanging baseline, changed from -
  tips     <- exp_root - obs_root       # expected - observed
  
  df <- data.frame(
    k,
    xmin = k - width/2,
    xmax = k + width/2,
    ymin = pmin(baseline, tips),      # bottom of bar
    ymax = pmax(baseline, tips)       # top of bar
    # pos  = tips >= baseline         # TRUE = observed > expected
  ) |> 
    dplyr::mutate(middle = (.data$xmin + .data$xmax)/2)
  
  ## add poisson or negative binomial model
  if(family == "poisson"){
    plot_lab <- "Poisson"
  }
  else if(family == "nbinomial") {
    plot_lab <- "Negative Binomial"
  }
  
  ## 5. plot
  ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$xmin,
                   xmax = .data$xmax,
                   ymin = .data$ymin,
                   ymax = .data$ymax),
      colour = "lightgray",
      fill = "gray"
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(
      x = "Count",
      y = expression(sqrt(count)),
      title = "Hanging Rootogram",
      subtitle = plot_lab
    ) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes(x = .data$middle, y = .data$ymax),
                       col = "red") +
    ggplot2::geom_point(ggplot2::aes(x = .data$middle, y = .data$ymax),
                        col = "red") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) 
}

#' Hanging Rootogram for Fitted ARD Model
#'
#' @param y ard matrix 
#' @param fit matrix of (estimated) means of each entry if poisson ARD model
#' @param dist whether ard model poisson or negative binomial
#' @param size if negative binomial, matrix of size estimate
#' @param prob if negative binomial, matrix of prob estimate
#' @param width width of bars
#' @param x_max the maximum x value to display
#'
#' @return a ggplot of the hanging rootogram
#' @export
#' @importFrom rlang .data
hang_rootogram_ard <- function(y,
                               fit = NULL,         # poisson rates
                               dist = "poisson",   # default distribution to fit
                               size = NULL,        # NB dispersion if needed
                               prob = NULL,        # NB size if needed 
                               width  = 0.9,       # bar width (0–1)
                               x_max = NULL){      
  
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
  fit_vec <- as.numeric(fit)
  size_vec <- as.numeric(size)
  prob_vec <- as.numeric(prob)
  
  
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
      if (dist == "poisson")
        sum(stats::dpois(j, lambda = fit_vec))
      else if (dist == "negbin") {
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
  if(dist == "poisson"){
    plot_lab <- "Poisson"
  }
  else if(dist == "negbin") {
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

#' Hanging Rootogram for Fitted ARD Model
#'
#' @param y ard matrix 
#' @param fit matrix of (estimated) means of each entry if poisson ARD model
#' @param width width of bars
#' @param x_max the maximum x value to display
#'
#' @return a ggplot of the hanging rootogram
#' @export
#' @importFrom rlang .data
hang_rootogram_ard <- function(y,
                               model_fit,         # fitted stan model
                               width  = 0.9,       # bar width (0–1)
                               x_max = NULL){      
  
  n_i <- nrow(y)
  
  family <- model_fit$family
  if (family == "poisson") {
    pois_lambda_est <- model_fit$mu
    fit_vec <- as.numeric(pois_lambda_est)
  } else if (family == "nbinomial") {
    prob_vec <- rep(model_fit$prob, each = n_i)
    size_vec <- model_fit$size
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

#' Dispersion Metric for Fitted ARD Model
#'
#' @param y ard matrix 
#' @param fit matrix of (estimated) means of each entry if poisson ARD model
#'
#' @return a ggplot of the hanging rootogram
#' @export
#' @importFrom rlang .data
dispersion_metric <- function(y, model_fit) {
  
  n_i <- nrow(y)
  n_k <- ncol(y)
  
  family <- model_fit$family
  if (family == "poisson") {
    pois_lambda_est <- matrix(model_fit$mu, nrow = n_i, ncol = n_k)
  } else {
    stop("Invalid family argument. Metrics assume a Poisson likelihood.",
         call. = FALSE)
  }
  
  # Initialize results storage
  dispersion_stats <- data.frame(
    column = 1:n_k,
    statistic = numeric(n_k),
    df = numeric(n_k),
    p_value = numeric(n_k),
    dispersion_ratio = numeric(n_k)
  )
  
  # Calculate dispersion test for each column
  for (k in 1:n_k) {
    y_k <- y[, k]
    mu_k <- pois_lambda_est[, k]
    
    # Pearson residuals
    pearson_resid <- (y_k - mu_k) / sqrt(mu_k)
    
    # Dispersion statistic (sum of squared Pearson residuals)
    disp_stat <- sum(pearson_resid^2)
    
    # Degrees of freedom (n - number of parameters)
    # For simple case, use n - 1; adjust if you know exact df
    df <- n_i - 1
    
    # P-value from chi-squared distribution
    p_val <- pchisq(disp_stat, df = df, lower.tail = FALSE)
    
    # Dispersion ratio (should be ~1 for Poisson)
    disp_ratio <- disp_stat / df
    
    dispersion_stats[k, ] <- c(k, disp_stat, df, p_val, disp_ratio)
  }
  
  # Create visualization
  plot_data <- data.frame(
    column = factor(dispersion_stats$column),
    dispersion_ratio = dispersion_stats$dispersion_ratio,
    p_value = dispersion_stats$p_value,
    significant = dispersion_stats$p_value < 0.05
  )
  
  # Dispersion ratio plot
  disp_plot <- ggplot2::ggplot(plot_data, 
                               ggplot2::aes(x = column, y = dispersion_ratio,
                                            fill = significant)) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", 
                        color = "red", linewidth = 1) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "gray70"),
                               labels = c("FALSE" = "Not significant", 
                                          "TRUE" = "Significant (p < 0.05)")) +
    ggplot2::labs(
      x = "Column",
      y = "Dispersion Ratio (χ²/df)",
      title = "Dispersion Test by Column",
      subtitle = "Ratio = 1 indicates Poisson fit; >1 indicates overdispersion",
      fill = "Dispersion Test"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)
    )
  
  # P-value plot
  pval_plot <- ggplot2::ggplot(plot_data,
                               ggplot2::aes(x = column, y = -log10(p_value),
                                            fill = significant)) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                        color = "red", linewidth = 1) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "gray70"),
                               labels = c("FALSE" = "Not significant",
                                          "TRUE" = "Significant (p < 0.05)")) +
    ggplot2::labs(
      x = "Column",
      y = "-log10(p-value)",
      title = "Dispersion Test P-values",
      subtitle = "Dashed line at p = 0.05",
      fill = "Dispersion Test"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)
    )
  
  # Combined plot
  combined_plot <- gridExtra::grid.arrange(disp_plot, pval_plot, ncol = 2)
  
  # Summary statistics
  summary_stats <- list(
    n_significant = sum(dispersion_stats$p_value < 0.05),
    n_columns = n_k,
    prop_significant = mean(dispersion_stats$p_value < 0.05),
    mean_dispersion_ratio = mean(dispersion_stats$dispersion_ratio),
    median_dispersion_ratio = median(dispersion_stats$dispersion_ratio)
  )
  
  return(list(
    dispersion_stats = dispersion_stats,
    plots = list(
      dispersion_plot = disp_plot,
      pvalue_plot = pval_plot,
      combined_plot = combined_plot
    ),
    summary = summary_stats
  ))
}

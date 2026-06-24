#' Tracy-Widom test for residual group correlation
#'
#' @param model_fit fitted model object
#' @param correction correction constant, either "none", "half"
#' @param plot a logical, whether to return a ggplot density plot of TW with observed statistic
#'
#' @return a list containing test statistic, p-value, and diagnostic plots
#' @export
#' @importFrom rlang .data
tw_group_corr_test <- function(model_fit,
                               correction = c("none", "half"),
                               plot = TRUE) {
  correction <- match.arg(correction)
  
  # Residual matrix
  resid_mat <- model_fit$rqr
  n_i <- nrow(resid_mat)
  n_k <- ncol(resid_mat)
  
  # Eigenvalues of covariance matrix
  S <- (1 / (n_i - 1)) * t(resid_mat) %*% resid_mat
  eigenvalues <- eigen(S)$values
  lambda_max <- max(eigenvalues)
  
  # Centering and scaling constants
  
  if (correction == "none") {
    mu_n <- (sqrt(n_i) + sqrt(n_k))^2 / n_i
    sigma_n <- (sqrt(n_i) + sqrt(n_k)) / n_i *
      (1 / sqrt(n_i) + 1 / sqrt(n_k))^(1 / 3)
  } else if (correction == "half") {
    mu_n <- (sqrt(n_i - 1 / 2) + sqrt(n_k - 1 / 2))^2 / n_i
    sigma_n <- (sqrt(n_i - 1 / 2) + sqrt(n_k - 1 / 2)) / n_i *
      (1 / sqrt(n_i - 1 / 2) + 1 / sqrt(n_k - 1 / 2))^(1 / 3)
  }
  
  # Tracy-Widom statistic and p-value
  
  tw_stat <- (lambda_max - mu_n) / sigma_n
  p_value <- 1 - RMTstat::ptw(tw_stat, beta = 1)
  
  # Tracy-Widom density plot with observed statistic
  if (plot) {
    tw_density_plot <- plot_tw_statistic(tw_stat)
  } else {
    tw_density_plot <- NULL
  }
  
  
  return(list(
    lambda_max = lambda_max,
    tw_statistic = tw_stat,
    p_value = p_value,
    mu_n = mu_n,
    sigma_n = sigma_n,
    all_eigenvalues = eigenvalues,
    n_i = n_i,
    n_k = n_k,
    tw_density_plot = tw_density_plot
  ))
}



#' Plot the TW Test Statistic
#'
#' @param tw_stat Tracey-Widom test statistic
#' @param beta order param of TW distribution
#' @param base_x_limits limits for axis 
#' @param base_size plot sizes
#' @param title_size title size
#' @param axis_title_size axis title size
#' @param axis_text_size axis text size
#' @param label_size label size
#'
#' @returns a ggplot
plot_tw_statistic <- function(tw_stat,
                              beta = 1,
                              base_x_limits = c(-5, 10),
                              base_size = 9,
                              title_size = 10,
                              axis_title_size = 8,
                              axis_text_size = 7,
                              label_size = 2.8) {
  base_x_min <- base_x_limits[1]
  base_x_max <- base_x_limits[2]
  
  x_min <- min(base_x_min, tw_stat)
  x_max <- max(base_x_max, tw_stat)
  
  x_pad <- 0.08 * (x_max - x_min)
  x_limits <- c(x_min - x_pad, x_max + x_pad)
  
  density_grid <- data.frame(
    x = seq(x_limits[1], x_limits[2], length.out = 1000)
  )
  density_grid$y <- RMTstat::dtw(density_grid$x, beta = beta)
  
  y_max <- max(density_grid$y, na.rm = TRUE)
  
  rel_pos <- (tw_stat - x_limits[1]) / diff(x_limits)
  
  if (rel_pos < 0.5) {
    label_x <- tw_stat + 0.08 * diff(x_limits)
    label_hjust <- 0
  } else {
    label_x <- tw_stat - 0.08 * diff(x_limits)
    label_hjust <- 1
  }
  
  label_x <- min(
    max(label_x, x_limits[1] + 0.03 * diff(x_limits)),
    x_limits[2] - 0.03 * diff(x_limits)
  )
  
  label_y <- 0.90 * y_max
  
  ggplot2::ggplot(density_grid, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$y),
      linewidth = 1,
      color = "black"
    ) +
    ggplot2::geom_vline(
      xintercept = tw_stat,
      color = "red",
      linewidth = 1
    ) +
    ggplot2::geom_label(
      data = data.frame(
        x = label_x,
        y = label_y,
        label = sprintf("T = %.2f", tw_stat)
      ),
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        label = .data$label
      ),
      inherit.aes = FALSE,
      hjust = label_hjust,
      vjust = 0.5,
      color = "red",
      label.size = 0.25,
      fill = "white",
      size = label_size
    ) +
    ggplot2::labs(
      title = "Tracy-Widom Density with Observed Statistic",
      x = expression(TW[1]~value),
      y = "Density"
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = c(0, 1.08 * y_max),
      clip = "off"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      plot.margin = ggplot2::margin(5.5, 20, 5.5, 20)
    )
}
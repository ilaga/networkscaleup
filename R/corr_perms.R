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
    tw_density_plot <- ggplot2::ggplot(data.frame(x = c(-5, 10)), ggplot2::aes(.data$x)) +
      ggplot2::stat_function(
        fun = RMTstat::dtw,
        args = list(beta = 1),
        linewidth = 1, color = "black"
      ) +
      ggplot2::geom_vline(xintercept = tw_stat, color = "red", linewidth = 1) +
      ggplot2::annotate("text",
        x = tw_stat,
        y = 0,
        label = sprintf(" T = %.2f", tw_stat),
        vjust = -1, color = "red"
      ) +
      ggplot2::labs(
        title = "Tracy-Widom Density with Observed Statistic",
        x = "TW_1 value",
        y = "Density"
      ) +
      ggplot2::theme_minimal()
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

# PCA-based residual test for residual correlation (permute each column independently)
#' Title
#'
#' @param ard ARD data
#' @param model_fit list containing fitted model and additional details
#' @param b number of replications to do
#'
#' @returns list containing plots, test statistics, etc
#' @export
#'
pca_group_corr_test <- function(ard,
                                model_fit,
                                b = 1000) {
  ## Obtain residuals
  resid_mat <- model_fit$rqr
  # PCA on observed data
  obs_pca <- stats::prcomp(resid_mat, center = TRUE, scale. = TRUE)
  obs_var <-
    obs_pca$sdev[1]^2 / sum(obs_pca$sdev^2) # first PC variance
  obs_var_all <-
    obs_pca$sdev^2 / sum(obs_pca$sdev^2) # all PCs (scree)

  # Permutation test
  n_cols <- ncol(resid_mat)
  perm_var <- numeric(b)
  perm_var_all <- matrix(NA, nrow = b, ncol = n_cols)

  for (i in 1:b) {
    # Permute each column independently
    perm_df <- apply(resid_mat, 2, sample)

    tmp_pca <- stats::prcomp(perm_df, center = TRUE, scale. = TRUE)
    perm_var[i] <- tmp_pca$sdev[1]^2 / sum(tmp_pca$sdev^2)
    perm_var_all[i, ] <- tmp_pca$sdev^2 / sum(tmp_pca$sdev^2)
  }

  # Two-sided p-value
  p_val <-
    mean(abs(perm_var - mean(perm_var)) >= abs(obs_var - mean(perm_var)))

  hist_plot <- ggplot2::ggplot(data.frame(perm_var), ggplot2::aes(x = perm_var)) +
    ggplot2::geom_histogram(binwidth = diff(range(perm_var)) / 30, fill = "gray80", color = "black") +
    ggplot2::geom_vline(xintercept = obs_var, color = "red", linewidth = 1) +
    ggplot2::labs(
      x = "Variance explained by first PC",
      y = "Count",
      title = "PCA Residual Test"
    ) +
    ggplot2::theme_minimal()


  perm_df <- as.data.frame(t(perm_var_all)) |>
    dplyr::mutate(PC = 1:dplyr::n()) |>
    tidyr::pivot_longer(
      cols = -PC,
      names_to = "Permutation",
      values_to = "Variance"
    )

  obs_df <- data.frame(
    PC = 1:length(obs_var_all),
    Variance = obs_var_all
  )

  ymax <- max(c(obs_var_all, perm_var_all), na.rm = TRUE)

  scree_plot <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = perm_df,
      ggplot2::aes(x = PC, y = Variance, group = Permutation),
      color = rgb(0, 0, 0, 0.2)
    ) +
    ggplot2::geom_line(
      data = obs_df,
      ggplot2::aes(x = PC, y = Variance),
      color = "red",
      linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = obs_df,
      ggplot2::aes(x = PC, y = Variance),
      color = "red"
    ) +
    ggplot2::labs(
      x = "PC index",
      y = "Proportion variance explained",
      title = "Scree plots (obs vs permuted)"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::guides(color = "none")

  return(
    list(
      test_stat = obs_var,
      bootstrap_stat = perm_var,
      p_value = p_val,
      obs_scree = obs_var_all,
      perm_scree = perm_var_all,
      hist_plot = hist_plot,
      scree_plot = scree_plot
    )
  )
}




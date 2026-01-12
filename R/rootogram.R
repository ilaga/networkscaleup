#' Hanging Rootogram for Fitted ARD Model
#'
#' @param ard ard matrix
#' @param model_fit fitted model object
#' @param width width of bars
#' @param x_max the maximum x value to display
#' @param by_group logical; if TRUE, create separate rootograms for each column (group)
#'
#' @return a ggplot of the hanging rootogram (single plot if by_group=FALSE, combined plot if by_group=TRUE)
#' @export
#' @importFrom rlang .data
hang_rootogram_ard <- function(ard,
                               model_fit,
                               width = 0.9,
                               x_max = NULL,
                               by_group = FALSE) {
  n_i <- nrow(ard)
  n_k <- ncol(ard)
  family <- model_fit$family
  if (family == "poisson") {
    pois_lambda_est <- model_fit$mu
    pois_lambda_mat <- matrix(pois_lambda_est, nrow = n_i)
  } else if (family == "nbinomial") {
    prob_vec <- rep(model_fit$prob, each = n_i)
    size_vec <- model_fit$size
    prob_vec_mat <- matrix(prob_vec, nrow = n_i)
    size_vec_mat <- matrix(size_vec, nrow = n_i)
  } else {
    stop("Invalid family argument. Must be one of poisson or nbinomial.",
      call. = FALSE
    )
  }

  # Create column name vector for group rootogram
  col_names <- colnames(ard)
  if (is.null(col_names)) {
    col_names <- names(ard)
  }

  if (!is.null(col_names)) {
    group_labels <- col_names
  } else {
    group_labels <- paste0("Group ", 1:n_k)
  }

  # Helper function to create a single rootogram
  create_rootogram <- function(y_vec, fit_vec = NULL, group_label = NULL,
                               size_vec = NULL, prob_vec = NULL) {
    if (is.null(x_max)) {
      k <- 0:max(y_vec, floor(max(y_vec) * 1.25))
    } else {
      k <- 0:floor(x_max)
    }

    # Observed bin counts
    obs_counts <- as.numeric(table(factor(y_vec, levels = k)))

    # Expected bin counts
    exp_counts <- vapply(
      k,
      function(j) {
        if (family == "poisson") {
          sum(stats::dpois(j, lambda = fit_vec))
        } else if (family == "nbinomial") {
          sum(stats::dnbinom(j, size = size_vec, prob = prob_vec))
        }
      },
      numeric(1)
    )

    # Square-root transform and bar coordinates
    obs_root <- sqrt(obs_counts)
    exp_root <- sqrt(exp_counts)
    baseline <- exp_root
    tips <- exp_root - obs_root

    df <- data.frame(
      k,
      xmin = k - width / 2,
      xmax = k + width / 2,
      ymin = pmin(baseline, tips),
      ymax = pmax(baseline, tips)
    ) |>
      dplyr::mutate(middle = (.data$xmin + .data$xmax) / 2)

    # Plot label
    if (family == "poisson") {
      plot_lab <- "Poisson"
    } else if (family == "nbinomial") {
      plot_lab <- "Negative Binomial"
    }

    # Add group label if provided
    title <- if (!is.null(group_label)) {
      paste0("Hanging Rootogram - ", group_label)
    } else {
      "Hanging Rootogram"
    }

    # Create plot
    ggplot2::ggplot(df) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax
        ),
        colour = "lightgray",
        fill = "gray"
      ) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::labs(
        x = "Count",
        y = expression(sqrt(count)),
        title = title,
        subtitle = plot_lab
      ) +
      ggplot2::theme_bw() +
      ggplot2::geom_line(ggplot2::aes(x = .data$middle, y = .data$ymax),
        col = "red"
      ) +
      ggplot2::geom_point(ggplot2::aes(x = .data$middle, y = .data$ymax),
        col = "red"
      ) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 6))
  }

  # If by_group = FALSE, create single rootogram with all data
  if (!by_group) {
    y_vec <- as.numeric(ard)

    if (family == "nbinomial") {
      return(create_rootogram(y_vec, size_vec = size_vec, prob_vec = prob_vec))
    } else {
      fit_vec <- as.numeric(pois_lambda_est)
      return(create_rootogram(y_vec, fit_vec))
    }
  } else {
    # If by_group = TRUE, create separate rootograms for each column
    plot_list <- list()
    for (k in 1:n_k) {
      y_vec_k <- ard[, k]

      if (family == "nbinomial") {
        size_vec_k <- size_vec_mat[, k]
        prob_vec_k <- prob_vec_mat[, k]
        plot_list[[k]] <- create_rootogram(y_vec_k,
          group_label = group_labels[k],
          size_vec = size_vec_k,
          prob_vec = prob_vec_k
        )
      } else {
        fit_vec_k <- pois_lambda_mat[, k]
        plot_list[[k]] <- create_rootogram(y_vec_k, fit_vec_k,
          group_label = group_labels[k]
        )
      }
    }

    combined_plot <- gridExtra::arrangeGrob(grobs = plot_list, ncol = min(3, n_k))
    plot(combined_plot)

    return(combined_plot)
  }
}


#' Dispersion Metric for Fitted ARD Model
#'
#' @param ard ard matrix
#' @param model_fit list of fitted model and details
#'
#' @return a ggplot of the hanging rootogram
#' @export
#' @importFrom rlang .data
dispersion_metric <- function(ard, model_fit) {
  n_i <- nrow(ard)
  n_k <- ncol(ard)
  
  p <- NCOL(model_fit$x_cov_local) + NCOL(model_fit$x_cov_global)
  
  family <- model_fit$family
  if (family == "poisson") {
    pois_lambda_est <- matrix(model_fit$mu, nrow = n_i, ncol = n_k)
  } else {
    stop("Invalid family argument. Metrics assume a Poisson likelihood.",
      call. = FALSE
    )
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
    y_k <- ard[, k]
    mu_k <- pois_lambda_est[, k]

    # Pearson residuals
    pearson_resid <- (y_k - mu_k) / sqrt(mu_k)

    # Dispersion statistic (sum of squared Pearson residuals)
    disp_stat <- sum(pearson_resid^2)

    # Degrees of freedom (n - number of parameters)
    df <- n_i - p

    # P-value from chi-squared distribution
    p_val <- stats::pchisq(disp_stat, df = df, lower.tail = FALSE)

    # Dispersion ratio (should be approx 1 for Poisson)
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
  disp_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$column,
      y = .data$dispersion_ratio,
      fill = .data$significant
    )
  ) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_hline(
      yintercept = 1, linetype = "dashed",
      color = "red", linewidth = 1
    ) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "coral", "FALSE" = "gray70"),
      labels = c(
        "FALSE" = "Not significant",
        "TRUE" = "Significant (p < 0.05)"
      )
    ) +
    ggplot2::labs(
      x = "Column",
      y = expression("Dispersion Ratio " * (frac(chi^2, df))),
      title = "Dispersion Test by Column",
      fill = "Dispersion Test"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9)
    )

  # Summary statistics
  summary_stats <- list(
    n_significant = sum(dispersion_stats$p_value < 0.05),
    n_columns = n_k,
    prop_significant = mean(dispersion_stats$p_value < 0.05),
    mean_dispersion_ratio = mean(dispersion_stats$dispersion_ratio),
    median_dispersion_ratio = stats::median(dispersion_stats$dispersion_ratio)
  )

  return(list(
    dispersion_stats = dispersion_stats,
    plots = disp_plot,
    summary = summary_stats
  ))
}

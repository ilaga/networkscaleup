#' Covariance plots
#'
#' Plots of the estimated covariance structure from a given fitted model
#'
#' @param model_fit output of a fitted model, including pearson residuals
#' @param ard ard matrix
#' @param x_cov covariate matrix
#' @param method the method to use
#' @param se whether to compute standard errors of estimates
#'
#' @return a list of ggplots, corresponding to covariance structure
#' @export
#'
cov_plots <- function(model_fit,
                      ard,
                      x_cov,
                      method = "lm",
                      se = F) {
  ## Grab family
  family <- model_fit$family

  ## Obtain residuals
  resid_mat <- model_fit$pearson_residuals
  alpha_est <- model_fit$alphas

  ## Convert ard to data.frame, if not already
  if (!inherits(ard, "data.frame")) {
    ard <- data.frame(ard)
  }

  ## Convert x_cov to data.frame, if not already
  if (!inherits(x_cov, "data.frame")) {
    x_cov <- data.frame(x_cov) ## Auto-names to X1:Xp if no names
  }

  ard_x <- data.frame(resid = resid_mat, cov = x_cov, alpha = alpha_est)
  ard_long <- ard_x |> 
    tidyr::pivot_longer(
      cols = starts_with("resid."),
      names_to = "Group",
      values_to = "resid"
    )


  ard_longer <- ard_long |> 
    tidyr::pivot_longer(
      cols = starts_with("cov."),
      names_to = "cov_names",
      values_to = "CovValue"
    ) |> 
    mutate(cov_label = str_remove(cov_names, "^cov\\."))


  ## Produce plot 1, group-specific plots
  gg1 <- ggplot2::ggplot(ard_longer, ggplot2::aes(
    x = CovValue,
    y = resid,
    col = Group,
    group = Group
  )) +
    ggplot2::geom_smooth(method = method, se = se) +
    ggplot2::facet_wrap(~cov_label, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Covariate Value", y = "Residual", title = "Residuals vs Covariates") +
    ggplot2::scale_color_discrete(guide = "none")


  ## Find endpoints to add covariate label
  label_df <- ard_longer |>
    split(.$cov_label) |>
    purrr::imap_dfr(~ {
      sm <- suppressMessages(
        ggplot2::ggplot_build(
          ggplot2::ggplot(.x, ggplot2::aes(x = CovValue, y = alpha)) +
            ggplot2::geom_smooth(method = method, se = se)
        )$data[[1]]
      )
      
      sm |>
        dplyr::filter(x == max(x, na.rm = TRUE)) |>
        dplyr::mutate(cov_label = .y)
    })

  # Produce plot 2, averaged over groups
  gg2 <- ggplot2::ggplot(ard_longer, 
                         ggplot2::aes(x = CovValue, y = alpha, color = cov_label)) +
    ggplot2::geom_smooth(method = method, se = se) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = x, y = y, label = cov_label, color = cov_label),
      hjust = -0.1,
      show.legend = FALSE
    ) +
    ggplot2::expand_limits(x = max(label_df$x) * 1.1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Covariate Value",
      y = "Estimated Respondent Effect",
      title = "Respondent Effect vs Covariate",
      color = "Covariate"
    )

  return(list(Group_plot = gg1, Respondent_plot = gg2))
}

#' Construct tibble from ARD matrix
#'
#' @param y ARD matrix
#'
#' @return a tibble of ARD, with columns for row/col index
#'
make_ard_tidy <- function(y){
  ard_df <- as.data.frame(y)
  colnames(ard_df) <- 1:ncol(y)
  
  long_ard_tidy <- tibble::as_tibble(ard_df) |>                    # <- the matrix
    # as_tibble(.name_repair = "universal") |>       # keep column names as-is
    tibble::rowid_to_column("row") |>                    # add a numeric row index
    tidyr::pivot_longer(
      cols      = -row,                            # everything except the row id
      names_to  = "col",
      values_to = "value"
    ) |> 
    dplyr::mutate(col = as.integer(col)) |> 
    dplyr::arrange(col, row)
  long_ard_tidy
}




#' Compute Pearson Residuals for ARD matrix and fitted model
#'
#' @param y ARD matrix y
#' @param model_fit estimated model 
#' 
#' @return a vector (column by column) of corresponding residuals from ARD matrix
#' @export
#'
#' @importFrom rlang .data
construct_pearson <- function(y, model_fit) {
  long_ard <- make_ard_tidy(y)
  n_i <- nrow(y)
  family <- model_fit$family
  family <- match.arg(family, c("poisson", "nbinomial"))
  if (family == "poisson") {
    pois_lambda_est <- model_fit$mu
    fit_vec <- as.numeric(pois_lambda_est)
  } else if (family == "nbinomial") {
    nb_prob_est <- model_fit$prob
    nb_size_est <- model_fit$size
    size_vec <- as.numeric(nb_size_est)
    prob_vec <- as.numeric(nb_prob_est)
    prob_vec <- rep(prob_vec, each = n_i)
  } else {
    stop("Invalid family argument. Must be one of poisson or nbinomial.",
         call. = FALSE)
  }
  
  ## transform matrix to vector
  y_vec <- as.numeric(y)
  if(family == "poisson") {
    long_ard |> 
      dplyr::mutate(est = fit_vec,
                    resid = (.data$value - .data$est)/sqrt(.data$est)) |> 
      dplyr::pull(.data$resid)
  }
  else if(family == "nbinomial") {
    long_ard |> 
      dplyr::mutate(size = size_vec,
                    prob = prob_vec,
                    est = .data$size * (1 - .data$prob)/.data$prob,
                    resid = (.data$value -
                               .data$est)/sqrt(.data$est/.data$prob)) |> 
      dplyr::pull(.data$resid)
  }
  else{
    stop("Invalid distribution")
  }
}

#' Compute Randomized Quantile Residuals for ARD Models
#'
#' @param y ard matrix
#' @param model_fit fitted model, along with required details
#'
#' @returns a vector of residuals (column by column)
#' @export
construct_rqr <- function(y, model_fit) {
  n_i <- nrow(y)
  family <- model_fit$family
  family <- match.arg(family, c("poisson", "nbinomial", "binomial"))
  if (family == "poisson") {
    pois_lambda_est <- model_fit$mu
    mu_vec <- as.numeric(pois_lambda_est)
  } else if (family == "nbinomial") {
    nb_prob_est <- model_fit$prob
    nb_size_est <- model_fit$size
    size_vec <- as.numeric(nb_size_est)
    prob_vec <- as.numeric(nb_prob_est)
    prob_vec <- rep(prob_vec, each = n_i)
  } else {
    stop("Invalid family argument. Must be one of poisson or nbinomial.",
         call. = FALSE)
  }
  ## TO DO: Add Binomial correctly here and extract the p
  y_vec <- as.numeric(y)
  rqr <- rep(NA, length(y_vec))
  
  if (family == "binomial") {
    stop("Not implemented yet, the probability not passed in below.")
    for (i in 1:length(y_vec)) {
      # Get CDF at y[i] and at y[i] - 1
      F_lower <- NA
      if (y_vec[i] == 0) {
        F_lower <- 0
      } else {
        F_lower <- stats::pbinom(y_vec[i] - 1, size = size_vec[i], prob = 0.5)
      }
      F_upper <- stats::pbinom(y_vec[i], size = size_vec[i], prob = 0.5)
      
      # Sample a uniform value between F_lower and F_upper
      u <- stats::runif(1, min = F_lower, max = F_upper)
      
      # Inverse standard normal transformation
      rqr[i] <- stats::qnorm(u)
    }
  } else if (family == "nbinomial") {
    for (i in 1:length(y_vec)) {
      ## to avoid some numerical issues
      rqr[i] <- rqr_nbinom_logs(y_vec[i], size_vec[i], prob_vec[i])
    }
  } else if (family == "poisson") {
    for (i in 1:length(y_vec)) {
      ## to avoid some numerical issues
      rqr[i] <- rqr_pois_logs(y_vec[i], mu_vec[i])
    }
  } else {
    stop("Invalid family")
  }
  
  return(rqr)
}



#' Construct heatmap of residuals
#'
#' @param ard_residuals a vector (column wise) of estimated residuals
#' @param y an ard matrix
#'
#' @return A ggplot of residual heatmap
#' @export
#'
residual_heatmap <- function(ard_residuals, y){
  long_ard <- make_ard_tidy(y)
  long_ard$residuals <- ard_residuals
  n_cols <- max(long_ard$col)
  n_rows <- max(long_ard$row)
  ggplot2::ggplot(long_ard, ggplot2::aes(y = row, x = col, fill = .data$residuals)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient2(
      low  = "red",       # negative
      mid  = "white",     # zero
      high = "blue",      # positive
      midpoint = 0
    ) +
    ggplot2::labs(x = "Column", y = "Row", fill = "Residual") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed(ratio = n_cols / n_rows)
}



#' Construction Residual (row/column) correlation matrix
#'
#' @param ard_residuals vector of residuals
#' @param y ard matrix y
#' @param type type of correlation to use (row or column)
#'
#' @return a ggplot of the specified correlation matrix
#' @export
#'
#' @importFrom rlang .data
residual_correlation <- function(ard_residuals, y,
                                 type = "column") {
  
  long_ard <- make_ard_tidy(y)
  long_ard$residuals <- ard_residuals
  n_cols <- max(long_ard$col)
  n_rows <- max(long_ard$row)
  resid_mat <- long_ard |>
    dplyr::select(-.data$value) |> 
    tidyr::pivot_wider(names_from = col,
                       values_from = .data$residuals) |> 
    dplyr::select(-row) |>                     # drop row id
    as.matrix()
  
  if(type == "column"){
    cors <- stats::cor(resid_mat, use = "pairwise.complete.obs",
                method = "pearson")
    cors_long <- cors |>
      as.data.frame() |>
      tibble::rownames_to_column("row") |>
      tidyr::pivot_longer(-row, names_to = "col", values_to = "corr") |> 
      dplyr::mutate(col = factor(col, levels = 1:n_cols),
                    row = factor(row, levels = n_cols:1))
    plot_label <- "Column Wise Residual Correlation"
    plot_axis <- ggplot2::element_text(angle = 45, hjust = 1)
  }
  if(type == "row"){
    if(nrow(y) > 500){
      stop("ARD too large for row-wise correlation plot", call. = FALSE)
    }
    cors <- stats::cor(t(resid_mat), use = "pairwise.complete.obs",
                method = "pearson")
    cors_long <- cors |>
      as.data.frame() |>
      tibble::rownames_to_column("row") |>
      tidyr::pivot_longer(-row, names_to = "col", values_to = "corr") |> 
      dplyr::mutate(col = readr::parse_number(col)) |> 
      dplyr::mutate(col = factor(col, levels = 1:n_rows),
                    row = factor(row, levels = n_rows:1))
    plot_label <- "Row Wise Residual Correlation"
    plot_axis <- ggplot2::element_blank()
  }
  
  ggplot2::ggplot(cors_long, ggplot2::aes(col, row, fill = .data$corr)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),          # full correlation range
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      name = "r"
    ) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = plot_label) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = plot_axis,
      axis.text.y = plot_axis,
      legend.key.height = ggplot2::unit(3, "mm"),
      legend.key.width  = ggplot2::unit(4, "mm"),
      legend.position   = "right"
    )
}


#' log computed uniform quantile
#'
#' @param logFl log of lower value 
#' @param logFu log of upper value
#'
#' @returns log value of uniform between Flower and Fupper
log_mix_uniform <- function(logFl, logFu) {
  u <- stats::runif(1)
  if(is.infinite(logFl)) {
    logu <- logFu + log(u)
  } else{
    a <- logFu - logFl
    logu <- logFl + base::log1p( u * (base::exp(a) - 1) )
  }
  logu
}


#' compute numerically stable Poisson rqr
#'
#' @param y observed value
#' @param mu mean value of poisson
#' @param eps precision parameter
#'
#' @returns appropriate randomized quantile residual
rqr_pois_logs <- function(y, mu, eps = 1e-12) {
  logFu <- stats::ppois(y, mu, log.p = TRUE)
  logFl <- stats::ppois(y - 1, mu, log.p = TRUE)
  logu  <- log_mix_uniform(logFl, logFu)
  # Clip in probability space *after* exponentiating
  u     <- base::exp(logu)
  u     <- base::pmin(base::pmax(u, eps), 1 - eps)
  stats::qnorm(u)
}

#' compute numerically stable negative binomial rqr
#'
#' @param y observed value
#' @param size size parameter
#' @param prob prob parameter
#' @param eps precision parameter
#'
#' @returns appropriate randomized quantile residual
rqr_nbinom_logs <- function(y, size, prob, eps = 1e-12) {
  logFu <- stats::pnbinom(y, size = size, prob = prob, log.p = TRUE)
  logFl <- stats::pnbinom(y - 1, size = size, prob = prob, log.p = TRUE)
  logu  <- log_mix_uniform(logFl, logFu)
  # Clip in probability space *after* exponentiating
  u     <- base::exp(logu)
  u     <- base::pmin(base::pmax(u, eps), 1 - eps)
  stats::qnorm(u)
}

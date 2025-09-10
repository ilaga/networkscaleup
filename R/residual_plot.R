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
#' @param dist poisson or negative binomial model 
#' @param fit (posterior) matrix means for fitted parameters
#' @param size matrix for negative binomial. vector of length n * k (column wise)
#' @param prob matrix for negative binomial. vector of length n * k (column wise)
#'
#' @return a vector (column by column) of corresponding residuals from ARD matrix
#' @export
#'
#' @importFrom rlang .data
construct_pearson <- function(y, dist = "poisson",
                              fit = NULL, size = NULL, prob = NULL) {
  long_ard <- make_ard_tidy(y)
  
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
  
  fit_vec <- as.numeric(fit)
  size_vec <- as.numeric(size)
  prob_vec <- as.numeric(prob)
  
  if(dist == "poisson") {
    long_ard |> 
      dplyr::mutate(est = fit_vec,
                    resid = (.data$value - .data$est)/sqrt(.data$est)) |> 
      dplyr::pull(.data$resid)
  }
  else if(dist == "negbin") {
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


#' Construct heatmap of ARD residuals
#'
#' @param ard_residuals vector (column by column) of fitted residuals
#' @param y ARD matrix y
#'
#' @return a heatmap of the residuals
#' @export
#'
residual_heatmap <- function(ard_residuals, y){
  
  if (!is.matrix(y)) {
    stop("ARD must be a matrix", call. = FALSE)
  }
  
  long_ard <- make_ard_tidy(y)
  long_ard$residuals <- ard_residuals
  n_cols <- max(long_ard$col)
  n_rows <- max(long_ard$row)
  resid_plot <- ggplot2::ggplot(long_ard,
                                ggplot2::aes(y = row, x = col,
                                             fill = .data$residuals)) +
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
    ggplot2::coord_fixed(ratio = n_cols / n_rows) +
    NULL
  # TO DO - get rid of this message when plotting
  resid_plot
}



#' Construction Residual (row/column) correlation matrix
#'
#' @param ard_residuals a vector (column by column) of ard residuals
#' @param y ard matrix y
#' @param type the type of correlation to be computed (row or column)
#'
#' @return a ggplot of the correlation matrix
#' @export
#'
#' @importFrom rlang .data
residual_correlation <- function(ard_residuals, y,
                                 type = "column") {
  
  if (!is.matrix(y)) {
    stop("ARD must be a matrix", call. = FALSE)
  }
  if (length(y) != length(ard_residuals)) {
    stop("ARD matrix must match size of the residual vector", call. = FALSE)
  }
  
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
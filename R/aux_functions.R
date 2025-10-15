library(ggplot2)
library(tidyverse)


#' Title
#'
#' @param y ARD matrix
#'
#' @return
#' @export
#'
#' @examples
make_ard_tidy <- function(y){
  ard_df <- as.data.frame(y)
  colnames(ard_df) <- 1:ncol(y)
  
  long_ard_tidy <- as_tibble(ard_df) |>                    # <- the matrix
    # as_tibble(.name_repair = "universal") |>       # keep column names as-is
    rowid_to_column("row") |>                    # add a numeric row index
    pivot_longer(
      cols      = -row,                            # everything except the row id
      names_to  = "col",
      values_to = "value"
    ) |> 
    mutate(col = as.integer(col)) |> 
    arrange(col, row)
  long_ard_tidy
}



#' Title
#'
#' @param y ARD matrix y
#' @param family 
#' @param fit (posterior) means for fitted parameters
#' @param size for negative binomial. vector of length n * k (column wise)
#' @param prob for negative binomial. vector of length n * k (column wise)
#'
#' @return
#' @export
#'
#' @examples
construct_pearson <- function(y, family = "poisson",
                              fit = NULL, size = NULL, prob = NULL) {
  long_ard <- make_ard_tidy(y)
  n_i <- nrow(y)
  n_k <- ncol(y)
  if(length(prob) != n_i * n_k & family == "nbinomial") {
    stop("You have not specified the probability vector for the negative binomial
         the correct way. Please check the documentation.")
  }
  if(family == "poisson") {
    long_ard |> 
      mutate(est = fit,
             resid = (value - est)/sqrt(est)) |> 
      pull(resid)
  }
  else if(family == "nbinomial") {
    long_ard |> 
      mutate(est = size * (1 -prob)/prob,
             resid = (value - est)/sqrt(est/prob)) |> 
      pull(resid)
  }
}


#' Title
#'
#' @param ard_residuals 
#' @param long_ard 
#'
#' @return
#' @export
#'
#' @examples
residual_heatmap <- function(ard_residuals, y){
  long_ard <- make_ard_tidy(y)
  long_ard$residuals <- ard_residuals
  n_cols <- max(long_ard$col)
  n_rows <- max(long_ard$row)
  ggplot(long_ard, aes(y = row, x = col, fill = residuals)) +
    geom_tile() +
    coord_fixed() +
    scale_fill_gradient2(
      low  = "red",       # negative
      mid  = "white",     # zero
      high = "blue",      # positive
      midpoint = 0
    ) +
    labs(x = "Column", y = "Row", fill = "Residual") +
    theme_minimal() +
    theme(
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    coord_fixed(ratio = n_cols / n_rows)
}



#' Construction Residual (row/column) correlation matrix
#'
#' @param ard_residuals 
#' @param y ard matrix y
#' @param type 
#'
#' @return
#' @export
#'
#' @examples
residual_correlation <- function(ard_residuals, y,
                                 type = "column") {
  long_ard <- make_ard_tidy(y)
  long_ard$residuals <- ard_residuals
  n_cols <- max(long_ard$col)
  n_rows <- max(long_ard$row)
  resid_mat <- long_ard |>
    select(-value) |> 
    pivot_wider(names_from = col,
                values_from = residuals) |> 
    select(-row) |>                     # drop row id
    as.matrix()
  
  if(type == "column"){
    cors <- cor(resid_mat, use = "pairwise.complete.obs",
                method = "pearson")
    cors_long <- cors |>
      as.data.frame() |>
      rownames_to_column("row") |>
      pivot_longer(-row, names_to = "col", values_to = "corr") |> 
      mutate(col = factor(col, levels = 1:n_cols),
             row = factor(row, levels = n_cols:1))
    plot_label <- "Column Wise Residual Correlation"
    plot_axis <- element_text(angle = 45, hjust = 1)
  }
  if(type == "row"){
    cors <- cor(t(resid_mat), use = "pairwise.complete.obs",
                method = "pearson")
    cors_long <- cors |>
      as.data.frame() |>
      rownames_to_column("row") |>
      pivot_longer(-row, names_to = "col", values_to = "corr") |> 
      mutate(col = parse_number(col)) |> 
      mutate(col = factor(col, levels = 1:n_rows),
             row = factor(row, levels = n_rows:1))
    plot_label <- "Row Wise Residual Correlation"
    plot_axis <- element_blank()
  }
  
  
  
  ggplot(cors_long, aes(col, row, fill = corr)) +
    geom_tile(colour = "white") +
    coord_fixed() +
    scale_fill_gradient2(
      limits = c(-1, 1),          # full correlation range
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      name = "r"
    ) +
    labs(x = NULL, y = NULL,
         title = plot_label) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = plot_axis,
      axis.text.y = plot_axis,
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(4, "mm"),
      legend.position   = "right"
    )
}
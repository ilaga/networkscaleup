#' Permutation distribution of correlation between estimated alpha and residuals
#' from corresponding model
#' 
#' This function is incomplete, in progress
#'
#' This computes the correlation between the alphas and each column 
#' of the residual matrix, and finds the range of these correlations,
#' comparing the value to those under permutations of the alpha row
#'
#' @param alphas vector of the estimated
#' @param resids matrix of pearson residuals for ARD matrix
#' @param b number of permutations
#' @param plot whether to show the plot (default) or not
#'
#' @returns a list of test statistics
#' @export
degree_corr_perm <- function(alphas, resids, b = 1000,
                             plot = TRUE) {
  est_df <- cbind(alphas, resids)
  est_corr <- stats::cor(est_df)
  est_corr_range <- as.numeric(stats::dist(range(est_corr[1, -1])))
  
  perm_corr_range <- rep(NA, length = b)
  for (i in 1:b) {
    est_df[, 1] <- est_df[sample(1:nrow(est_df)), 1]
    tmp_corr <- stats::cor(est_df)
    perm_corr_range[i] <- as.numeric(stats::dist(range(tmp_corr[1, -1])))
  }
  
  if (plot) {
    graphics::hist(perm_corr_range,
         xlim = c(
           min(c(perm_corr_range, est_corr_range)),
           max(c(perm_corr_range, est_corr_range))
         ),
         xlab = c("Correlation Range"), main = "Histogram of Correlation Ranges"
    )
    graphics::abline(v = est_corr_range, col = "red")
  }
  return(list(
    test_stat = est_corr_range,
    bootstrap_stat = perm_corr_range
  ))
}

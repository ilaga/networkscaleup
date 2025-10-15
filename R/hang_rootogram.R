library(ggplot2)
library(dplyr)
library(scales)


#### TO DO
## add argument to cap the range of the plots
## add a label to the plot pointing out that it assumes poisson or negative
##


#' Hanging Rootogram for Fitted ARD Model
#'
#' @param y vector of ard entries
#' @param fit (posterior) means of each entry if poisson ARD model
#' @param family whether ard model poisson or negative binomial
#' @param size if negative binomial, size estimate
#' @param prb if negative binomial, prob estimate
#' @param width width of bars
#'
#' @return a ggplot of the hanging rootogram
#' @export
#'
#' @examples
hang_rootogram_ard <- function(y,
                               fit,
                               family = c("poisson", "nbinomial"),
                               # default
                               size = NULL,
                               # NB dispersion if needed
                               prob = NULL,
                               width = 0.9,
                               x_max = NULL) {
  # bar width (0–1)
  
  family <- match.arg(family)
  
  if (is.null(x_max)) {
    ## 1. support (integer counts)
    k <- 0:max(y, floor(max(y) * 1.5))
  } else {
    k <- 0:max(y, floor(x_max))
  }
  # generous upper limit
  
  if (family == "poisson") {
    lambda <- fit$summary(variables = "mu")$estimate
  } else if (family == "nbinomial") {
    prob <- rep(fit$summary(variables = "inv_omegas")$estimate, each = n_i)
    size <- fit$summary(variables = "par1")$estimate
  }
  
  
  ## 2. observed bin counts
  obs_counts <- as.numeric(table(factor(y, levels = k)))
  
  ## 3. expected bin counts
  exp_counts <- vapply(k, function(j) {
    if (family == "poisson") {
      sum(dpois(j, lambda = lambda))
    } else if (family == "nbinomial") {
      sum(dnbinom(j, size = size, prob = prob))
    }
  }, numeric(1))
  
  ## 4. square-root transform and bar coordinates
  obs_root <- sqrt(obs_counts)
  exp_root <- sqrt(exp_counts)
  baseline <- exp_root # hanging baseline, changed from -
  tips <- exp_root - obs_root # expected - observed
  
  df <- data.frame(
    k,
    xmin = k - width / 2,
    xmax = k + width / 2,
    ymin = pmin(baseline, tips),
    # bottom of bar
    ymax = pmax(baseline, tips) # top of bar
    # pos  = tips >= baseline         # TRUE = observed > expected
  ) |>
    dplyr::mutate(middle = (xmin + xmax) / 2)
  
  ## 5. plot
  ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      colour = "lightgray",
      fill = "gray"
    ) +
    geom_hline(yintercept = 0) +
    labs(x = "Count",
         y = expression(sqrt(count)),
         title = "Hanging Rootogram",) +
    theme_bw() +
    geom_line(aes(x = middle, y = ymax), col = "red") +
    geom_point(aes(x = middle, y = ymax), col = "red") +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = breaks_pretty(n = 6))
}

#' Scale raw log degree and log prevalence estimates
#'
#' This function scales estimates from either the overdispersed model or from
#' the correlated models. Several scaling options are available.
#'
#' @param log_degree The matrix of estimated raw log degrees from either the
#'   overdispersed or correlated models.
#' @param log_prevalence The matrix of estimates raw log prevalences from either
#'   the overdispersed or correlated models.
#' @param scaling An character vector providing the name of scaling procedure
#'   should be performed in order to transform estimates to degrees and
#'   subpopulation sizes. Scaling options are `overdispersed`, `all` (the
#'   default), `weighted`, or `weighted_sq` (`weighted` and `weighted_sq` are
#'   only available if `Correlation` is provided. Further details are provided
#'   in the Details section.
#' @param known_sizes The known subpopulation sizes corresponding to a subset of
#'   the columns of \code{ard}.
#' @param known_ind The indices that correspond to the columns of \code{ard}
#'   with known_sizes. By default, the function assumes the first \code{n_known}
#'   columns, where \code{n_known} corresponds to the number of
#'   \code{known_sizes}.
#' @param Correlation The estimated correlation matrix used to calculate scaling
#'   weights. Required if `scaling = weighted` or `scaling = weighted_sq`.
#' @param G1_ind If `scaling = overdispersed`, a vector of indices corresponding
#'   to the subpopulations that belong to the primary scaling groups, i.e. the
#'   collection of rare girls' names in Zheng, Salganik, and Gelman (2006). By
#'   default, all known_sizes are used. If G2_ind and B2_ind are not provided,
#'   `C = C_1`, so only G1_ind are used. If G1_ind is not provided, no scaling
#'   is performed.
#' @param G2_ind If `scaling = overdispersed`, a vector of indices corresponding
#'   to the subpopulations that belong to the first secondary scaling groups,
#'   i.e. the collection of somewhat popular girls' names.
#' @param B2_ind If `scaling = overdispersed`, a vector of indices corresponding
#'   to the subpopulations that belong to the second secondary scaling groups,
#'   i.e. the collection of somewhat popular boys' names.
#' @param N The known total population size.
#'
#' @details The `scaling` options are described below: \describe{\item{NULL}{No
#'   scaling is performed} \item{overdispersed}{The scaling procedure outlined
#'   in Zheng et al. (2006) is performed. In this case, at least `Pg1_ind` must
#'   be provided. See \link[networkscaleup]{overdispersedStan} for more
#'   details.} \item{all}{All subpopulations with known sizes are used to scale
#'   the parameters, using a modified scaling procedure that standardizes the
#'   sizes so each population is weighted equally. Additional details are
#'   provided in Laga et al. (2022+).} \item{weighted}{All subpopulations with
#'   known sizes are weighted according their correlation with the unknown
#'   subpopulation size. Additional details are provided in Laga et al. (2022+)}
#'   \item{weighted_sq}{Same as `weighted`, except the weights are squared,
#'   providing more relative weight to subpopulations with higher correlation.}}
#'
#' @return The named list containing the scaled log degree, degree, log
#'   prevalence, and size estimates
#' @references Zheng, T., Salganik, M. J., and Gelman, A. (2006). How many
#'   people do you know in prison, \emph{Journal of the American Statistical
#'   Association}, \bold{101:474}, 409--423
#'
#'   Laga, I., Bao, L., and Niu, X (2022+). A Correlated Network Scaleup Model:
#'   Finding the Connection Between Subpopulations, arxiv preprint:
#'   <https://doi.org/10.48550/arXiv.2109.10204>
#' @import rstan
#' @export
#'
#' @examples
#' # Analyze an example ard data set using Zheng et al. (2006) models
#' # Note that in practice, both warmup and iter should be much higher
#' data(example_data)
#'
#' overdisp.est = overdispersed(example_data$ard,
#' known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
#' known_ind = c(1, 3, 4),
#' G1_ind = 1,
#' G2_ind = 2,
#' B2_ind = 2,
#' N = example_data$N,
#' warmup = 250,
#' iter = 250)
#'
#' # Compare size estimates
#' data.frame(true = example_data$subpop_sizes,
#' basic = colMeans(overdisp.est$betas))
#'
#' # Compare degree estimates
#' plot(example_data$degrees, colMeans(overdisp.est$alphas))
#'
#' # Look at overdispersion parameter
#' colMeans(overdisp.est$omegas)
scaling <-
  function(log_degree,
           log_prevalence,
           scaling = c("all", "overdispersed", "weighted", "weighted_sq"),
           known_sizes = NULL,
           known_ind = NULL,
           Correlation = NULL,
           G1_ind = NULL,
           G2_ind = NULL,
           B2_ind = NULL,
           N = NULL) {
    ## Extract dimensions
    iter = nrow(log_degree) ## Number of MCMC samples
    N_i = ncol(log_degree) ## Number of respondents
    N_k = ncol(log_prevalence) ## Number of subpopulations

    scaling = match.arg(scaling)

    alphas = log_degree
    betas = log_prevalence



    if (is.null(refresh)) {
      ## By default, refresh every 10% of iterations
      refresh = round(iter / 10)
    }

    known_prevalences = known_sizes / N
    if (!is.null(G1_ind)) {
      Pg1 = sum(known_prevalences[G1_ind])
    }
    if (!is.null(G2_ind)) {
      Pg2 = sum(known_prevalences[G2_ind])
    }
    if (is.null(B2_ind)) {
      Pb2 = sum(known_prevalences[B2_ind])
    }


    if (scaling == "overdispersed") {
      if (is.null(G1_ind)) {
        stop("G1_ind cannot be null for scaling option \'overdispersed\'")
      }

      if (is.null(G2_ind) |
          is.null(B2_ind)) {
        ## Perform scaling with only main
        for (ind in 1:iter) {
          C1 = log(sum(exp(log_prevalence[ind, G1_ind]) / Pg1))
          C = C1

          alphas[ind, ] = alphas[ind, ] + C
          betas[ind, ] = betas[ind, ] - C
        }

      } else{
        ## Perform scaling with secondary groups
        for (ind in 1:iter) {
          C1 = log(sum(exp(log_prevalence[ind, G1_ind]) / Pg1))
          C2 = log(sum(exp(log_prevalence[ind, B2_ind]) / Pb2)) -
            log(sum(exp(log_prevalence[ind, G2_ind]) / Pg2))
          C = C1 + 1 / 2 * C2

          alphas[ind, ] = alphas[ind, ] + C
          betas[ind, ] = betas[ind, ] - C
        }

      }
    } else if (scaling == "all") {
      for (ind in 1:iter) {
        C = log(mean(exp(log_prevalence[ind, known_ind]) / known_prevalences[known_ind]))

        alphas[ind,] = alphas[ind,] + C
        betas[ind,] = betas[ind,] - C
      }
    } else if (scaling == "weighted") {
      if (is.null(Correlation)) {
        stop("Correlation cannot be null for scaling option \'weighted\'")
      }

      for (k in 1:N_k) {
        scale.weights = Correlation[k, known_ind]
        ## Set negative weights to 0
        scale.weights[scale.weights < 0] = 0
        scale.weights = scale.weights / sum(scale.weights, na.rm = T) * sum(!is.na(scale.weights))
        for (ind in 1:iter) {
          C = log(mean(
            exp(log_prevalence[ind, known_ind]) * scale.weights / known_prevalences[known_ind],
            na.rm = T
          ))

          betas[ind, k] = betas[ind, k] - C
        }
        print(k)
      }

      ## Scale degrees separately using all
      for (ind in 1:iter) {
        C = log(mean(exp(log_prevalence[ind, known_ind]) / known_prevalences[known_ind]))

        alphas[ind,] = alphas[ind,] + C
      }


    } else if (scaling == "weighted_sq") {
      if (is.null(Correlation)) {
        stop("Correlation cannot be null for scaling option \'weighted_sq\'")
      }

      for (k in 1:N_k) {
        scale.weights = Correlation[k, known_ind]
        ## Set negative weights to 0
        scale.weights[scale.weights < 0] = 0
        scale.weights = scale.weights ^ 2
        scale.weights = scale.weights / sum(scale.weights, na.rm = T) * sum(!is.na(scale.weights))
        for (ind in 1:iter) {
          C = log(mean(
            exp(log_prevalence[ind, known_ind]) * scale.weights / known_prevalences[known_ind],
            na.rm = T
          ))

          betas[ind, k] = betas[ind, k] - C
        }
        print(k)
      }

      ## Scale degrees separately using all
      for (ind in 1:iter) {
        C = log(mean(exp(log_prevalence[i, known_ind]) / known_prevalences[known_ind]))

        alphas[ind,] = alphas[ind,] + C
      }
    }


    return_list = list(log_degree = alphas,
                       log_prevalence = betas)

    return(return_list)


  }

#' Fit Killworth models to ARD. This function estimates the degrees and
#' population sizes using the plug-in MLE and MLE estimator.
#'
#' @param ard The `n_i x n_k` matrix of non-negative ARD integer responses,
#'   where the `(i,k)th` element corresponds to the number of people that
#'   respondent `i` knows in subpopulation `k`.
#' @param known_sizes The known subpopulation sizes corresponding to a subset of
#'   the columns of \code{ard}.
#' @param known_ind The indices that correspond to the columns of \code{ard}
#'   with known_sizes. By default, the function assumes the first \code{n_known}
#'   columns, where \code{n_known} corresponds to the number of
#'   \code{known_sizes}.
#' @param N The known total population size.
#' @param model A character string corresponding to either the plug-in MLE
#'   (PIMLE) or the MLE (MLE). The function assumes MLE by default.
#'
#' @return A named list with the estimated degrees and sizes.
#' @references Killworth, P. D., Johnsen, E. C., McCarty, C., Shelley, G. A.,
#'   and Bernard, H. R. (1998). A Social Network Approach to Estimating
#'   Seroprevalence in the United States, \emph{Social Networks}, \bold{20},
#'   23--50
#'
#'   Killworth, P. D., McCarty, C., Bernard, H. R., Shelley, G. A., and Johnsen,
#'   E. C. (1998). Estimation of Seroprevalence, Rape and Homelessness in the
#'   United States Using a Social Network Approach, \emph{Evaluation Review},
#'   \bold{22}, 289--308
#'
#'   Laga, I., Bao, L., and Niu, X. (2021). Thirty Years of the Network Scale-up
#'   Method, \emph{Journal of the American Statistical Association},
#'   \bold{116:535}, 1548--1559
#' @export
#'
#' @examples
#' # Analyze an example ard data set using the killworth function
#' data(example_data)
#'
#' ard = example_data$ard
#' subpop_sizes = example_data$subpop_sizes
#' N = example_data$N
#'
#' mle.est = killworth(ard,
#' known_sizes = subpop_sizes[c(1, 2, 4)],
#' known_ind = c(1, 2, 4),
#' N = N, model = "MLE")
#'
#' pimle.est = killworth(ard,
#' known_sizes = subpop_sizes[c(1, 2, 4)],
#' known_ind = c(1, 2, 4),
#' N = N, model = "PIMLE")
#'
#' ## Compare estimates with the truth
#' plot(mle.est$degrees, example_data$degrees)
#'
#' data.frame(true = subpop_sizes[c(3, 5)],
#' mle = mle.est$sizes,
#' pimle = pimle.est$sizes)
killworth <-
  function(ard,
           known_sizes = NULL,
           known_ind = 1:length(known_sizes),
           N = NULL,
           model = c("MLE", "PIMLE")) {
    ## Extract dimensions
    N_i = nrow(ard) ## Number of respondents
    N_k = ncol(ard) ## Number of subpopulations


    n_known = length(known_sizes)
    n_unknown = N_k - n_known
    unknown_ind = (1:N_k)[-known_ind]


    if (!(model %in% c("MLE", "PIMLE"))) {
      stop("model is not one of \'MLE\' or \'PIMLE\'")
    }

    if (length(known_sizes) != length(known_ind)) {
      stop("known_sizes and known_ind must the same length")
    }

    ## Estimate degrees
    d.est = N * rowSums(ard[, known_ind]) / sum(known_sizes)


    if (n_unknown == 1) {
      ## If only 1 unknown subpopulation
      if (model == "MLE") {
        N.est = N * sum(ard[, unknown_ind]) / sum(d.est)
      } else{
        if (min(d.est) == 0) {
          pos.ind = which(d.est > 0)
          warning("Estimated a 0 degree for at least one respondent. Ignoring response for PIMLE")
        } else{
          pos.ind = 1:N_i
        }
        N.est = N * mean(ard[pos.ind, unknown_ind] / d.est[pos.ind])
      }
    } else{
      ## If multiple unknown subpopulations
      N.est = rep(NA, n_unknown)
      if (model == "MLE") {
        for (k in 1:n_unknown) {
          N.est[k] = N * sum(ard[, unknown_ind[k]]) / sum(d.est)
        }
      } else{
        if (min(d.est) == 0) {
          pos.ind = which(d.est > 0)
          warning("Estimated a 0 degree for at least one respondent. Ignoring response for PIMLE")
        } else{
          pos.ind = 1:N_i
        }
        for (k in 1:n_unknown) {
          N.est[k] = N * mean(ard[pos.ind, unknown_ind[k]] / d.est[pos.ind])
        }
      }
    }

    return.list = list(degrees = d.est,
                       sizes = N.est)

    return(return.list)
  }

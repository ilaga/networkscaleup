#' Fit ARD using the transmission model

#' This function fits the ARD using the transmission and accompanying models
#' proposed in Teo et al. (2019) in Stan.
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
#' @param model A character vector denoting which of the three models should be
#'   fit, either 'basic', 'transmission', and 'barrier'. More details of these
#'   models are provided below.
#' @param x A matrix with dimensions `n_i x n_unknown`, where `n_unknown` refers
#'   to the number of unknown subpopulation sizes. In the language of Teo et al.
#'   (2019), these represent the individual's perception of each hidden
#'   population.
#' @param z A matrix with dimensions `n_i x p`, where `p` is the number of
#'   demographic covariates used. This matrix represents the demographic
#'   information about the respondents in order to capture the barrier effects.
#' @param chains A positive integer specifying the number of Markov chains.
#' @param cores A positive integer specifying the number of cores to use to run
#'   the Markov chains in parallel.
#' @param warmup A positive integer specifying the total number of samples for
#'   each chain (including warmup). Matches the usage in \link[rstan]{stan}.
#' @param iter A positive integer specifying the number of warmup samples for
#'   each chain. Matches the usage in \link[rstan]{stan}.
#' @param thin A positive integer specifying the interval for saving posterior
#'   samples. Default value is 1 (i.e. no thinning).
#' @param return_fit A logical indicating whether the fitted Stan model should
#'   be returned instead of the \code{rstan::extracted} and scaled parameters.
#'   This is FALSE by default.
#' @param ... Additional arguments to be passed to \link[rstan]{stan}.
#'
#' @details This function fits either the basic (\code{"basic"}), transmission
#'   (\code{"transmission"}), or transmission and barrier effects
#'   (\code{"barrier"}) models proposed in Teo et al. (2019). The basic model
#'   does not include any covariates or social acceptability rating (\code{z}
#'   and \code{x}, respectively). The transmission model includes the social
#'   acceptability through \code{z}. The barrier model includes both the social
#'   acceptability and the participants' demographics through \code{x}. The
#'   models are coded in Stan and fit using RStan.
#'
#' @return Either the full fitted Stan model if \code{return_fit = TRUE}, else a
#'   named list with the estimated parameters extracted using
#'   \link[rstan]{extract} (the default). The estimated parameters are named as
#'   follows, with additional descriptions as needed:
#'
#'   \describe{\item{lambda}{} \item{alpha}{} \item{tau}{Standard deviation of
#'   alpha} \item{S_H}{Estimated hidden subpopulation sizes}}
#' @references Teo, A.K.J., Prem, K, Chen, M.I.C., Roellin, A, Wong, M.L., La,
#'   H.H., Cook, A. (2019). Estimating the size of key populations for HIV in
#'   Singapore using the network scale-up method, \emph{Epidemiology},
#'   \bold{95}, 602--607
#' @export
#'
#' @examples
#' # Analyze an example ard data set using Teo et al. (2019) models
#' # Note that in practice, both warmup and iter should be much higher
#' data(example_data)
#'
#' basic.est = transmissionStan(example_data$ard,
#' known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
#' known_ind = c(1, 2, 4),
#' N = example_data$N,
#' model = "basic",
#' chains = 1,
#' cores = 1,
#' warmup = 250,
#' iter = 500)
#'
#' transmission.est = transmissionStan(example_data$ard,
#' known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
#' known_ind = c(1, 2, 4),
#' N = example_data$N,
#' model = "transmission",
#' x = x,
#' chains = 1,
#' cores = 1,
#' warmup = 250,
#' iter = 500)
#'
#' barrier.est = transmissionStan(example_data$ard,
#' known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
#' known_ind = c(1, 2, 4),
#' N = example_data$N,
#' model = "barrier",
#' x = x,
#' z = z,
#' chains = 1,
#' cores = 1,
#' warmup = 250,
#' iter = 500)
#'
#' # Compare size estimates
#' data.frame(true = example_data$subpop_sizes[c(3, 5)],
#' basic = colMeans(basic.est$S_H),
#' transmission = colMeans(transmission.est$S_H),
#' barrier = colMeans(barrier.est$S_H))
#'
#' # Look at z slope parameters
#' colMeans(transmission.est$beta)
#'
#' # Compare to z slope parameters from the barrier effect model
#' colMeans(barrier.est$beta)
#'
#' # Look at x slope parameters
#' colMeans(barrier.est$gamma)
transmissionStan <-
  function(ard,
           known_sizes = NULL,
           known_ind = NULL,
           N = NULL,
           model = c("basic", "transmission", "barrier"),
           x = NULL,
           z = NULL,
           chains = 3,
           cores = 1,
           warmup = 1000,
           iter = 1500,
           thin = 1,
           return_fit = FALSE,
           ...) {
    N_i = nrow(ard)
    N_k = ncol(ard)

    ## Check dimensions of x
    if (is.null(x)) {
      if ((nrow(x) != N_i) | (ncol(x) != N_k)) {
        stop("Dimensions of x do not match dimensions of ard")
      }
    }

    ## Check dimensions of z
    if (is.null(z)) {
      if ((nrow(z) != N_i)) {
        stop("Dimensions of z do not match dimensions of ard")
      } else{
        z_size = ncol(z)
      }
    }

    model = match.arg(model)

    ## reorder ard, if necessary (i.e. not in order of known then unknown)
    unknown_ind = (1:N_k)[-known_ind]
    ard = ard[, c(known_ind, unknown_ind)]

    n_known = length(known_ind)
    n_unknown = N_k - n_known



    if (model == "basic") {
      stan_data = list(
        N = N,
        n_i = N_i,
        n_k = N_k,
        n_known = n_known,
        n_unknown = n_unknown,
        S_k = known_sizes,
        y = ard
      )

      ## Fit model
      transmission_fit = rstan::sampling(
        object = stanmodels$Teo_basic_model,
        data = stan_data,
        chains = chains,
        cores = cores,
        iter = iter,
        warmup = warmup,
        refresh = 100,
        thin = thin,
        ...
      )

    } else if (model == "transmission") {
      stan_data = list(
        N = N,
        n_i = N_i,
        n_k = N_k,
        n_known = n_known,
        n_unknown = n_unknown,
        x = x,
        S_k = known_sizes,
        y = ard
      )



      ## Fit model
      transmission_fit = rstan::sampling(
        object = stanmodels$Teo_transmission_model,
        data = stan_data,
        chains = chains,
        cores = cores,
        iter = iter,
        warmup = warmup,
        refresh = 100,
        thin = thin,
        ...
      )

    } else if (model == "barrier") {
      stan_data = list(
        N = N,
        n_i = N_i,
        n_k = N_k,
        n_known = n_known,
        n_unknown = n_unknown,
        z_size = z_size,
        x = x,
        z = z,
        S_k = known_sizes,
        y = ard
      )


      ## Fit model
      transmission_fit = rstan::sampling(
        object = stanmodels$Teo_barrier_model,
        data = stan_data,
        chains = chains,
        cores = cores,
        iter = iter,
        warmup = warmup,
        refresh = 100,
        thin = thin,
        ...
      )
    } else{
      stop(
        "A valid model argument was not provided. Must be one of \'basic\', \'transmission\', or \'barrier\'."
      )
    }

    draws = rstan::extract(transmission_fit)

    ## Return values

    if (return_fit) {
      return(transmission_fit)
    } else{
      return(draws)
    }

  }

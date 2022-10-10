#' Fit ARD using the uncorrelated or correlated model in Stan

#' This function fits the ARD using either the uncorrelated or correlated model
#' in Laga et al. (2021) in Stan. The population size estimates and degrees are
#' scaled using a post-hoc procedure.
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
#' @param model A character vector denoting which of the two models should be
#'   fit, either 'uncorrelated' or 'correlated'. More details of these models
#'   are provided below. The function decides which covariate model is needed
#'   based on the covariates provided below.
#' @param scaling An optional character vector providing the name of scaling
#'   procedure should be performed in order to transform estimates to degrees
#'   and subpopulation sizes. If `NULL`, the parameters will be returned
#'   unscaled. Alternatively, scaling may be performed independently using the
#'   \link[networkscaleup]{scaling} function. Scaling options are `NULL`,
#'   `overdispersed`, `all`, `weighted`, or `weighted_sq` (`weighted` and
#'   `weighted_sq` are only available if `model = "correlated"`. Further details
#'   are provided in the Details section.
#' @param x A matrix with dimensions `n_i x n_unknown`, where `n_unknown` refers
#'   to the number of unknown subpopulation sizes. In the language of Teo et al.
#'   (2019), these represent the individual's perception of each hidden
#'   population.
#' @param z_global A matrix with dimensions `n_i x p_global`, where `p_global`
#'   is the number of demographic covariates used. This matrix represents the
#'   demographic information about the respondents in order to capture the
#'   barrier effects.
#' @param z_subpop A matrix with dimensions `n_i x p_subpop`, where `p_subpop`
#'   is the number of demographic covariates used. This matrix represents the
#'   demographic information about the respondents in order to capture the
#'   barrier effects.
#' @param G1_ind A vector of indices denoting the columns of `ard` that
#'   correspond to the primary scaling groups, i.e. the collection of rare
#'   girls' names in Zheng, Salganik, and Gelman (2006). By default, all
#'   known_sizes are used. If G2_ind and B2_ind are not provided, `C = C_1`, so
#'   only G1_ind are used. If G1_ind is not provided, no scaling is performed.
#' @param G2_ind A vector of indices denoting the columns of `ard` that
#'   correspond to the subpopulations that belong to the first secondary scaling
#'   groups, i.e. the collection of somewhat popular girls' names.
#' @param B2_ind A vector of indices denoting the columns of `ard` that
#'   correspond to the subpopulations that belong to the second secondary
#'   scaling groups, i.e. the collection of somewhat popular boys' names.
#' @param chains A positive integer specifying the number of Markov chains.
#' @param cores A positive integer specifying the number of cores to use to run
#'   the Markov chains in parallel.
#' @param warmup A positive integer specifying the total number of samples for
#'   each chain (including warmup). Matches the usage in \link[rstan]{stan}.
#' @param iter A positive integer specifying the number of warmup samples for
#'   each chain. Matches the usage in \link[rstan]{stan}.
#' @param thin A positive integer specifying the interval for saving posterior
#'   samples. Default value is 1 (i.e. no thinning).
#' @param return_fit A logical indicating whether the fitted `stanfit` object
#'   should be return. Defaults to `FALSE`.
#' @param ... Additional arguments to be passed to \link[rstan]{stan}.
#'
#' @details This function currently fits a variety of models proposed in Laga et
#'   al. (2022+). The user may provide any combination of `x`, `z_global`, and
#'   `z_subpop`. Additionally, the user may choose to fit a uncorrelated version
#'   of the model, where the correlation matrix is equal to the identity matrix.
#'
#'   The `scaling` options are described below: \describe{\item{NULL}{No scaling
#'   is performed} \item{overdispersed}{The scaling procedure outlined in Zheng
#'   et al. (2006) is performed. In this case, at least `Pg1_ind` must be
#'   provided. See \link[networkscaleup]{overdispersedStan} for more details.}
#'   \item{all}{All subpopulations with known sizes are used to scale the
#'   parameters, using a modified scaling procedure that standardizes the sizes
#'   so each population is weighted equally. Additional details are provided in
#'   Laga et al. (2022+).} \item{weighted}{All subpopulations with known sizes
#'   are weighted according their correlation with the unknown subpopulation
#'   size. Additional details are provided in Laga et al. (2022+)}
#'   \item{weighted_sq}{Same as `weighted`, except the weights are squared,
#'   providing more relative weight to subpopulations with higher correlation.}}
#'
#' @return Either the full fitted Stan model if \code{return_fit = TRUE}, else a
#'   named list with the estimated parameters extracted using
#'   \link[rstan]{extract} (the default). The estimated parameters are named as
#'   follows (if estimated in the corresponding model), with additional
#'   descriptions as needed:
#'
#'   \describe{\item{delta}{Raw delta parameters} \item{sigma_delta}{Standard
#'   deviation of delta} \item{rho}{Log prevalence, if scaled, else raw rho
#'   parameters} \item{mu_rho}{Mean of rho} \item{sigma_rho}{Standard deviation
#'   of rho} \item{alpha}{Slope parameters corresponding to z}
#'   \item{beta_global}{Slope parameters corresponding to x_global}
#'   \item{beta_subpop}{Slope parameters corresponding to x_subpop}
#'   \item{tau_N}{Standard deviation of random effects b}
#'   \item{Corr}{Correlation matrix, if `Correlation = TRUE`}}
#'
#'   If scaled, the following additional parameters are included:
#'   \describe{\item{log_degrees}{Scaled log degrees} \item{degree}{Scaled
#'   degrees} \item{log_prevalences}{Scaled log prevalences}
#'   \item{sizes}{Subpopulation size estimates}}
#' @references Laga, I., Bao, L., and Niu, X (2021). A Correlated Network
#'   Scaleup Model: Finding the Connection Between Subpopulations
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_data)
#'
#' x = example_data$x
#' z_global = example_data$z[,1:2]
#' z_subpop = example_data$z[,3:4]
#'
#' basic_corr_est = correlatedStan(example_data$ard,
#'      known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
#'      known_ind = c(1, 2, 4),
#'      N = example_data$N,
#'      model = "correlated",
#'      scaling = "weighted",
#'      chains = 1,
#'      cores = 1,
#'      warmup = 50,
#'      iter = 100)
#'
#' cov_uncorr_est = correlatedStan(example_data$ard,
#'      known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
#'      known_ind = c(1, 2, 4),
#'      N = example_data$N,
#'      model = "uncorrelated",
#'      scaling = "all",
#'      x = x,
#'      z_global = z_global,
#'      z_subpop = z_subpop,
#'      chains = 1,
#'      cores = 1,
#'      warmup = 50,
#'      iter = 100)
#'
#' cov_corr_est = correlatedStan(example_data$ard,
#'      known_sizes = example_data$subpop_sizes[c(1, 2, 4)],
#'      known_ind = c(1, 2, 4),
#'      N = example_data$N,
#'      model = "correlated",
#'      scaling = "all",
#'      x = x,
#'      z_subpop = z_subpop,
#'      chains = 1,
#'      cores = 1,
#'      warmup = 50,
#'      iter = 100)
#'
#' # Compare size estimates
#' round(data.frame(true = example_data$subpop_sizes,
#'      corr_basic = colMeans(basic_corr_est$sizes),
#'      uncorr_x_zsubpop_zglobal = colMeans(cov_uncorr_est$sizes),
#'      corr_x_zsubpop = colMeans(cov_corr_est$sizes)))
#'
#' # Look at z slope parameters
#' colMeans(cov_uncorr_est$beta_global)
#' colMeans(cov_corr_est$beta_subpop)
#' colMeans(cov_uncorr_est$beta_subpop)
#'
#' # Look at x slope parameters
#' colMeans(cov_uncorr_est$alpha)
#' colMeans(cov_corr_est$alpha)
#' }
correlatedStan <-
  function(ard,
           known_sizes = NULL,
           known_ind = NULL,
           N = NULL,
           model = c("correlated", "uncorrelated"),
           scaling = c("all", "overdispersed", "weighted", "weighted_sq"),
           x = NULL,
           z_global = NULL,
           z_subpop = NULL,
           G1_ind = NULL,
           G2_ind = NULL,
           B2_ind = NULL,
           chains = 3,
           cores = 1,
           warmup = 1000,
           iter = 1500,
           thin = 1,
           return_fit = FALSE,
           ...) {
    N_i = nrow(ard)
    N_k = ncol(ard)

    model = match.arg(model)
    scaling = match.arg(scaling)

    ## Check dimensions of x
    if (!is.null(x)) {
      if ((nrow(x) != N_i) | (ncol(x) != N_k)) {
        stop("Dimensions of x do not match dimensions of ard")
      }
    }

    ## Check for scaling method
    if (model == "uncorrelated" &
        (scaling == "weighted" | scaling == "weighted_sq")) {
      stop("Model must be `correlated` to using `weighted` or `weighted_sq` scaling")
    }

    ## Check dimensions of z
    if (!is.null(z_global)) {
      if ((nrow(z_global) != N_i)) {
        stop("Dimensions of z_global do not match dimensions of ard")
      } else{
        z_global_size = ncol(z_global)
      }
    }

    if (!is.null(z_subpop)) {
      if ((nrow(z_subpop) != N_i)) {
        stop("Dimensions of z_subpop do not match dimensions of ard")
      } else{
        z_subpop_size = ncol(z_subpop)
      }
    }


    ## Set model of 16 possible combinations
    if (model == "correlated") {
      if (is.null(x)) {
        ## No level of respect
        if (is.null(z_global)) {
          ## No global covariates
          if (is.null(z_subpop)) {
            ## No subpop covariates
            ## Basic model
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Correlated_basic,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          } else{
            ## Includes only zsubpop
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              z_subpop_size = z_subpop_size,
              z_subpop = z_subpop,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Correlated_zsubpop,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          }

        } else{
          ## Does include global covariates
          if (is.null(z_subpop)) {
            ## No subpop covariates
            ## Only global
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              z_global_size = z_global_size,
              z_global = z_global,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Correlated_zglobal,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          } else{
            ## Includes zglobal and zsubpop
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              z_global_size = z_global_size,
              z_global = z_global,
              z_subpop_size = z_subpop_size,
              z_subpop = z_subpop,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Correlated_zsubpop_zglobal,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          }

        }
      } else{
        ## Includes x
        if (is.null(z_global)) {
          ## No global covariates
          if (is.null(z_subpop)) {
            ## No subpop covariates
            ## Only x
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              x = x,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Correlated_x,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          } else{
            ## Includes x and z_subpop
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              x = x,
              z_subpop_size = z_subpop_size,
              z_subpop = z_subpop,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Correlated_x_zsubpop,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          }

        } else{
          ## Does include global covariates
          if (is.null(z_subpop)) {
            ## No subpop covariates
            ## Includes x and zglobal
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              x = x,
              z_global_size = z_global_size,
              z_global = z_global,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Correlated_x_zglobal,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          } else{
            ## Includes x, zglobal, and zsubpop
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              z_global_size = z_global_size,
              z_global = z_global,
              z_subpop_size = z_subpop_size,
              z_subpop = z_subpop,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Correlated_x_zsubpop_zglobal,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          }

        }
      }
    } else if (model == "uncorrelated") {
      if (is.null(x)) {
        ## No level of respect
        if (is.null(z_global)) {
          ## No global covariates
          if (is.null(z_subpop)) {
            ## No subpop covariates
            ## Basic model
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Uncorrelated_basic,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          } else{
            ## Includes only zsubpop
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              z_subpop_size = z_subpop_size,
              z_subpop = z_subpop,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Uncorrelated_zsubpop,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          }

        } else{
          ## Does include global covariates
          if (is.null(z_subpop)) {
            ## No subpop covariates
            ## Only global
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              z_global_size = z_global_size,
              z_global = z_global,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Uncorrelated_zglobal,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          } else{
            ## Includes zglobal and zsubpop
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              z_global_size = z_global_size,
              z_global = z_global,
              z_subpop_size = z_subpop_size,
              z_subpop = z_subpop,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Uncorrelated_zsubpop_zglobal,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          }

        }
      } else{
        ## Includes x
        if (is.null(z_global)) {
          ## No global covariates
          if (is.null(z_subpop)) {
            ## No subpop covariates
            ## Only x
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              x = x,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Uncorrelated_x,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          } else{
            ## Includes x and z_subpop
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              x = x,
              z_subpop_size = z_subpop_size,
              z_subpop = z_subpop,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Uncorrelated_x_zsubpop,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          }

        } else{
          ## Does include global covariates
          if (is.null(z_subpop)) {
            ## No subpop covariates
            ## Includes x and zglobal
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              x = x,
              z_global_size = z_global_size,
              z_global = z_global,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Uncorrelated_x_zglobal,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          } else{
            ## Includes x, zglobal, and zsubpop
            stan_data = list(
              N = N,
              n_i = N_i,
              n_k = N_k,
              z_global_size = z_global_size,
              z_global = z_global,
              z_subpop_size = z_subpop_size,
              z_subpop = z_subpop,
              y = ard
            )
            ## Fit model
            model_fit = rstan::sampling(
              object = stanmodels$Uncorrelated_x_zsubpop_zglobal,
              data = stan_data,
              chains = chains,
              cores = cores,
              iter = iter,
              warmup = warmup,
              thin = thin,
              ...
            )
          }

        }
      }
    } else{
      stop("Invalid model choice")
    }


    ## Extract draws
    ## Exclude eps and L_Omega  (if correlated) for memory
    if (model == "correlated") {
      draws = rstan::extract(model_fit,
                             pars = c("eps", "L_Omega"),
                             include = FALSE)
    } else{
      draws = rstan::extract(model_fit, pars = c("eps"), include = FALSE)
    }


    ## Perform scaling procedure
    delta = draws$delta
    sigma_delta = draws$sigma_delta
    log_degrees = matrix(NA, nrow = nrow(delta), ncol = ncol(delta))
    for (i in 1:nrow(log_degrees)) {
      log_degrees[i, ] = delta[i, ] * sigma_delta[i]
    }

    if (!is.null(scaling)) {
      if ((scaling == "weighted") | (scaling == "weighted_sq")) {
        ## First get point estimate for correlation matrix
        Correlation = draws$Corr
        Correlation = apply(Correlation, c(2, 3), mean)
        scaling_res = networkscaleup::scaling(
          log_degrees,
          draws$rho,
          scaling = scaling,
          known_sizes = known_sizes,
          known_ind = known_ind,
          Correlation = Correlation,
          N = N
        )


      } else if (scaling == "all") {
        scaling_res = networkscaleup::scaling(
          log_degrees,
          draws$rho,
          scaling = scaling,
          known_sizes = known_sizes,
          known_ind = known_ind,
          N = N
        )

      } else if (scaling == "overdispersed") {
        scaling_res = networkscaleup::scaling(
          log_degrees,
          draws$rho,
          scaling = scaling,
          known_sizes = known_sizes,
          known_ind = known_ind,
          G1_ind = G1_ind,
          G2_ind = G2_ind,
          B2_ind = B2_ind,
          N = N
        )
      }

      draws$log_degrees = scaling_res$log_degrees
      draws$degrees = exp(draws$log_degrees)
      draws$log_prevalences = scaling_res$log_prevalences
      draws$sizes = exp(draws$log_prevalences) * N
    }





    ## Return values

    if (return_fit) {
      return(model_fit)
    } else{
      return(draws)
    }

  }

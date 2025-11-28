#' 
#' library(cmdstanr)
#' library(here)
#' 
#' 
#' #' Find MAP estimates of basic Poisson and Negative Binomial models using Stan optimization
#' #'
#' #' @param ard n_i by n_k ARD matrix
#' #' @param x_cov n_i by p covariate matrix with columns names
#' #' @param cov_names Vector of column names to include
#' #' @param family Distribution to fit, either "poisson" or "nbinomial"
#' #'
#' #' @return Stan fit
#' #' @export
#' 
#' fit_stan_optim = function(ard,
#'                           x_cov = NULL,
#'                           cov_names = NULL,
#'                           family = c("poisson", "nbinomial"),
#'                           ...) {
#'   ## Grab family
#'   family = match.arg(family, c("poisson", "nbinomial"))
#'   
#'   ## Only does subpop covariates at the moment
#'   if(is.null(x_cov)){
#'     stan_data = list(
#'       y = ard,
#'       n_i = nrow(ard),
#'       n_k = ncol(ard)
#'     )
#'     if(family == "poisson"){
#'       mod = cmdstan_model(here("Stan_Files/Poisson.stan"))
#'     }else{
#'       mod = cmdstan_model(here("Stan_Files/Overdispersed.stan"))
#'     }
#'   }else{
#'     stan_data = list(
#'       y = ard,
#'       n_i = nrow(ard),
#'       n_k = ncol(ard),
#'       z_subpop_size = ncol(x_cov), #length(cov_names),
#'       z_subpop = x_cov#[,cov_names]
#'     )
#'     if(family == "poisson"){
#'       mod = cmdstan_model(here("Stan_Files/Poisson_zsubpop.stan"))
#'     }else{
#'       mod = cmdstan_model(here("Stan_Files/Overdispersed_zsubpop.stan"))
#'     }
#'   }
#'   
#'   fit = mod$optimize(data = stan_data, ...)
#'   
#' }

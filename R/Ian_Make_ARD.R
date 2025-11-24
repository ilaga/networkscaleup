# library(trialr)
# 
# 
# #1. Poisson
# make_Poisson <- function(n_i = 500,
#                          n_k = 20,
#                          N = 1000000,
#                          p = 9,
#                          p_nonzero = 3,
#                          seed = NULL) {
#   if (!is.null(seed)) { # Only set seed if provided
#     set.seed(seed)
#   }
#   nk_prev = runif(n_k, 0.01, 0.15)
#   nk_size = round(nk_prev * N)
#   d = round(1 + rlnorm(n_i, 5, 0.4))
#   betas = log(nk_size / N)
#   alphas = log(d)
#   ard = matrix(NA, nrow = n_i, ncol = n_k)
#   if (p == 0) {
#     for (k in 1:n_k) {
#       ard[, k] = rpois(n_i, lambda = exp(alphas + betas[k]))
#     }
#   } else{
#     x_cov = matrix(runif(n_i * p), nrow = n_i, ncol = p)
#     x_beta = matrix(0, nrow = p, ncol = n_k)
#     x_beta[sample(p, p_nonzero), ] = runif(p_nonzero * n_k, -2, 2)
#     for (k in 1:n_k) {
#       ard[, k] = rpois(n_i, lambda = exp(alphas + betas[k] + x_cov %*% x_beta[, k]))
#     }
#   }
#   #The return gives all the information need to recreate this data in and out of function.
#   if (p == 0) {
#     return(
#       list(
#         ard = ard,
#         prev = nk_prev,
#         size = nk_size,
#         d = d,
#         alphas = alphas,
#         betas = betas,
#         n_i = n_i,
#         n_k = n_k,
#         seed = seed
#       )
#     )
#   } else{
#     return(
#       list(
#         ard = ard,
#         prev = nk_prev,
#         size = nk_size,
#         d = d,
#         x_cov = x_cov,
#         x_beta = x_beta,
#         n_i = n_i,
#         n_k = n_k,
#         seed = seed
#       )
#     )
#   }
#   
# }
# 
# 
# #2. Negative Binomial
# 
# make_NB <- function(n_i = 500,
#                     n_k = 20,
#                     N = 1000000,
#                     p = 9,
#                     p_nonzero = 3,
#                     omega_range = c(1, 5),
#                     seed = NULL) {
#   if (!is.null(seed)) { # Only set seed if provided
#     set.seed(seed)
#   }
#   nk_prev = runif(n_k, 0.01, 0.15)
#   nk_size = round(nk_prev * N)
#   d = round(1 + rlnorm(n_i, 5, 0.4))
#   omega = runif(n_k, omega_range[1], omega_range[2])
#   betas = log(nk_size / N)
#   alphas = log(d)
#   ard = matrix(NA, nrow = n_i, ncol = n_k)
#   if (p == 0) {
#     for (k in 1:n_k) {
#       ard[, k] = rnbinom(n_i,
#                          size = exp(alphas + betas[k]) / (omega[k] - 1),
#                          prob = 1 / omega[k])
#     }
#   } else{
#     x_cov = matrix(runif(n_i * p), nrow = n_i, ncol = p)
#     x_beta = matrix(0, nrow = p, ncol = n_k)
#     x_beta[sample(p, p_nonzero), ] = runif(p_nonzero * n_k, -2, 2)
#     for (k in 1:n_k) {
#       ard[, k] = rnbinom(
#         n_i,
#         size = exp(alphas + betas[k] + x_cov %*% x_beta[, k]) / (omega[k] - 1),
#         prob = 1 / omega[k]
#       )
#     }
#   }
#   if (p == 0) {
#     return(
#       list(
#         ard = ard,
#         prev = nk_prev,
#         size = nk_size,
#         omega = omega,
#         d = d,
#         n_i = n_i,
#         n_k = n_k,
#         seed = seed
#       )
#     )
#   } else{
#     return(
#       list(
#         ard = ard,
#         prev = nk_prev,
#         size = nk_size,
#         omega = omega,
#         d = d,
#         alphas = alphas,
#         betas = betas,
#         x_cov = x_cov,
#         x_beta = x_beta,
#         n_i = n_i,
#         n_k = n_k,
#         seed = seed
#       )
#     )
#   }
# }
# 
# # 3. Correlated
# make_Correlated <- function(n_i = 500,
#                             n_k = 20,
#                             N = 1000000,
#                             eta = 2,
#                             p = 9,
#                             p_nonzero = 3,
#                             seed = NULL) {
#   if (!is.null(seed)) { # Only set seed if provided
#     set.seed(seed)
#   }
#   nk_prev = runif(n_k, 0.01, 0.15)
#   nk_size = round(nk_prev * N)
#   d = round(1 + rlnorm(n_i, 5, 0.4))
#   tau.N = runif(n_k, min = 0.5, max = 1.5)
#   mu = log(1 / sqrt(1 + tau.N ^ 2))
#   tau = sqrt(log(1 + tau.N ^ 2))
#   Omega = trialr::rlkjcorr(1, n_k, eta = eta)
#   L.omega = t(chol(Omega))
#   eps = matrix(rnorm(n_i * n_k), nrow = n_i, ncol = n_k)
#   bias = matrix(NA, nrow = n_i, ncol = n_k)
#   ard = matrix(NA, nrow = n_i, ncol = n_k)
#   for (i in 1:n_i) {
#     bias[i,] = mu + diag(tau) %*% L.omega %*% eps[i,]
#   }
#   if (p == 0) {
#     for (i in 1:n_i) {
#       ard[i,] = rpois(n_k, nk_prev * d[i] * exp(bias[i,]))
#     }
#   } else{
#     x_cov = matrix(runif(n_i * p), nrow = n_i, ncol = p)
#     x_beta = matrix(0, nrow = p, ncol = n_k)
#     x_beta[sample(p, p_nonzero), ] = runif(p_nonzero * n_k, -2, 2)
#     for (k in 1:n_k) {
#       ard[, k] = rpois(n_i, nk_prev[k] * d * exp(bias[, k] + x_cov %*% x_beta[, k]))
#     }
#   }
#   if (p == 0) {
#     return(
#       list(
#         ard = ard,
#         prev = nk_prev,
#         size = nk_size,
#         omega = Omega,
#         d = d,
#         eta = eta,
#         n_i = n_i,
#         n_k = n_k,
#         seed = seed
#       )
#     )
#   } else{
#     return(
#       list(
#         ard = ard,
#         prev = nk_prev,
#         size = nk_size,
#         omega = Omega,
#         d = d,
#         eta = eta,
#         x_cov = x_cov,
#         x_beta = x_beta,
#         n_i = n_i,
#         n_k = n_k,
#         seed = seed
#       )
#     )
#   }
# }
# 
# 
# # 4. Degree-Correlated
# make_Deg_Correlated <- function(n_i = 500,
#                                 n_k = 20,
#                                 N = 1000000,
#                                 eta = 2,
#                                 p = 9,
#                                 p_nonzero = 3,
#                                 seed = NULL) {
#   if (!is.null(seed)) { # Only set seed if provided
#     set.seed(seed)
#   }
#   nk_prev = runif(n_k, 5.01, 5.15)
#   nk_size = round(nk_prev * N)
#   # d = round(1 + rlnorm(n_i, 5, 0.4))
#   tau.N = runif(n_k, min = 0.5, max = 1.5)
#   mu = c(0, log(1 / sqrt(1 + tau.N ^ 2)))
#   tau = c(0.4, sqrt(log(1 + tau.N ^ 2)))
#   Omega = trialr::rlkjcorr(1, n_k + 1, eta = eta)
#   L.omega = t(chol(Omega))
#   eps = matrix(rnorm(n_i * (n_k + 1)), nrow = n_i, ncol = n_k + 1)
#   bias = matrix(NA, nrow = n_i, ncol = n_k + 1)
#   ard = matrix(NA, nrow = n_i, ncol = n_k)
#   for (i in 1:n_i) {
#     bias[i,] = mu + diag(tau) %*% L.omega %*% eps[i,]
#   }
#   if (p == 0) {
#     for (i in 1:n_i) {
#       ard[i,] = rpois(n_k, nk_prev * exp(bias[i, 1] + bias[i,-1]))
#     }
#   } else{
#     x_cov = matrix(runif(n_i * p), nrow = n_i, ncol = p)
#     x_beta = matrix(0, nrow = p, ncol = n_k)
#     x_beta[sample(p, p_nonzero), ] = runif(p_nonzero * n_k, -2, 2)
#     for (k in 1:n_k) {
#       ard[, k] = rpois(n_i, nk_prev[k] * exp(bias[, 1] + bias[, k + 1] + x_cov %*% x_beta[, k]))
#     }
#   }
#   if (p == 0) {
#     return(
#       list(
#         ard = ard,
#         prev = nk_prev,
#         size = nk_size,
#         omega = Omega,
#         d = d,
#         eta = eta,
#         n_i = n_i,
#         n_k = n_k,
#         seed = seed
#       )
#     )
#   } else{
#     return(
#       list(
#         ard = ard,
#         prev = nk_prev,
#         size = nk_size,
#         omega = Omega,
#         d = d,
#         eta = eta,
#         x_cov = x_cov,
#         x_beta = x_beta,
#         n_i = n_i,
#         n_k = n_k,
#         seed = seed
#       )
#     )
#   }
# }

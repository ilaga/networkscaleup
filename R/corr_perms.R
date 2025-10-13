degree_corr_perm <- function(model_fit,
                             ard,
                             b = 1000,
                             distribution = c("poisson", "nbinomial"),
                             plot = TRUE,
                             attr = c("range", "min", "max", "mean", "median")) {
  # Allow multiple attr values
  attr <- match.arg(attr, several.ok = TRUE)
  
  
  distribution <- match.arg(distribution)
  
  ## Obtain residuals
  if (distribution == "poisson") {
    pois_lambda_est <- model_fit$summary(variables = "mu")$estimate
    alphas <- model_fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(data.frame(value = c(ard)),
                                   family = "poisson",
                                   fit = pois_lambda_est)
    resids <- matrix(resid_vec,
                     nrow = nrow(ard),
                     ncol = ncol(ard))
  } else if (distribution == "nbinomial") {
    nb_prob_est <- model_fit$summary(variables = "inv_omegas")$estimate
    nb_size_est <- model_fit$summary(variables = "par1")$estimate
    alphas <- model_fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(
      data.frame(value = c(ard)),
      family = "negbin",
      size = nb_size_est,
      prob = rep(nb_prob_est, each = n_i)
    )
    resids <- matrix(resid_vec,
                     nrow = nrow(ard),
                     ncol = ncol(ard))
  } else {
    stop("Invalid distribution argument. Must be one of poisson or nbinomial.")
  }
  
  est_df <-
    cbind(alphas, resids)         # Combine degree estimates and residuals
  est_corr <- cor(est_df)                 # Empirical correlation
  
  # Helper: compute statistic from a vector
  compute_stat <- function(x, fun_name) {
    fun <- switch(
      fun_name,
      range  = function(y)
        diff(range(y)),
      min    = min,
      max    = max,
      mean   = mean,
      median = median,
      stop("Unknown attr called.")
    )
    fun(x)
  }
  
  results <- list()
  
  for (a in attr) {
    # Observed statistic
    est_stat <- compute_stat(est_corr[1,-1], a)
    
    # Permutation test
    perm_stat <- numeric(b)
    for (i in 1:b) {
      est_df[, 1] <- sample(est_df[, 1])   # permute degree estimates
      tmp_corr <- cor(est_df)
      perm_stat[i] <- compute_stat(tmp_corr[1,-1], a)
    }
    
    results[[a]] <-
      list(test_stat = est_stat, bootstrap_stat = perm_stat)
  }
  
  # Plotting
  if (plot) {
    n_attr <- length(attr)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    if (n_attr > 1) {
      ncol <- ceiling(sqrt(n_attr))
      nrow <- ceiling(n_attr / ncol)
      par(mfrow = c(nrow, ncol))
    }
    
    for (a in attr) {
      est_stat <- results[[a]]$test_stat
      perm_stat <- results[[a]]$bootstrap_stat
      
      hist(
        perm_stat,
        xlim = range(c(perm_stat, est_stat)),
        xlab = paste("Correlation", a),
        main = paste(
          "Histogram of Correlation",
          if (a == "range")
            "Ranges"
          else
            paste0(a, "s")
        )
      )
      abline(v = est_stat,
             col = "red",
             lwd = 2)
    }
  }
  
  return(results)
}



# Multivariate residual correlation test
## Doesn't really work
multivariate_resid_test <- function(model_fit,
                                    ard,
                                    b = 1000,
                                    plot = TRUE,
                                    distribution = c("poisson", "nbinomial")) {
  
  distribution <- match.arg(distribution)
  
  ## Obtain residuals
  if (distribution == "poisson") {
    pois_lambda_est <- model_fit$summary(variables = "mu")$estimate
    alphas <- model_fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(data.frame(value = c(ard)),
                                   family = "poisson",
                                   fit = pois_lambda_est)
    resids <- matrix(resid_vec,
                     nrow = nrow(ard),
                     ncol = ncol(ard))
  } else if (distribution == "nbinomial") {
    nb_prob_est <- model_fit$summary(variables = "inv_omegas")$estimate
    nb_size_est <- model_fit$summary(variables = "par1")$estimate
    alphas <- model_fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(
      data.frame(value = c(ard)),
      family = "negbin",
      size = nb_size_est,
      prob = rep(nb_prob_est, each = n_i)
    )
    resids <- matrix(resid_vec,
                     nrow = nrow(ard),
                     ncol = ncol(ard))
  } else {
    stop("Invalid distribution argument. Must be one of poisson or nbinomial.")
  }
  
  est_df <- cbind(alphas, resids)
  est_corr <- cor(est_df)
  
  # Off-diagonal sum of squares as test statistic
  obs_stat <- sum((est_corr[1,-1]) ^ 2)
  
  perm_stat <- numeric(b)
  for (i in 1:b) {
    tmp <- est_df
    tmp[, 1] <- sample(tmp[, 1])       # permute degree estimates
    tmp_corr <- cor(tmp)
    perm_stat[i] <- sum((tmp_corr[1,-1]) ^ 2)
  }
  
  p_val <-
    mean(abs(perm_stat - mean(perm_stat)) >= abs(obs_stat - mean(perm_stat)))
  
  if (plot) {
    hist(perm_stat, xlab = "Sum of squared correlations", main = "Multivariate Residual Test")
    abline(v = obs_stat,
           col = "red",
           lwd = 2)
  }
  
  return(list(
    test_stat = obs_stat,
    bootstrap_stat = perm_stat,
    p_value = p_val
  ))
}


# PCA-based residual test for alpha correlation
pca_degree_corr_test <- function(model_fit,
                               ard,
                               b = 1000,
                               distribution = c("poisson", "nbinomial")) {
  distribution <- match.arg(distribution)
  
  ## Obtain residuals
  if (distribution == "poisson") {
    pois_lambda_est <- model_fit$summary(variables = "mu")$estimate
    alphas <- model_fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(data.frame(value = c(ard)),
                                   family = "poisson",
                                   fit = pois_lambda_est)
    resids <- matrix(resid_vec,
                     nrow = nrow(ard),
                     ncol = ncol(ard))
  } else if (distribution == "nbinomial") {
    nb_prob_est <- model_fit$summary(variables = "inv_omegas")$estimate
    nb_size_est <- model_fit$summary(variables = "par1")$estimate
    alphas <- model_fit$summary(variables = "alphas")$estimate
    resid_vec <- construct_pearson(
      data.frame(value = c(ard)),
      family = "negbin",
      size = nb_size_est,
      prob = rep(nb_prob_est, each = n_i)
    )
    resids <- matrix(resid_vec,
                     nrow = nrow(ard),
                     ncol = ncol(ard))
  } else {
    stop("Invalid distribution argument. Must be one of poisson or nbinomial.")
  }
  
  
  est_df <- cbind(alphas, resids)
  
  # PCA on observed data
  obs_pca <- prcomp(est_df, center = TRUE, scale. = TRUE)
  obs_var <-
    obs_pca$sdev[1] ^ 2 / sum(obs_pca$sdev ^ 2)   # fraction of variance explained by first PC
  obs_var_all <-
    obs_pca$sdev ^ 2 / sum(obs_pca$sdev ^ 2)  # all PCs (scree)
  
  # Permutation test
  perm_var <- numeric(b)              # first PC variance
  perm_var_all <-
    matrix(NA, nrow = b, ncol = ncol(est_df)) # all PCs for scree
  
  for (i in 1:b) {
    tmp = est_df
    tmp[, 1] <- sample(est_df[, 1])   # permute alphas
    tmp_pca <- prcomp(tmp, center = TRUE, scale. = TRUE)
    
    perm_var[i] <- tmp_pca$sdev[1] ^ 2 / sum(tmp_pca$sdev ^ 2)
    perm_var_all[i,] <- tmp_pca$sdev ^ 2 / sum(tmp_pca$sdev ^ 2)
  }
  
  # Two-sided p-value
  p_val <-
    mean(abs(perm_var - mean(perm_var)) >= abs(obs_var - mean(perm_var)))
  
  hist_plot = ggplot(data.frame(perm_var), aes(x = perm_var)) +
    geom_histogram(binwidth = diff(range(perm_var)) / 30, fill = "gray80", color = "black") +
    geom_vline(xintercept = obs_var, color = "red", linewidth = 1) +
    labs(
      x = "Variance explained by first PC",
      y = "Count",
      title = "PCA Residual Test"
    ) +
    theme_minimal()
    

  perm_df <- as.data.frame(t(perm_var_all)) |>
    mutate(PC = 1:n()) |>
    pivot_longer(
      cols = -PC,
      names_to = "Permutation",
      values_to = "Variance"
    )
  
  obs_df <- data.frame(
    PC = 1:length(obs_var_all),
    Variance = obs_var_all
  )
  
  ymax <- max(c(obs_var_all, perm_var_all), na.rm = TRUE)
  
  scree_plot = ggplot() +
    geom_line(
      data = perm_df,
      aes(x = PC, y = Variance, group = Permutation),
      color = rgb(0, 0, 0, 0.2)
    ) +
    geom_line(
      data = obs_df,
      aes(x = PC, y = Variance),
      color = "red",
      linewidth = 1.2
    ) +
    geom_point(
      data = obs_df,
      aes(x = PC, y = Variance),
      color = "red"
    ) +
    labs(
      x = "PC index",
      y = "Proportion variance explained",
      title = "Scree plots (obs vs permuted)"
    ) +
    coord_cartesian(ylim = c(0, ymax)) +
    theme_minimal() +
    theme(legend.position = "top") +
    guides(color = "none")
  
  return(
    list(
      test_stat = obs_var,
      bootstrap_stat = perm_var,
      p_value = p_val,
      obs_scree = obs_var_all,
      perm_scree = perm_var_all,
      hist_plot = hist_plot,
      scree_plot = scree_plot
    )
  )
}










# PCA-based residual test for residual correlation (permute each column independently)
pca_group_corr_test <- function(model_fit,
                           ard,
                           b = 1000,
                           distribution = c("poisson", "nbinomial")) {
  distribution <- match.arg(distribution)
  
  ## Obtain residuals
  if (distribution == "poisson") {
    pois_lambda_est <- model_fit$summary(variables = "mu")$estimate
    resid_vec <- construct_pearson(data.frame(value = c(ard)),
                                   family = "poisson",
                                   fit = pois_lambda_est)
    resids <- matrix(resid_vec,
                     nrow = nrow(ard),
                     ncol = ncol(ard))
  } else if (distribution == "nbinomial") {
    nb_prob_est <- model_fit$summary(variables = "inv_omegas")$estimate
    nb_size_est <- model_fit$summary(variables = "par1")$estimate
    resid_vec <- construct_pearson(
      data.frame(value = c(ard)),
      family = "negbin",
      size = nb_size_est,
      prob = rep(nb_prob_est, each = n_i)
    )
    resids <- matrix(resid_vec,
                     nrow = nrow(ard),
                     ncol = ncol(ard))
  } else {
    stop("Invalid distribution argument. Must be one of poisson or nbinomial.")
  }
  
  est_df <- resids
  
  # PCA on observed data
  obs_pca <- prcomp(est_df, center = TRUE, scale. = TRUE)
  obs_var <-
    obs_pca$sdev[1] ^ 2 / sum(obs_pca$sdev ^ 2)        # first PC variance
  obs_var_all <-
    obs_pca$sdev ^ 2 / sum(obs_pca$sdev ^ 2)       # all PCs (scree)
  
  # Permutation test
  n_cols <- ncol(est_df)
  perm_var <- numeric(b)
  perm_var_all <- matrix(NA, nrow = b, ncol = n_cols)
  
  for (i in 1:b) {
    # Permute each column independently
    perm_df <- apply(est_df, 2, sample)
    
    tmp_pca <- prcomp(perm_df, center = TRUE, scale. = TRUE)
    perm_var[i] <- tmp_pca$sdev[1] ^ 2 / sum(tmp_pca$sdev ^ 2)
    perm_var_all[i,] <- tmp_pca$sdev ^ 2 / sum(tmp_pca$sdev ^ 2)
  }
  
  # Two-sided p-value
  p_val <-
    mean(abs(perm_var - mean(perm_var)) >= abs(obs_var - mean(perm_var)))
  
  hist_plot = ggplot(data.frame(perm_var), aes(x = perm_var)) +
    geom_histogram(binwidth = diff(range(perm_var)) / 30, fill = "gray80", color = "black") +
    geom_vline(xintercept = obs_var, color = "red", linewidth = 1) +
    labs(
      x = "Variance explained by first PC",
      y = "Count",
      title = "PCA Residual Test"
    ) +
    theme_minimal()
  
  
  perm_df <- as.data.frame(t(perm_var_all)) |>
    mutate(PC = 1:n()) |>
    pivot_longer(
      cols = -PC,
      names_to = "Permutation",
      values_to = "Variance"
    )
  
  obs_df <- data.frame(
    PC = 1:length(obs_var_all),
    Variance = obs_var_all
  )
  
  ymax <- max(c(obs_var_all, perm_var_all), na.rm = TRUE)
  
  scree_plot = ggplot() +
    geom_line(
      data = perm_df,
      aes(x = PC, y = Variance, group = Permutation),
      color = rgb(0, 0, 0, 0.2)
    ) +
    geom_line(
      data = obs_df,
      aes(x = PC, y = Variance),
      color = "red",
      linewidth = 1.2
    ) +
    geom_point(
      data = obs_df,
      aes(x = PC, y = Variance),
      color = "red"
    ) +
    labs(
      x = "PC index",
      y = "Proportion variance explained",
      title = "Scree plots (obs vs permuted)"
    ) +
    coord_cartesian(ylim = c(0, ymax)) +
    theme_minimal() +
    theme(legend.position = "top") +
    guides(color = "none")
  
  return(
    list(
      test_stat = obs_var,
      bootstrap_stat = perm_var,
      p_value = p_val,
      obs_scree = obs_var_all,
      perm_scree = perm_var_all,
      hist_plot = hist_plot,
      scree_plot = scree_plot
    )
  )
}


# Parametric bootstrap overdispersion test for an ARD matrix
dispersion_test_pb_matrix <- function(ard, mu_hat, B = 1000) {
  # ard: n x K observed counts
  # mu_hat: n x K fitted Poisson means (same dimensions as ard)
  
  stopifnot(all(dim(ard) == dim(mu_hat)))
  n <- nrow(ard)
  K <- ncol(ard)
  
  # helper: Pearson chi-square statistic
  pearson_Q <- function(y, mu)
    sum((y - mu) ^ 2 / mu)
  
  # --- observed stats ---
  Q_obs <- sapply(1:K, function(k)
    pearson_Q(ard[, k], mu_hat[, k]))
  
  # omnibus (sum across groups)
  Q_obs_total <- sum(Q_obs)
  
  # --- bootstrap reference ---
  Q_boot <- matrix(NA, nrow = B, ncol = K)
  Q_boot_total <- numeric(B)
  
  for (b in 1:B) {
    # simulate new data under Poisson null
    Yb <-
      matrix(rpois(n * K, lambda = c(mu_hat)), nrow = n, ncol = K)
    Q_boot[b,] <-
      sapply(1:K, function(k)
        pearson_Q(Yb[, k], mu_hat[, k]))
    Q_boot_total[b] <- sum(Q_boot[b,])
  }
  
  # --- p-values (one-sided for overdispersion) ---
  pvals <- colMeans(t(Q_boot) >= Q_obs)
  pval_total <- mean(Q_boot_total >= Q_obs_total)
  
  return(
    list(
      stat = Q_obs,
      stat_total = Q_obs_total,
      ref = Q_boot,
      ref_total = Q_boot_total,
      pvals = pvals,
      pval_total = pval_total
    )
  )
}




# Stratified-permutation overdispersion test for ARD matrix (Option B, full matrix)
dispersion_test_stratperm_matrix <-
  function(Y,
           MU,
           B = 2000,
           n_bins = 10,
           seed = NULL) {
    # Y: n x K observed counts
    # MU: n x K fitted Poisson means (same dims as Y)
    # B: number of permutations
    # n_bins: number of quantile bins for MU within each column
    # seed: optional RNG seed for reproducibility
    
    if (!is.null(seed))
      set.seed(seed)
    stopifnot(all(dim(Y) == dim(MU)))
    
    n <- nrow(Y)
    K <- ncol(Y)
    
    # Pearson dispersion statistic for one column
    pearson_Q_col <- function(y, mu)
      sum((y - mu) ^ 2 / mu)
    
    # Observed statistics
    Q_obs <- numeric(K)
    for (k in seq_len(K))
      Q_obs[k] <- pearson_Q_col(Y[, k], MU[, k])
    Q_obs_total <- sum(Q_obs)
    
    # Precompute bin indices for each column
    bins_list <- vector("list", K)
    for (k in seq_len(K)) {
      mu_k <- MU[, k]
      # if all MU identical, single bin
      if (all(mu_k == mu_k[1])) {
        bins_list[[k]] <- rep(1L, n)
      } else {
        probs <- seq(0, 1, length.out = n_bins + 1)
        breaks <- unique(quantile(mu_k, probs = probs, na.rm = TRUE))
        # ensure at least two breaks
        if (length(breaks) == 1) {
          bins_list[[k]] <- rep(1L, n)
        } else {
          bins_list[[k]] <-
            as.integer(cut(mu_k, breaks = breaks, include.lowest = TRUE))
        }
      }
    }
    
    # Prepare storage
    Q_boot <- matrix(NA_real_, nrow = B, ncol = K)
    Q_boot_total <- numeric(B)
    
    # Permutation loop
    for (b in seq_len(B)) {
      Yb <- matrix(NA_integer_, nrow = n, ncol = K)
      
      for (k in seq_len(K)) {
        idx_grp <- split(seq_len(n), bins_list[[k]])
        y_col <- Y[, k]
        y_perm_col <- integer(n)
        # permute within each bin
        for (gid in names(idx_grp)) {
          ids <- idx_grp[[gid]]
          if (length(ids) == 1) {
            y_perm_col[ids] <- y_col[ids]
          } else {
            y_perm_col[ids] <- sample(y_col[ids], length(ids), replace = FALSE)
          }
        }
        Yb[, k] <- y_perm_col
      }
      
      # compute Q per column (compare to MU of same column)
      for (k in seq_len(K))
        Q_boot[b, k] <- pearson_Q_col(Yb[, k], MU[, k])
      Q_boot_total[b] <- sum(Q_boot[b,])
    }
    
    # p-values (right tail = test for overdispersion)
    pvals <- numeric(K)
    for (k in seq_len(K))
      pvals[k] <- mean(Q_boot[, k] >= Q_obs[k])
    pval_total <- mean(Q_boot_total >= Q_obs_total)
    
    # return
    list(
      stat = Q_obs,
      stat_total = Q_obs_total,
      ref = Q_boot,
      ref_total = Q_boot_total,
      pvals = pvals,
      pval_total = pval_total,
      bins_list = bins_list
    )
  }

#############################################################
## Simulation 1:
## Model: Poisson likelihood, 6 covariates, 2 local, 1 global, no correlation
## Output: Visualizations for model choice
##
#############################################################


library(ggplot2)
library(tidyverse)

source("Make_ARD.R")
source("fit_stan_optim.R")
source("hang_rootogram.R")
source("corr_plot.R")
source("corr_perms.R")
source("cov_plots.R")


set.seed(2)
sim_dat = make_ARD(
  p = 6,
  p_local_nonzero = 2,
  p_global_nonzero = 1,
  distribution = "poisson"
)


colnames(sim_dat$ard) =
  paste0("g", 1:ncol(sim_dat$ard))

n_i = nrow(sim_dat$ard)
n_k = ncol(sim_dat$ard)

sim_dat$x_beta_global
sim_dat$x_beta_local


# Check covariate structure first -----------------------------------------

## Get basic model fit
pois_fit_list = fit_stan_optim(sim_dat$ard, family = "poisson", init = 0)

cov_p_list = cov_plots(
  pois_fit_list,
  sim_dat$ard,
  x_cov = sim_dat$x_cov,
  family = "poisson",
  se = T
)

local_cov_plot = cov_p_list$Group_plot
global_cov_plot = cov_p_list$Respondent_plot

local_cov_plot # Suggests x3 and x6
global_cov_plot # Suggests x5 and x6, a little x4

# Check correlation structure ---------------------------------------------

## Refit model with chosen covariate structure


pois_fit_corr_list = fit_stan_optim(
  sim_dat$ard,
  x_cov_global = as.matrix(sim_dat$x_cov[, c(5)]),
  x_cov_local = as.matrix(sim_dat$x_cov[, c(3, 6)]),
  family = "poisson",
  init = 0
)

## Run correlation tests and visualizations

degree_corr_res = pca_degree_corr_test(pois_fit_corr_list, sim_dat$ard)
group_corr_res = pca_group_corr_test(pois_fit_corr_list, sim_dat$ard)

group_corr_hist = group_corr_res$hist_plot + ggtitle("Group Correlation PCA Residual Test")
group_corr_scree = group_corr_res$scree_plot + ggtitle("Group Correlation PCA Scree Plot")
degree_corr_hist = degree_corr_res$hist_plot + ggtitle("Degree Correlation PCA Residual Test")
degree_corr_scree = degree_corr_res$scree_plot + ggtitle("Degree Correlation PCA Scree Plot")




# Check distribution ------------------------------------------------------


## No changes, so look at rootogram

## Poisson rootogram
pois_root = hang_rootogram_ard(sim_dat$ard,
                   pois_fit_corr_list) + ggtitle("Poisson Hanging Rootogram")

## Looks great

nb_fit_corr_list = fit_stan_optim(
  sim_dat$ard,
  x_cov_global = as.matrix(sim_dat$x_cov[, c(5)]),
  x_cov_local = as.matrix(sim_dat$x_cov[, c(3, 6)]),
  family = "nbinomial",
  init = 0
)
# Trouble converging because variance \approx mean


## nbinomial rootogram
negbin_root = hang_rootogram_ard(sim_dat$ard, nb_fit_corr_list, family = "nbinomial") + ggtitle("nbinomial Hanging Rootogram")
## Look almost identical



# Save results ------------------------------------------------------------



ggsave("./Figures/sim1_local_cov.pdf", local_cov_plot, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim1_global_cov.pdf", global_cov_plot, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim1_group_corr_hist.pdf", group_corr_hist, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim1_group_corr_scree.pdf", group_corr_scree, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim1_degree_corr_hist.pdf", degree_corr_hist, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim1_degree_corr_scree.pdf", degree_corr_scree, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim1_pois_root.pdf", pois_root, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim1_nbinomial_root.pdf", negbin_root, units = c("in"), width = 12, height = 10)









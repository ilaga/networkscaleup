#############################################################
## Simulation 4:
## Model: Poisson likelihood, 6 covariates, 2 local, 2 global, group correlation
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
  p_global_nonzero = 2,
  group_corr = T,
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

local_cov_plot # Suggests 1 and 5
global_cov_plot # Suggests 4 and 6


# Check correlation structure ---------------------------------------------

## Refit model with chosen covariate structure


pois_fit_corr_list = fit_stan_optim(
  sim_dat$ard,
  x_cov_global = as.matrix(sim_dat$x_cov[, c(4, 6)]),
  x_cov_local = as.matrix(sim_dat$x_cov[, c(1, 5)]),
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


## Has correlation, so can't check distribution




# Save results ------------------------------------------------------------



ggsave("./Figures/sim4_local_cov.pdf", local_cov_plot, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim4_global_cov.pdf", global_cov_plot, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim4_group_corr_hist.pdf", group_corr_hist, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim4_group_corr_scree.pdf", group_corr_scree, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim4_degree_corr_hist.pdf", degree_corr_hist, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim4_degree_corr_scree.pdf", degree_corr_scree, units = c("in"), width = 12, height = 10)







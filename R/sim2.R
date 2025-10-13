#############################################################
## Simulation 2:
## Model: Negative Binomial likelihood, 6 covariates, 0 local, 0 global, group correlation
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
  p_local_nonzero = 0,
  p_global_nonzero = 0,
  group_corr = T,
  distribution = "nbinomial"
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

local_cov_plot # Suggests nothing
global_cov_plot # Suggests nothing


# Check correlation structure ---------------------------------------------

## Covariates unchanged, use old model

## Run correlation tests and visualizations

degree_corr_res = pca_degree_corr_test(pois_fit_list, sim_dat$ard)
group_corr_res = pca_group_corr_test(pois_fit_list, sim_dat$ard)

group_corr_hist = group_corr_res$hist_plot + ggtitle("Group Correlation PCA Residual Test")
group_corr_scree = group_corr_res$scree_plot + ggtitle("Group Correlation PCA Scree Plot")
degree_corr_hist = degree_corr_res$hist_plot + ggtitle("Degree Correlation PCA Residual Test")
degree_corr_scree = degree_corr_res$scree_plot + ggtitle("Degree Correlation PCA Scree Plot")



# Check distribution ------------------------------------------------------

## Has correlation, do can't check distribution



# Save results ------------------------------------------------------------



ggsave("./Figures/sim2_local_cov.pdf", local_cov_plot, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim2_global_cov.pdf", global_cov_plot, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim2_group_corr_hist.pdf", group_corr_hist, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim2_group_corr_scree.pdf", group_corr_scree, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim2_degree_corr_hist.pdf", degree_corr_hist, units = c("in"), width = 12, height = 10)
ggsave("./Figures/sim2_degree_corr_scree.pdf", degree_corr_scree, units = c("in"), width = 12, height = 10)









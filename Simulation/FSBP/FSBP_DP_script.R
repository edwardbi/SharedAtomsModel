rm(list=ls())

##################################################################
# Command line arguments
##################################################################

args = commandArgs(TRUE)

#if(length(args) != 5){
#  stop("Usage: SAM1_m1_script.R <seed_n> <work_dir> <input_dat> <sim_data_num> <out_dir>")
#}

seed_n = 12345 #as.numeric(args[1])
work_dir = "/Users/dehuabi/Desktop/PAM_Extra_Sim/FSBP_Sim/" #args[2]
input_dat = "./data/data_300.RData" #args[3]
#sim_data_num = as.numeric(args[1])
sim_data_num = 1
out_dir = "./results/300/" #args[5]

setwd(work_dir)
set.seed(seed_n)
print(sim_data_num)

##################################################################
# Pcakages and Sources
##################################################################

source("./code/FSBP_slice.R")
library(dirichletprocess)

# Read in simulation data
sim_data_set <- readRDS(input_dat)
sim_data <- sim_data_set[[sim_data_num]]
d_vec <- sim_data$d_vec
c_vec <- sim_data$c_vec

# Hyperparameters of FSBP
mu0 <- 0
nu0 <- 0.1
a0 <- 3
b0 <- 1
a1 <- 0.5
b1 <- 0.5
a2 <- 3
b2 <- 3
N_iter <- 5000
burn_in <- 5000
n <- 450

# FSBP model

fsbp_fit <- MCMC_FSBP(d_vec, mu0, nu0, a0, b0, a1, b1, a2, b2, N_iter, burn_in, iter_num = 1000)
fsbp_metrics_result <- comp_metrics(fsbp_fit, c_vec, burn_in, N_iter, n)

fsbp_post_mcmc_dir <- paste0(out_dir, "FSBP_MCMC_", seed_n, "_", sim_data_num, ".RData")
saveRDS(fsbp_fit, file = fsbp_post_mcmc_dir)
fsbp_post_optc_dir <- paste0(out_dir, "FSBP_optcluster_", seed_n, "_", sim_data_num, ".RData")
saveRDS(fsbp_metrics_result, file = fsbp_post_optc_dir)

# DP model

DP_fit <- DirichletProcessGaussian(d_vec, g0Priors = c(mu0, nu0, a0, b0), alphaPriors = c(a2,b2))
DP_fit <- Fit(DP_fit, 10000, progressBar = T)

rand_idx_DP <- adjustedRandIndex(DP_fit$clusterLabels, c_vec)

all_mat_DP <- matrix(rep(NA, (n*n)), ncol = n)
for(i in 1:n){
  for(i1 in 1:n){
    if(DP_fit$clusterLabels[i] == DP_fit$clusterLabels[i1]){
      all_mat_DP[i,i1] <- 1
    }else{
      all_mat_DP[i,i1] <- 0
    }
  }
}
all_mat_true <- matrix(rep(NA, (n*n)), ncol = n)
for(i in 1:n){
  for(i1 in 1:n){
    if(c_vec[i] == c_vec[i1]){
      all_mat_true[i,i1] <- 1
    }else{
      all_mat_true[i,i1] <- 0
    }
  }
}

ndf_val_DP <- nfd(all_mat_DP, all_mat_true)
dp_metrics_result <- list(rand_idx = rand_idx_DP, nfd = ndf_val_DP, clusters = DP_fit$numberClusters)

dp_post_mcmc_dir <- paste0(out_dir, "DP_MCMC_", seed_n, "_", sim_data_num, ".RData")
saveRDS(DP_fit, file = dp_post_mcmc_dir)
dp_post_optc_dir <- paste0(out_dir, "DP_optcluster_", seed_n, "_", sim_data_num, ".RData")
saveRDS(dp_metrics_result, file = dp_post_optc_dir)

# Output results
print("FSBP:")
print(fsbp_metrics_result)
print("DP:")
print(dp_metrics_result)

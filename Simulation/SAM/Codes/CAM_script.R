rm(list=ls())

##################################################################
# Command line arguments
##################################################################

args = commandArgs(TRUE)

#if(length(args) != 5){
#  stop("Usage: CAM_script.R <seed_n> <work_dir> <input_dat> <sim_data_num> <out_dir>")
#}

seed_n = as.numeric(args[1])
work_dir = args[2]
input_dat = args[3]
sim_data_num = as.numeric(args[4])
out_dir = args[5]
all_flag = as.logical(args[6]) #T

setwd(work_dir)
set.seed(seed_n)
print(sim_data_num)

##################################################################
# Pcakages and Sources
##################################################################
library(Rcpp)
library(tidyverse)
library(mcclust)
library(binhf)
sourceCpp("./codes/CAM/newCAM.cpp")
source("./codes/CAM/newCAM.R")
source("./codes/comp_stats.R")
source("./codes/post_mcmc_ls.R")

# Read in simulation data
sim_data_set <- readRDS(input_dat)
sim_data <- sim_data_set[[sim_data_num]]
d_mat <- sim_data$d_mat
c_mat <- sim_data$c_mat
num_p_j <- sim_data$num_p_j

d_input <- as.vector(t(d_mat))
d_input <- d_input[!is.na(d_input)]
g_input <- rep(NA, length(d_input))#as.vector(t(c_mat))
for(j in 1:length(num_p_j)){
  if(j == 1){
    start_idx <- 1
  }else{
    start_idx <- sum(num_p_j[1:(j-1)])+1
  }
  end_idx <- sum(num_p_j[1:j])
  g_input[start_idx:end_idx] <- rep(j, num_p_j[j])
}

# Hyperparameters of SAM 1
a0 <- 3
b0 <- 1

a1 <- 3
b1 <- 3

a2 <- 3
b2 <- 3

M_iter <- 5000
burnin <- 5000

K <- 10

# CAM and save data
if(all_flag == T){
  CAM_result <- CAM(y_obser = d_input, y_group = g_input, K0 = 3, L0 = K,
                    prior = list(# hyperparameters NIG
                      m0=0, k0=0.1, a0=a0, b0=b0,
                      # hyperparameters alpha and beta
                      a_alpha=a1, b_alpha = b1, a_beta =a2, b_beta = b2),
                    nsim = M_iter, burn_in = burnin, thinning = 1,verbose = 1,
                    fixedAB = F, kappa=0.5, cheap=F, seed=seed_n)
  post_mcmc_dir <- paste0(out_dir, "CAM_MCMC_univ_", seed_n, "_", sim_data_num, ".RData")
  saveRDS(CAM_result, file = post_mcmc_dir)
}else{
  post_mcmc_dir <- paste0(out_dir, "CAM_MCMC_univ_", seed_n, "_", sim_data_num, ".RData")
  CAM_result <- readRDS(post_mcmc_dir)
}

# Rand index and NFD
cls <- CAM_result$Csi_ij
psm <- comp.psm(cls)
if(all_flag == T){
  optmal_cluster <- minVI(psm, cls, method="all", include.greedy = TRUE)
}else{
  optmal_cluster <- minVI(psm, cls)#, method="all", include.greedy = TRUE)
}
summary_result <- comp_randidx_ndf(num_p_j, c_mat, psm, optmal_cluster, all_flag = all_flag)
post_optc_dir <- paste0(out_dir, "CAM_univ_optcluster_", seed_n, "_", sim_data_num, ".RData")
saveRDS(optmal_cluster, file = post_optc_dir)
print(summary_result)
if(all_flag == T){
  clus <- optmal_cluster$cl[1,]
}else{
  clus <- optmal_cluster$cl
}

num_c <- length(unique(clus))
uniq_g1 <- unique(clus[1:num_p_j[1]])
uniq_g2 <- unique(clus[(num_p_j[1]+1):sum(num_p_j)])
#uniq_g3 <- unique(clus[(sum(num_p_j[1:2])+1):sum(num_p_j)])

# Number of clusters:
print("Point estimates:")
print(paste0("Total number of clusters: ", num_c, ", g1 clusters: ", length(uniq_g1), ", g2 clusters: ",length(uniq_g2)))
print(paste0("Number of shared clusters:",length(intersect(uniq_g1,uniq_g2))))
#print(paste0("Number of shared clusters g1g3:",length(intersect(uniq_g1,uniq_g3))))
#print(paste0("Number of shared clusters g2g3:",length(intersect(uniq_g2,uniq_g3))))
#print(paste0("Number of shared clusters all:",length(intersect(intersect(uniq_g1,uniq_g2),uniq_g3))))
print(paste("Number of unique clusters in g1:",length(setdiff(uniq_g1, uniq_g2))))
print(paste("Number of unique clusters in g2:",length(setdiff(uniq_g2, uniq_g1))))
#print(paste("Number of unique clusters in g3:",length(setdiff(uniq_g3, union(uniq_g1,uniq_g2)))))

inc_prob_point_est <- function(c_vec, clus, c_id){
  unique_clus <- unique(clus)
  inc_prob <- rep(NA, length(unique_clus))
  for(i in 1:length(unique_clus)){
    inc_prob[i] <- sum(clus[which(c_vec == c_id)] == i)/length(which(c_vec == c_id))
  }
  return(inc_prob)
}

print("Inclusion probability:")
for(i1 in unique(as.vector(t(c_mat)))){
  print(paste("Clusters in truth:",i1))
  print(inc_prob_point_est(as.vector(t(c_mat)), clus, i1))
}

for(i2 in unique(clus)){
  print(paste("Clusters in est:",i2))
  print(inc_prob_point_est(clus, as.vector(t(c_mat)), i2))
}


# Label switching
post_OMG <- list()
post_mu <- list()
post_z <- list()
for(iter in 1:M_iter){
  post_OMG_temp <- t(CAM_result$OMEGA[[iter]])
  post_Z_j <- sort(CAM_result$Z_j[iter,])
  post_OMG[[iter]] <- post_OMG_temp[post_Z_j,]
  mus <- CAM_result$THETA[[iter]]
  post_mu[[iter]] <- mus[,1]
  post_z[[iter]] <- CAM_result$Csi_ij[iter,]
}
label_sw_result <- process_mcmc(d_mat, num_p_j, post_z, post_OMG, 
                                post_mu, 0, M_iter, max_K = 50, mat_flag = F)
relabel_result <- relabel_mcmc(M_iter, label_sw_result$K, label_sw_result$mcmc_par, label_sw_result$J, 
                               label_sw_result$cov_p, label_sw_result$prem_label)
post_sw_dir <- paste0(out_dir, "CAM_univ_label_sw_", seed_n, "_", sim_data_num, ".RData")
saveRDS(label_sw_result$prem_label, file = post_sw_dir)


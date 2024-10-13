rm(list=ls())

##################################################################
# Command line arguments
##################################################################

args = commandArgs(TRUE)

#if(length(args) != 5){
#  stop("Usage: SAM1_m1_script.R <seed_n> <work_dir> <input_dat> <sim_data_num> <out_dir>")
#}

seed_n = as.numeric(args[1])
work_dir = args[2]
input_dat = args[3]
sim_data_num = as.numeric(args[4])
out_dir = args[5]
all_flag = as.logical(args[6]) #T#F

setwd(work_dir)
set.seed(seed_n)
print(sim_data_num)

##################################################################
# Pcakages and Sources
##################################################################
source("./codes/PAM_slice.R")
source("./codes/comp_stats.R")
source("./codes/post_mcmc_ls.R")

# Read in simulation data
sim_data_set <- readRDS(input_dat)
sim_data <- sim_data_set[[sim_data_num]]
d_mat <- sim_data$d_mat
c_mat <- sim_data$c_mat
num_p_j <- sim_data$num_p_j

# Hyperparameters of SAM 1
mu0 <- 0
nu0 <- 0.1

a <- 0.5
b <- 0.5

a0 <- 3
b0 <- 1

a1 <- 3
b1 <- 3

a2 <- 3
b2 <- 3

M_iter <- 5000
burnin <- 5000
K <- 10

# SAM 1 and save data
if(all_flag == T){
  post_result <- MCMC_SAM1(d_mat, mu0, nu0, a, b, a0, b0, a1, b1, a2, b2, 0.3, 1, 0.5, K, M_iter, burnin, iter_num = 1000)
  post_mcmc_dir <- paste0(out_dir, "PAM_MCMC_univ_", seed_n, "_", sim_data_num, ".RData")
  saveRDS(post_result, file = post_mcmc_dir)
}else{
  post_mcmc_dir <- paste0(out_dir, "PAM_MCMC_univ_", seed_n, "_", sim_data_num, ".RData")
  post_result <- readRDS(post_mcmc_dir)
}

# Rand index and NFD
clust_result <- comp_psm_optclust(num_p_j, burnin, M_iter, post_result$all_z_mat, all_flag = all_flag)
summary_result <- comp_randidx_ndf(num_p_j, c_mat, clust_result$psm, clust_result$optclust, all_flag = all_flag)
post_optc_dir <- paste0(out_dir, "PAM_univ_optcluster_", seed_n, "_", sim_data_num, ".RData")
saveRDS(clust_result$optclust, file = post_optc_dir)
print(summary_result)
optmal_cluster <- clust_result$optclust
if(all_flag == T){
  clus <- optmal_cluster$cl[1,]
}else{
  clus <- optmal_cluster$cl
}

num_c <- length(unique(clus))
uniq_g1 <- unique(clus[1:num_p_j[1]])
uniq_g2 <- unique(clus[(num_p_j[1]+1):sum(num_p_j[1:2])])
#uniq_g3 <- unique(clus[(sum(num_p_j[1:2])+1):sum(num_p_j)])

# Number of clusters:
print("Point estimates:")
print(paste0("Total number of clusters: ", num_c, ", g1 clusters: ", length(uniq_g1), ", g2 clusters: ",length(uniq_g2)))
print(paste0("Number of shared clusters all:",length(intersect(uniq_g1,uniq_g2))))
print(paste("Number of unique clusters in g1:",length(setdiff(uniq_g1, uniq_g2))))
print(paste("Number of unique clusters in g2:",length(setdiff(uniq_g2, uniq_g1))))

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
#label_sw_result <- process_mcmc(d_mat, num_p_j, post_result$all_z_mat, post_result$all_pi_p_mat, 
#                                post_result$all_phi_vec, burnin, M_iter, max_K = 50)
#relabel_result <- relabel_mcmc(M_iter, label_sw_result$K, label_sw_result$mcmc_par, label_sw_result$J, 
#                               label_sw_result$cov_p, label_sw_result$prem_label)
#post_sw_dir <- paste0(out_dir, "PAM_univ_label_sw_", seed_n, "_", sim_data_num, ".RData")
#saveRDS(label_sw_result$prem_label, file = post_sw_dir)

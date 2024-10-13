rm(list=ls())

##################################################################
# Command line arguments
##################################################################

args = commandArgs(TRUE)

#if(length(args) != 5){
#  stop("Usage: HDP_script.R <seed_n> <work_dir> <input_dat> <sim_data_num> <out_dir>")
#}

seed_n = as.numeric(args[1])
work_dir = args[2]
input_dat = args[3]
sim_data_num = as.numeric(args[4])
out_dir = args[5]#args[5]

setwd(work_dir)
set.seed(seed_n)
print(sim_data_num)

##################################################################
# Pcakages and Sources
##################################################################
source("./codes/HDP.R")
source("./codes/comp_stats.R")
source("./codes/post_mcmc_ls.R")

# Read in simulation data
sim_data_set <- readRDS(input_dat)
sim_data <- sim_data_set[[sim_data_num]]
d_mat <- sim_data$d_mat
c_mat <- sim_data$c_mat
num_p_j <- sim_data$num_p_j

M_iter <- 1000
burnin <- 1000

mu0 <- 0
lambda0 <- 0.1
a0 <- 3
b0 <- 1

a1 <- 3
b1 <- 3

a2 <- 3
b2 <- 3

# HDP and save data
post_result <- mcmc_alg(M_iter, burnin, d_mat, 1, 1, c(mu0, lambda0, a0, b0), a1, a2, b1, b2)
post_mcmc_dir <- paste0(out_dir, "HDP_MCMC_univ_", seed_n, "_", sim_data_num, ".RData")
saveRDS(post_result, file = post_mcmc_dir)

# Rand index and NFD
clust_result <- comp_psm_optclust(num_p_j, burnin, M_iter, post_result$z_mat_list)
summary_result <- comp_randidx_ndf(num_p_j, c_mat, clust_result$psm, clust_result$optclust)
post_optc_dir <- paste0(out_dir, "HDP_univ_optcluster_", seed_n, "_", sim_data_num, ".RData")
saveRDS(clust_result$optclust, file = post_optc_dir)
print(summary_result)
optmal_cluster <- clust_result$optclust
clus <- optmal_cluster$cl[1,]

num_c <- length(unique(clus))
uniq_g1 <- unique(clus[1:num_p_j[1]])
uniq_g2 <- unique(clus[(num_p_j[1]+1):sum(num_p_j[1:2])])
uniq_g3 <- unique(clus[(sum(num_p_j[1:2])+1):sum(num_p_j)])

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




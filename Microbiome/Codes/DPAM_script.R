rm(list=ls())

##################################################################
# Command line arguments
##################################################################

args = commandArgs(TRUE)

if(length(args) != 5){
  stop("Usage: SAM1_m1_script.R <seed_n> <work_dir> <input_dat> <sim_data_num> <out_dir>")
}

seed_n = as.numeric(args[1])
work_dir = args[2]

#setwd(work_dir)
setwd("/Users/dehuabi/Desktop/SAM 1 Final/Microbiome/")
#set.seed(seed_n)
set.seed(12345)
##################################################################
# Pcakages and Sources
##################################################################
source("./Codes/DPAM_slice.R")
source("./Codes/comp_stats.R")
library(phyloseq)

dietswap_data <- readRDS("./Data/Icondietswap.RData")
otu_table_data <- dietswap_data@otu_table@.Data
tax_table_data <- dietswap_data@tax_table@.Data
sam_table_data <- dietswap_data@sam_data@.Data

#saveRDS(list(otu = otu_table_data, tax = tax_table_data, samp = sam_table_data), 
#        file = "dietwap.RData")

samples_all <- sam_table_data[[1]]
samples <- unique(samples_all)
j_idx <- rep(NA, length(samples))
count <- 1
for(pat in samples){
  j_idx[count] <- which(samples_all == pat)[1]
  count <- count + 1
}

first_otu_table <- otu_table_data[,j_idx]
non_zero_idxs <- c()
for(i in 1:130){
  row_sum <- sum(first_otu_table[i,])
  if(row_sum != 0){
    non_zero_idxs <- c(non_zero_idxs, i)
  }
}

CAM_otu_table <- first_otu_table[non_zero_idxs, ]
tax_names <- tax_table_data[,1]
tax_names <- tax_names[non_zero_idxs]
col_avg <- rep(NA, 38)
for(j in 1:38){
  col_avg[j] <- mean(CAM_otu_table[,j])
}

pat_nation <- sam_table_data[[3]]
pat_nation <- pat_nation[j_idx]
AF_samples <- which(pat_nation == "AFR")
AM_samples <- which(pat_nation == "AAM")

select_pat_mat <- CAM_otu_table[,c(5,22,13,14)]

non_zero_idxs_sel <- c()
for(i in 1:119){
  row_sum <- sum(select_pat_mat[i,])
  if(row_sum != 0){
    non_zero_idxs_sel <- c(non_zero_idxs_sel, i)
  }
}
select_pat_mat <- select_pat_mat[non_zero_idxs_sel,]
tax_names_sel <- tax_names[non_zero_idxs_sel]
d_mat <- matrix(rep(NA, 4*length(select_pat_mat[,1])), nrow = 4)
for(j in 1:length(select_pat_mat[1,])){
  d_mat[j,] <- select_pat_mat[,j]
}

saveRDS(tax_names_sel, file = "tax_names_1sim.RData")

id_list <- c(5, 22, 13, 14)
belonging <- c("AF", "AF", "AA", "AA")
par(mfrow = c(2,2))
for(j in 1:4){
  hist(d_mat[j,], freq = F, nclass = 50, main = paste0("Subject ",id_list[j]," of ",belonging[j]),
       xlab = "OTU counts")
}

# Running MCMC
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

M_iter <- 10000
burnin <- 10000
K <- 10

post_result <- MCMC_SAM1(d_mat, col_avg, mu0, nu0, a, b, a0, b0, a1, b1, a2, b2, 0.3, 1, 0.5, K, M_iter, burnin)

#saveRDS(post_result, file = "DSAM_microbiome_1sim.RData")

# optimal clusters
opt_clust <- comp_psm_optclust(rep(109,4), burnin, M_iter, post_result$all_z_mat)
#saveRDS(opt_clust, file = "DSAM_microbiome_1sim_opt_clusters.RData")

optimal_cluster <- opt_clust$optclust$cl[1,]
print(length(unique(optimal_cluster)))

for(j_idx in 1:4){
  opt_c_j <- optimal_cluster[((j_idx-1)*109+1):(j_idx*109)]
  print(sort(unique(opt_c_j)))
}


K <- 50
all_pi_p_1 <- matrix(rep(NA, M_iter*K), ncol = K)
all_pi_p_2 <- matrix(rep(NA, M_iter*K), ncol = K)
all_pi_p_3 <- matrix(rep(NA, M_iter*K), ncol = K)
all_pi_p_4 <- matrix(rep(NA, M_iter*K), ncol = K)
all_phi <- matrix(rep(NA, M_iter*K), ncol = K)
all_sig <- matrix(rep(NA, M_iter*K), ncol = K)
max_length <- 0

for(iter in (burnin+1):(burnin+M_iter)){
  #print(iter)
  pi_p_mat <- post_result$all_pi_p_mat[[iter]]
  all_pi_p_1[iter-burnin,1:length(pi_p_mat[1,])] <- pi_p_mat[1,]
  all_pi_p_2[iter-burnin,1:length(pi_p_mat[1,])] <- pi_p_mat[2,]
  all_pi_p_3[iter-burnin,1:length(pi_p_mat[1,])] <- pi_p_mat[3,]
  all_pi_p_4[iter-burnin,1:length(pi_p_mat[1,])] <- pi_p_mat[4,]
  if(length(pi_p_mat[1,]) > max_length){
    max_length <- length(pi_p_mat[1,])
  }
  all_phi[iter-burnin,1:length(pi_p_mat[1,])] <- post_result$all_phi_vec[[iter]]
  all_sig[iter-burnin,1:length(pi_p_mat[1,])] <- post_result$all_sig_vec[[iter]]
}

library(label.switching)

z_mat_mcmc <- matrix(rep(NA, M_iter*4*109), nrow = M_iter)
mcmc_par <- array(dim = c(M_iter, max_length, 6))

for(iter in (burnin+1):(burnin+M_iter)){
  z_mat_mcmc_val <- as.vector(t(post_result$all_z_mat[[iter]]))
  z_mat_mcmc_val <- z_mat_mcmc_val[!is.na(z_mat_mcmc_val)]
  z_mat_mcmc[iter-burnin,] <- z_mat_mcmc_val
  pi_p_mat <- post_result$all_pi_p_mat[[iter]]
  mcmc_par[iter-burnin,1:length(pi_p_mat[1,]),1] <- pi_p_mat[1,]
  mcmc_par[iter-burnin,1:length(pi_p_mat[1,]),2] <- pi_p_mat[2,]
  mcmc_par[iter-burnin,1:length(pi_p_mat[1,]),3] <- pi_p_mat[3,]
  mcmc_par[iter-burnin,1:length(pi_p_mat[1,]),4] <- pi_p_mat[4,]
  mcmc_par[iter-burnin,1:length(pi_p_mat[1,]),5] <- post_result$all_phi_vec[[iter]]
  mcmc_par[iter-burnin,1:length(pi_p_mat[1,]),6] <- post_result$all_sig_vec[[iter]]
}

d_mat_vec <- as.vector(t(d_mat))
d_mat_vec <- d_mat_vec[!is.na(d_mat_vec)]

ls <- label.switching(method=c("ECR"),#, "DATA-BASED"),#, "ECR-ITERATIVE-1"), 
                      zpivot=z_mat_mcmc[1,],
                      z = z_mat_mcmc,K = max_length, data = d_mat_vec, mcmc = mcmc_par)#as.vector(d_mat))
prem_ls <- ls$permutations[["ECR"]]

#saveRDS(prem_ls, file = "DSAM_microbiome_1sim_label_switch.RData")

# ConstWeight Function
const_weight <- function(w_p){
  
  m_wp <- 1 - w_p
  m_wp2 <- shift(cumprod(m_wp), 1)
  m_wp2[1] <- 1
  
  w <- w_p*m_wp2
  return(w)
  
}

de_const_weight <- function(w){
  w_p <- rep(NA, length(w))
  o_m_w_p <- w_p
  for(k in 1:length(w_p)){
    if(k == 1){
      w_p[1] <- w[1]
      o_m_w_p[1] <- 1 - w_p[1]
    }else{
      w_p[k] <- w[k]/(prod(o_m_w_p[1:(k-1)]))
      o_m_w_p[k] <- 1 - w_p[k]
    }
  }
  return(w_p)
}

all_phi_swap <- all_phi[,1:max_length]
all_sig_swap <- all_sig[,1:max_length]
all_pi_p_1_swap <- all_pi_p_1[,1:max_length]
all_pi_p_2_swap <- all_pi_p_2[,1:max_length]
all_pi_p_3_swap <- all_pi_p_3[,1:max_length]
all_pi_p_4_swap <- all_pi_p_4[,1:max_length]
all_pi_1_swap <- all_pi_p_1[,1:max_length]
all_pi_2_swap <- all_pi_p_2[,1:max_length]
all_pi_3_swap <- all_pi_p_3[,1:max_length]
all_pi_4_swap <- all_pi_p_4[,1:max_length]
for(i in 1:M_iter){
  all_phi_swap[i,] <- all_phi_swap[i,prem_ls[i,]]
  all_sig_swap[i,] <- all_sig_swap[i,prem_ls[i,]]
  
  prob_all_1 <- const_weight(all_pi_p_1_swap[i,])
  prob_all_2 <- const_weight(all_pi_p_2_swap[i,])
  prob_all_3 <- const_weight(all_pi_p_3_swap[i,])
  prob_all_4 <- const_weight(all_pi_p_4_swap[i,])
  
  #print(prob_all_1)
  
  prob_1_swap <- prob_all_1[prem_ls[i,]]
  prob_2_swap <- prob_all_2[prem_ls[i,]]
  prob_3_swap <- prob_all_3[prem_ls[i,]]
  prob_4_swap <- prob_all_4[prem_ls[i,]]
  
  all_pi_1_swap[i,] <- prob_1_swap
  all_pi_2_swap[i,] <- prob_2_swap
  all_pi_3_swap[i,] <- prob_3_swap
  all_pi_4_swap[i,] <- prob_4_swap
  
}

x <- seq(1,M_iter)

par(mfrow = c(5,6))
for(i in 1:30){
  plot(x, all_pi_1_swap[,i], type = "l")
}
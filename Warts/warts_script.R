rm(list = ls())

setwd("/Users/dehuabi/Desktop/SAM 1 Final/multi sim/")
#set.seed(1111)
#set.seed(1004)
#set.seed(1009)
#set.seed(1010)
#set.seed(1013)
set.seed(1014)

warts <- read.csv("wart_dataset.csv")

d_mat <- array(dim = c(2, 71, 4))
d_mat[1,,] <- as.matrix(warts[warts$Treatment == 1,1:4])
d_mat[2,1:48,] <- as.matrix(warts[warts$Treatment == 0,1:4])

library(binhf)
source("MPAM_slice.R")
library(mclust)
library(mcclust)
library(mcclust.ext)

# SAM 1 MCMC
M_iter <- 10000
burnin <- 10000

mu0 <- c(0,0,0,0)
sig0 <- diag(nrow = 4)

a <- 0.5
b <- 0.5

kappa0 <- 0.1
nu0 <- 4

a1 <- 3
b1 <- 3

a2 <- 3
b2 <- 3

post_result <- MCMC_SAM1(d_mat, mu0, sig0, a, b, kappa0, nu0, a1, b1, a2, b2, 
                         0.5, 0.3, 0.5, 10, M_iter, burnin)

saveRDS(post_result, file = "post_warts_mcmc.RData")
post_result <- readRDS("post_warts_mcmc.RData")

cls <- matrix(rep(NA, M_iter*(71+48)), nrow = M_iter)
for(iter in (burnin+1):(burnin+M_iter)){
  z_mat_est <- post_result$all_z_mat[[iter]]
  z_mat_vec_est <- as.vector(t(z_mat_est))
  z_mat_vec_est <- z_mat_vec_est[!is.na(z_mat_vec_est)]
  cls[iter-burnin,] <- z_mat_vec_est
}
psm <- comp.psm(cls)
optmal_cluster <- minVI(psm, method="greedy")

optmal_cluster2 <- minVI(psm, cls, method="all")

saveRDS(optmal_cluster, file = "post_warts_mcmc_opt.RData")

print(sort(unique(optmal_cluster$cl)))
print(sort(unique(optmal_cluster$cl[1:71])))
print(sort(unique(optmal_cluster$cl[72:119])))

opt_c <- optmal_cluster2$cl[1,]
print(sort(unique(opt_c)))
print(sort(unique(opt_c[1:71])))
print(sort(unique(opt_c[72:119])))

K <- 50

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

all_pi_p_1 <- matrix(rep(NA, M_iter*K), ncol = K)
all_pi_p_2 <- matrix(rep(NA, M_iter*K), ncol = K)
all_pi_1 <- matrix(rep(NA, M_iter*K), ncol = K)
all_pi_2 <- matrix(rep(NA, M_iter*K), ncol = K)
all_phi1 <- matrix(rep(NA, M_iter*K), ncol = K)
all_phi2 <- matrix(rep(NA, M_iter*K), ncol = K)
all_phi3 <- matrix(rep(NA, M_iter*K), ncol = K)
all_phi4 <- matrix(rep(NA, M_iter*K), ncol = K)
all_sig <- array(dim = c(M_iter, K, 4, 4))
max_length <- 0

for(iter in (burnin+1):(burnin+M_iter)){
  pi_p_mat <- post_result$all_pi_p_mat[[iter-burnin]]
  all_pi_p_1[iter-burnin,1:length(pi_p_mat[1,])] <- pi_p_mat[1,]
  all_pi_p_2[iter-burnin,1:length(pi_p_mat[1,])] <- pi_p_mat[2,]
  all_pi_1[iter-burnin,1:length(pi_p_mat[1,])] <- const_weight(pi_p_mat[1,])
  all_pi_2[iter-burnin,1:length(pi_p_mat[1,])] <- const_weight(pi_p_mat[2,])
  if(length(pi_p_mat[1,]) > max_length){
    max_length <- length(pi_p_mat[1,])
  }
  all_sig[iter-burnin, 1:length(pi_p_mat[1,]),,] <- post_result$all_sig_mat[[iter-burnin]]
  all_phi <- post_result$all_phi_vec[[iter-burnin]]
  all_phi1[iter-burnin,1:length(pi_p_mat[1,])] <- all_phi[,1]
  all_phi2[iter-burnin,1:length(pi_p_mat[1,])] <- all_phi[,2]
  all_phi3[iter-burnin,1:length(pi_p_mat[1,])] <- all_phi[,3]
  all_phi4[iter-burnin,1:length(pi_p_mat[1,])] <- all_phi[,4]
}

x <- seq(1, M_iter)

par(mfrow = c(5,5))
for(i in 1:25){
  plot(x, all_pi_1[,i], type = "l")
}

par(mfrow = c(5,5))
for(i in 1:25){
  plot(x, all_pi_2[,i], type = "l")
}

par(mfrow = c(5,5))
for(i in 1:25){
  plot(x, all_phi1[,i], type = "l")
}

par(mfrow = c(5,5))
for(i in 1:25){
  plot(x, all_phi1[,i], type = "l")
}

par(mfrow = c(5,5))
for(i in 1:25){
  plot(x, all_phi2[,i], type = "l")
}

par(mfrow = c(5,5))
for(i in 1:25){
  plot(x, all_phi3[,i], type = "l")
}

par(mfrow = c(5,5))
for(i in 1:25){
  plot(x, all_phi4[,i], type = "l")
}

for(i in c(1,2,4,5,6,7,10)){#c(2,3,5,6,8,9,10)){
  print(c(mean(all_phi1[1000:M_iter,i], na.rm = T), mean(all_phi2[1000:M_iter,i], na.rm = T), 
          mean(all_phi3[1000:M_iter,i], na.rm = T), mean(all_phi4[1000:M_iter,i], na.rm = T)))
}

for(i in c(1,2,4,5,6,7,10,11)){#c(2,3,5,6,8,9,10)){
  print(c(mean(all_pi_1[1000:M_iter,i], na.rm = T), mean(all_pi_2[1000:M_iter,i], na.rm = T)))
}

for(i in c(1,2,4,5,6,7,11)){#c(2,3,5,6,8,9,10)){
  print(c(apply(all_sig[1000:M_iter,i,,], na.rm = T)))
}

library(scatterplot3d)

d_warts_all <- warts
d_warts_all$cluster <- opt_c #optmal_cluster$cl
d_warts_all$cluster[d_warts_all$cluster == 11] <- 2
d_warts_all$Treatment[d_warts_all$Treatment == 0] <- rep(2,48)

colors <- c("red", "orange", NA, "black", "blue", "green", "purple", NA, NA, "cyan")
colors <- colors[as.numeric(d_warts_all$cluster)]

shapes <- c(16, 17)
shapes <- shapes[as.numeric(d_warts_all$Treatment)]

scatterplot3d(d_warts_all[,1:3], pch = shapes, color=colors, 
              grid=TRUE, box=FALSE)

par(mfrow = c(1,1))

scatterplot3d(d_warts_all[,c(1,2,4)], pch = shapes, color=colors, 
              grid=TRUE, box=FALSE, angle = 45)

#[1] "red"    "black"  "green"  "blue"   "purple" "orange" "cyan"  

plt <- scatterplot3d(d_warts_all[,c(1,3,4)], pch = shapes, color=colors, 
              grid=TRUE, box=FALSE, angle = 35, xlab = "Age", ylab = "NW", zlab = "Area")#,
              #main = "Responders posterior cluster membership")

legend(plt$xyz.convert(48, 20, 50), col= c("green", "red", "blue", "black", "purple", "orange", "cyan"), bg="white", lty=rep(1,7), lwd=2, yjust=0, bty = "n",
       legend = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7"), cex = 1.1)
legend(plt$xyz.convert(48, 20, 540), pch= c(16,17), bg="white", lty=c(0,0), lwd=2, yjust=0, bty = "n",
       legend = c("  Immuno   ", "  Cryo   "), cex = 1.1)


scatterplot3d(d_warts_all[,2:4], pch = shapes, color=colors, 
              grid=TRUE, box=FALSE, angle = 45)




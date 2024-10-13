rm(list = ls())

setwd("/Users/dehuabi/Desktop/PAM_Extra_Sim/FSBP_Sim/")

# Extra Scenario 1
n_sim <- 30

clusters <- 5
gap <- 3
sd <- sqrt(0.6)
n <- 150

true_prob <- rep(1/clusters, clusters)
true_mean <- seq(0,(clusters-1))*gap

all_data <- list()

for(s in 1:n_sim){
  set.seed(12344+s)
  # Generate observations from N(0,1)
  d_vec <- rep(NA, n)
  c_vec <- sample(seq(1,clusters), n, replace = T, prob = true_prob)
  for(k in unique(c_vec)){
    d_vec[which(c_vec == k)] <- rnorm(length(which(c_vec == k)), mean = true_mean[k], sd = sqrt(0.6))
  }
  
  save_data <- list(d_vec = d_vec, c_vec = c_vec)
  all_data[[s]] <- save_data
  
}

saveRDS(all_data, paste0("./data/data_",n,".RData"))



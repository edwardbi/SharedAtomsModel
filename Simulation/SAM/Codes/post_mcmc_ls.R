library(label.switching)

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

process_mcmc <- function(d_mat, num_p_j, post_z_mat, post_w_mat, post_mu, burnin, M_iter, 
                         max_K = 50, mat_flag = T, univ = T){
  
  J <- length(num_p_j)
  if(univ == T){
    cov_p <- 1
  }else{
    cov_p <- dim(post_mu[[1]])[3]
  }
  
  if(univ == T){
    d_mat_vec <- as.vector(d_mat)
    d_mat_vec <- d_mat_vec[!is.na(d_mat_vec)]
  }else{
    d_mat_vec <- matrix(rep(NA, cov_p*sum(num_p_j)), nrow = sum(num_p_j))
    for(i in 1:length(num_p_j)){
      if(i == 1){
        d_mat_vec[1:num_p_j[1],] <- d_mat[i,1:num_p_j[i],]
      }else{
        d_mat_vec[(sum(num_p_j[1:(i-1)])+1):sum(num_p_j[1:i]),] <- d_mat[i,1:num_p_j[i],]
      }
    }
  }
  
  z_mat_mcmc <- matrix(rep(NA, M_iter*sum(num_p_j)), nrow = M_iter)
  mcmc_par <- array(dim = c(M_iter, max_K, J+cov_p))
  max_length <- 0
  
  for(iter in (burnin+1):(burnin+M_iter)){
    z_mat_est <- post_z_mat[[iter]]
    if(mat_flag == T){
      z_mat_mcmc_val <- as.vector(t(z_mat_est))
      z_mat_mcmc_val <- z_mat_mcmc_val[!is.na(z_mat_mcmc_val)]
    }else{
      z_mat_mcmc_val <- z_mat_est
    }
    pi_p_mat <- post_w_mat[[iter]]
    if(length(pi_p_mat[1,]) > max_length){
      max_length <- length(pi_p_mat[1,])
    }
    all_phi <- post_mu[[iter]]
    z_mat_mcmc[iter-burnin,] <- z_mat_mcmc_val
    for(l in 1:J){
      mcmc_par[iter-burnin,1:length(pi_p_mat[1,]),l] <- pi_p_mat[l,]
    }
    if(univ == T){
      mcmc_par[iter-burnin,1:length(pi_p_mat[1,]),(J+1)] <- all_phi
    }else{
      for(l2 in 1:cov_p){
        mcmc_par[iter-burnin,1:length(pi_p_mat[1,]),(l2+J)] <- all_phi[,l2]
      }
    }
  }
  
  if(max_K > max_length){
    mcmc_par <- mcmc_par[,1:max_length,]
  }
  
  ls <- label.switching(method=c("ECR"), zpivot=z_mat_mcmc[1,], z = z_mat_mcmc, 
                        K = max_length, data = d_mat_vec, mcmc = mcmc_par)
  prem_ls <- ls$permutations[["ECR"]]
  
  return(list(z_mat_mcmc = z_mat_mcmc, mcmc_par = mcmc_par, prem_label = prem_ls, J = J, cov_p = cov_p, K = max_length))
  
}

relabel_mcmc <- function(M_iter, K, mcmc_par, J, cov_p, prem_label, w_p = T){
  
  mcmc_par_swap <- mcmc_par
  for(i in 1:M_iter){
    for(j in 1:J){
      if(w_p == T){
        prob_j <- const_weight(mcmc_par[i,,j])
      }else{
        prob_j <- mcmc_par[i,,j]
      }
      mcmc_par_swap[i,,j] <- prob_j[prem_label[i,]]
    }
    for(p in 1:cov_p){
      mcmc_par_swap[i,,(J+p)] <- mcmc_par[i,prem_label[i,],(J+p)]
    }
  }
  
  return(mcmc_par_swap)
  
}



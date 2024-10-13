library(cluster)
library(invgamma)
library(Ckmeans.1d.dp)
library(binhf)
library(mclust)
library(mcclust)
library(mcclust.ext)

set.seed(12345)

# ConstWeight Function
const_weight <- function(w_p){
  
  m_wp <- 1 - w_p
  m_wp2 <- shift(cumprod(m_wp), 1)
  m_wp2[1] <- 1
  
  w <- w_p*m_wp2
  return(w)
  
}

xi_k <- function(kappa,k){
  return((1-kappa)*kappa^(k-1))
}

stoc_trunc <- function(kappa, c_vec){
  
  u_i <- runif(length(c_vec), 0, xi_k(kappa,c_vec))
  K_c <- 1 + floor((log(u_i) - log(1 - kappa))/log(kappa))
  J_c <- max(c_vec)
  K_s <- max(c(K_c, J_c))
  xi_c <- xi_k(kappa, 1:K_s)
  
  return(list(u_vec = u_i, K_s = K_s, xi = xi_c))
  
}

rnorminvgamma <- function(n, mu, nu, alpha, beta){
  Y <- rinvgamma(n, alpha, beta)
  if(n == 1){
    X <- rnorm(1, mu, sd = sqrt(Y/nu))
    return(data.frame(sig = sqrt(Y), mu = X))
  }else{
    X <- rep(NA, n)
    for(i in 1:n){
      X[i] <- rnorm(1, mu, sd = sqrt(Y[i]/nu))
    }
    return(data.frame(sig = sqrt(Y), mu = X))
  }
}

f_phi <- function(x_i, phi_k, sig_k){
  return(dnorm(x_i, mean = phi_k, sd = sig_k))
}

sample_zi <- function(x_i, u_i, pi_p, xi_c, phi, sig, p){
  
  # Construct pi from pi_p
  pi <- const_weight(p*pi_p)
  
  # Find p(z_ji = k|...)
  k_prob <- rep(0, length(pi))
  for(k in 1:length(pi_p)){
    phi_k <- phi[k]
    sig_k <- sig[k]
    pi_k <- pi[k]
    if(pi_k > 0 && (u_i < xi_c[k])){
      k_prob[k] <- log(pi_k) + log(f_phi(x_i, phi_k, sig_k)) - log(xi_c[k])
    }else{
      k_prob[k] <- -Inf
    }
  }
  
  if(!is.finite(max(k_prob))){
    updated_zi <- sample(seq(1,length(pi_p)),1)
  }else{
    k_prob <- exp(k_prob)#exp(k_prob - max(k_prob))
    updated_zi <- sample(seq(1,length(pi_p)),1,prob = k_prob)
  }
  
  return(updated_zi)
  
}

sample_pi_p <- function(z_vec, pi_p, p, alpha, beta, stepsize){
  
  log_post_pi_k_prop <- function(z_vec, k, pi_p_k, alpha, beta){
    
    m_k <- sum(z_vec == k)
    sum_m_s <- sum(z_vec > k)
    
    log_ppkp <- (m_k+alpha-1)*log(p*pi_p_k) + sum_m_s*log(1-p*pi_p_k) + (beta-1)*(1-pi_p_k)
    
    return(log_ppkp)
  }
  
  proposal <- function(pi_p_k_old, stepsize){
    pi_p_k_prop <- runif(1, min = pi_p_k_old - stepsize, max = pi_p_k_old + stepsize)
    if(pi_p_k_prop < 0){
      pi_p_k_prop <- -pi_p_k_prop
    }
    if(pi_p_k_prop > 1){
      pi_p_k_prop <- 2 - pi_p_k_prop
    }
    return(pi_p_k_prop)
  }
  
  MH_step <- function(pi_p, k, p, alpha, beta, stepsize){
    pi_p_k <- pi_p[k]
    # proposal for beta_k_p
    proposal_pi_p_k <- proposal(pi_p_k, stepsize)
    
    # acceptance ratio
    accept_ratio <- log_post_pi_k_prop(z_vec, k, proposal_pi_p_k, alpha, beta) - 
      log_post_pi_k_prop(z_vec, k, pi_p_k, alpha, beta)
    
    ll_pi_p <- log_post_pi_k_prop(z_vec, k, proposal_pi_p_k, alpha, beta)
    ll_pi_c <- log_post_pi_k_prop(z_vec, k, pi_p_k, alpha, beta)
    if(is.infinite(ll_pi_p) && is.infinite(ll_pi_c)){
      accept_ratio <- -Inf
    }
    
    if(log(runif(1)) <= accept_ratio){
      return(list(pi_p_k = proposal_pi_p_k, acc = 1))
    }else{
      return(list(pi_p_k = pi_p_k, acc = 0))
    }
  }
  
  updated_pi_p <- rep(NA, length(pi_p))
  accept_list <- rep(NA, length(pi_p))
  for(k in 1:length(pi_p)){
    updated_pi_p_k_result <- MH_step(pi_p, k, p, alpha, beta, stepsize)
    updated_pi_p_val <- updated_pi_p_k_result$pi_p_k
    if(updated_pi_p_val == 1){
      updated_pi_p_val <- updated_pi_p_val - 1e-4
    }
    updated_pi_p[k] <- updated_pi_p_val
    accept_list[k] <- updated_pi_p_k_result$acc
  }
  
  return(list(pi_p_upd = updated_pi_p, acc_list = accept_list))
  
}

sample_p <- function(p, z_vec, pi_p_vec, a1, b1, stepsize){
  
  log_ll <- function(z_vec, pi_p_vec, p){
    pi_vec <- const_weight(p*pi_p_vec)
    uniq_ks <- unique(z_vec)
    log_prod <- 0
    for(k in uniq_ks){
      log_prod <- log_prod + (sum(z_vec == k))*log(pi_vec[k])
    }
    return(log_prod)
  }
  
  proposal_p <- function(p_old, stepsize){
    p_prop <- runif(1, min = p_old - stepsize, max = p_old + stepsize)
    if(p_prop < 0){
      p_prop <- -p_prop
    }
    if(p_prop > 1){
      p_prop <- 2 - p_prop
    }
    return(p_prop)
  }
  
  p_prop <- proposal_p(p, stepsize)
  
  # acceptance ratio
  accept_ratio <- log_ll(z_vec, pi_p_vec, p_prop) + dbeta(p_prop, a1, b1, log = T) -
    log_ll(z_vec, pi_p_vec, p) - dbeta(p, a1, b1, log = T)
  
  ll_p_p <- log_ll(z_vec, pi_p_vec, p_prop)
  ll_p_c <- log_ll(z_vec, pi_p_vec, p)
  if(is.infinite(ll_p_p) && is.infinite(ll_p_c)){
    accept_ratio <- -Inf
  }
  
  if(log(runif(1)) <= accept_ratio){
    return(list(p = p_prop, acc = 1))
  }else{
    return(list(p = p, acc = 0))
  }
  
}

sample_phi_k <- function(k, z_vec, x_vec, mu0, nu0, a0, b0){
  
  sum_x_i <- sum(x_vec[which(z_vec == k)])
  n <- sum(z_vec == k)
  
  if(n > 0){
    x_bar <- sum_x_i/n
    diff_x_i_x_bar <- 0
    for(i in 1:length(x_vec)){
      if(z_vec[i] == k){
        diff_x_i_x_bar <- diff_x_i_x_bar + (x_vec[i] - x_bar)^2
      }
    }
    
    nu0_upd <- nu0 + n
    a0_upd <- a0 + n/2
    b0_upd <- b0 + diff_x_i_x_bar/2 + ((n*nu0)/nu0_upd)*((x_bar - mu0)^2)/2
    mu0_upd <- ((nu0*mu0) + sum_x_i)/(nu0_upd)
    return(rnorminvgamma(1, mu0_upd, nu0_upd, a0_upd, b0_upd))
  }else{
    return(rnorminvgamma(1, mu0, nu0, a0, b0))
  }
}

sample_beta <- function(z_vec, beta, a2, b2){
  
  n <- length(z_vec)
  c_s <- length(unique(z_vec))
  
  eta <- rbeta(1, beta+1, n)
  w <- rbinom(1, 1, prob = (a2 + c_s - 1)/(n*(b2 - log(eta))))
  if(w == 1){
    beta_upd <- rgamma(1, a2+c_s, b2-log(eta))
  }else{
    beta_upd <- rgamma(1, a2+c_s-1, b2-log(eta))
  }
  
  return(beta_upd)
}

gibbs_FSBP <- function(z_vec, pi_p_vec, phi_vec, sig_vec, p_val, beta_val, x_vec, alpha, kappa, mu0, nu0, a0, b0, a1, b1, a2, b2, stepsize1, stepsize2){
  
  z_vec_upd <- z_vec
  pi_p_vec_upd <- pi_p_vec
  phi_vec_upd <- phi_vec
  sig_vec_upd <- sig_vec 
  p_val_upd <- p_val
  beta_val_upd <- beta_val
  
  # Stochastic Truncation and sample of U_ji
  stoc_trunc_iter <- stoc_trunc(kappa, z_vec_upd)
  u_vec <- stoc_trunc_iter$u_vec
  max_K <- stoc_trunc_iter$K_s
  xi_c <- stoc_trunc_iter$xi
  
  # Update all
  prev_K <- length(pi_p_vec_upd)
  if(prev_K > max_K){
    pi_p_vec_upd <- pi_p_vec_upd[1:max_K]
    phi_vec_upd <- phi_vec_upd[1:max_K]
    sig_vec_upd <- sig_vec_upd[1:max_K]
  }
  if(prev_K < max_K){
    pi_p_vec_upd <- c(pi_p_vec_upd, rbeta(max_K - prev_K, alpha, beta_val_upd))
    rand_NIG <- rnorminvgamma(max_K - prev_K, mu0, nu0, a0, b0)
    phi_vec_upd <- c(phi_vec_upd, rand_NIG$mu)
    sig_vec_upd <- c(sig_vec_upd, rand_NIG$sig)
  }
  
  # Update all pi_p_k
  pi_p_vec_upd_result <- sample_pi_p(z_vec_upd, pi_p_vec_upd, p_val_upd, alpha, beta_val_upd, stepsize1)
  pi_p_vec_upd <- pi_p_vec_upd_result$pi_p_upd
  
  p_val_upd_result <- sample_p(p_val_upd, z_vec_upd, pi_p_vec_upd, a1, b1, stepsize2)
  p_val_upd <- p_val_upd_result$p
  
  # Update all z_is
  for(i in 1:length(x_vec)){
    x_i <- x_vec[i]
    u_i <- u_vec[i]
    z_vec_upd[i] <- sample_zi(x_i, u_i, pi_p_vec_upd, xi_c, phi_vec_upd, sig_vec_upd, p_val_upd)
  }
  
  for(k in 1:max_K){
    sampled_mu_sig <- sample_phi_k(k, z_vec_upd, x_vec, mu0, nu0, a0, b0)
    sig_vec_upd[k] <- sampled_mu_sig$sig
    phi_vec_upd[k] <- sampled_mu_sig$mu
  }
  
  beta_val_upd <- sample_beta(z_vec_upd, beta_val_upd, a1, b1)
  
  return(list(z_vec_upd = z_vec_upd, pi_p_vec_upd = pi_p_vec_upd, 
              phi_vec_upd = phi_vec_upd, sig_vec_upd = sig_vec_upd,
              p_val_upd = p_val_upd, beta_val_upd = beta_val_upd,
              acc = pi_p_vec_upd_result$acc_list))
}

init_vecs <- function(x_vec, n_clust, alpha){
  
  data_vec <- x_vec
  
  z_vec <- kmeans(data_vec, n_clust)$cluster
  
  p <- 0.5
  
  beta <- 1
  
  # Init phi_k
  phi_vec <- rep(NA, n_clust)
  sig_vec <-rep(1, n_clust)
  for(k in 1:n_clust){
    all_data_k <- data_vec[z_vec == k]
    phi_vec[k] <- mean(all_data_k)
  }
  
  # Init beta_k'
  pi_p_vec <- rbeta(n_clust, alpha, beta)
  
  return(list(z_vec = z_vec, pi_p_vec = pi_p_vec, phi_vec = phi_vec, sig_vec = sig_vec, p_val = p, beta_val = beta))
}

MCMC_FSBP <- function(x_vec, mu0, nu0, a0, b0, a1, b1, a2, b2, N_iter, burn_in,
                      alpha = 1, stepsize1 = 0.3, stepsize2 = 0.1, kappa = 0.5, init_clust = 10, iter_num = 100){
  
  init_params <- init_vecs(x_vec, init_clust, alpha)
  
  z_vec_list <- list()
  pi_p_vec_list <- list()
  phi_vec_list <- list()
  sig_vec_list <- list()
  p_list <- rep(NA, burn_in+N_iter+1)
  beta_list <- p_list
  acc_list <- list()
  
  z_vec_list[[1]] <- init_params$z_vec
  pi_p_vec_list[[1]] <- init_params$pi_p_vec
  phi_vec_list[[1]] <- init_params$phi_vec
  sig_vec_list[[1]] <- init_params$sig_vec
  p_list[1] <- init_params$p_val
  beta_list[1] <- init_params$beta_val
  
  for(iter in 1:(burn_in + N_iter)){
    if(iter %% iter_num == 0){
      print(iter) 
    }
    
    gibbs_result <- gibbs_FSBP(z_vec_list[[iter]], pi_p_vec_list[[iter]], phi_vec_list[[iter]], sig_vec_list[[iter]], 
                               p_list[iter], beta_list[iter], x_vec, alpha, kappa, mu0, nu0, a0, b0, a1, b1, a2, b2, stepsize1, stepsize2)
    
    #print(gibbs_result)
    
    z_vec_list[[iter+1]] <- gibbs_result$z_vec_upd
    pi_p_vec_list[[iter+1]] <- gibbs_result$pi_p_vec_upd
    phi_vec_list[[iter+1]] <- gibbs_result$phi_vec_upd
    sig_vec_list[[iter+1]] <- gibbs_result$sig_vec_upd
    p_list[iter+1] <- gibbs_result$p_val
    beta_list[iter+1] <- gibbs_result$beta_val
    acc_list[[iter+1]] <- gibbs_result$acc
  }
  
  return(list(all_z_vec = z_vec_list, all_pi_p_vec = pi_p_vec_list, 
              all_phi_vec = phi_vec_list, all_sig_vec = sig_vec_list,
              all_p = p_list, all_beta = beta_list, all_acc = acc_list))
  
}

nfd <- function(mat_a, mat_b){
  p <- length(mat_a[1,])
  sum_val <- 0
  for(i in 1:p){
    for(j in 1:p){
      sum_val <- sum_val + ((mat_a[i,j] - mat_b[i,j])^2)/(p^2)
    }
  }
  return(sum_val)
}

comp_metrics <- function(fsbp_fit, c_vec, burn_in, N_iter, n){
  
  cls <- matrix(rep(NA, N_iter*n), nrow = N_iter)
  for(iter in (burn_in+1):(burn_in+N_iter)){
    z_vec_est <- fsbp_fit$all_z_vec[[iter]]
    cls[iter-burn_in,] <- z_vec_est
  }
  psm <- comp.psm(cls)
  
  optmal_cluster <- minVI(psm, cls)
  
  rand_idx <- adjustedRandIndex(optmal_cluster$cl, c_vec)
  
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
  ndf_val <- nfd(psm, all_mat_true)
  
  return(list(opt_clus = optmal_cluster, rand_idx = rand_idx, nfd = ndf_val, clusters = length(unique(optmal_cluster$cl))))
  
}







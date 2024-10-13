library(DirichletReg)
library(copula)

stirling <- function(nn){
  
  if (!exists("maxnn", where=parent.frame(2))){
    assign("maxnn", 1, envir=parent.frame(2))
    assign("allss", list(1), envir=parent.frame(2))
  }
  
  
  maxnn.local <- get("maxnn", envir=parent.frame(2))
  
  #only calculate this if needed
  if (nn > maxnn.local){
    allss.local <- get("allss", envir=parent.frame(2))
    
    allss.local[(length(allss.local) + 1):nn] <- 0
    
    for (mm in (maxnn.local + 1):nn){
      allss.local[[mm]] <- c(allss.local[[mm - 1]] * (mm - 1), 0) +
        c(0, allss.local[[mm - 1]])
      mss <- max(allss.local[[mm]])
      allss.local[[mm]] <- allss.local[[mm]] / mss
    }
    
    assign("maxnn", nn, envir=parent.frame(2))
    assign("allss", allss.local, envir=parent.frame(2))
  }
  
  assign("nn", nn, envir=parent.frame(2))
  ss <- eval(quote(allss[[nn]]), envir=parent.frame(2))
  return(ss)
}

f_k_rji <- function(idx, data, z_mat, k, mu0, lambda0, a0, b0){
  x_ji <- data[idx[1],idx[2]]
  x_k <- c()
  n_k <- 0
  for(j in 1:length(z_mat[,1])){
    for(i in 1:length(z_mat[1,])){
      if(j == idx[1] && i == idx[2]){
        next
      }
      val <- z_mat[j,i]
      if(is.na(val)){
        next
      }
      if(val == k){
        x_k <- c(x_k,data[j,i])
        n_k <- n_k + 1
      }
    }
  }
  
  mean_x_k <- mean(x_k)
  mean_x_k_x <- mean(c(x_k,x_ji))
  
  lambda_n <- lambda0 + n_k
  a_n <- a0 + n_k/2
  b_n <- b0 + 0.5*sum((x_k-mean_x_k)^2) + ((n_k*lambda0)/(n_k+lambda0))*((mean_x_k-mu0)^2)/2
  
  if(n_k == 0){
    lambda_n <- lambda0
    a_n <- a0
    b_n <- b0
  }
  
  lambda_n_s <- lambda_n + 1
  a_n_s <- a_n + 1/2
  b_n_s <- b0 + 0.5*(sum((x_k-mean_x_k_x)^2) + (x_ji-mean_x_k_x)^2) + (((n_k+1)*lambda0)/(n_k+1+lambda0)) *
    ((mean_x_k_x-mu0)^2)/2
  
  term1 <- sqrt(lambda_n/lambda_n_s)
  term2 <- exp(a_n*log(b_n)-a_n_s*log(b_n_s))
  term3 <- exp(lgamma(a_n_s) - lgamma(a_n))
  
  #print(c(term1, term2, term3))
  
  ratio <- term1*term2*term3*(1/(2*sqrt(pi)))
  
  return(ratio)
  
}

# Function to update z matrix
sample_z_ji <- function(idx, z_mat, beta_vec, alpha, data, hyperp){
  
  # obtain the current z_ji and z_j.
  z_ji <- z_mat[idx[1],idx[2]]
  z_j <- z_mat[idx[1],]
  z_j <- z_j[!is.na(z_j)]
  
  k_vec_prob <- rep(0, length(beta_vec))
  
  # Compute p(z_ji = k|rest) if k previously used
  for(k in 1:(length(beta_vec)-1)){
    if(k %in% z_j){
      n_jdk_rmji <- sum(z_j == k)
      if(z_ji == k){
        n_jdk_rmji <- n_jdk_rmji - 1
      }
    }else{
      n_jdk_rmji <- 0
    }
    k_vec_prob[k] <- (n_jdk_rmji + alpha*beta_vec[k])*
      f_k_rji(idx, data, z_mat, k, hyperp[1], hyperp[2], hyperp[3], hyperp[4])
  }
  # Compute p(z_ji = k_new|rest)
  k_vec_prob[length(k_vec_prob)] <- 
    alpha*beta_vec[length(beta_vec)]*
    f_k_rji(idx, data, z_mat, length(beta_vec), hyperp[1], hyperp[2], hyperp[3], hyperp[4])
  
  samp_k <- sample(seq(1,length(beta_vec)), size = 1, prob = k_vec_prob)
  return(samp_k)
}

# Add an additional colume to z matrix
add_col_m_mat <- function(z_mat, m_mat){
  max_K <- max(as.vector(z_mat),na.rm = T)
  #print(max_K)
  #print(m_mat)
  if(max_K > length(m_mat[1,])){
    m_mat_upd <- matrix(rep(0, length(m_mat[,1])*(length(m_mat[1,])+1)), ncol = length(m_mat[1,])+1)
    m_mat_upd[,1:(length(m_mat_upd[1,])-1)] <- m_mat
    for(j in 1:length(z_mat[,1])){
      if(max_K %in% z_mat[j,]){
        m_mat_upd[j,length(m_mat_upd[1,])] <- 1
      }
    }
    return(m_mat_upd)
  }else{
    return(m_mat)
  }
}

# If a cluster is not in all J in z_mat, relabel z_mat and remove this cluster in m_mat
bookkeeping_z_m_mats <- function(z_mat, m_mat){
  unique_k <- unique(as.vector(z_mat),na.rm = T)
  unique_k <- sort(unique_k[!is.na(unique_k)])
  max_k <- max(as.vector(z_mat),na.rm = T)
  #print(unique_k)
  #print(max_k)
  if(length(unique_k) < max_k){
    m_mat_upd <- m_mat[,unique_k]
    empty_clust <- seq(1,max_k)[-unique_k]
    z_mat_upd <- z_mat
    for(j in 1:length(z_mat[,1])){
      for(i in 1:length(z_mat[1,])){
        if(!is.na(z_mat_upd[j,i])){
          for(c in sort(empty_clust, decreasing=TRUE)){
            if(z_mat_upd[j,i] > c){
              z_mat_upd[j,i] <- z_mat_upd[j,i] - 1
            }
          }
        }
      }
    }
    return(list(z_mat = z_mat_upd, m_mat = m_mat_upd))
  }else if(length(m_mat[1,]) > length(unique_k)){
    m_mat_upd <- m_mat[,unique_k]
    return(list(z_mat = z_mat, m_mat = m_mat_upd))
  }else{
    return(list(z_mat = z_mat, m_mat = m_mat))
  }
}

# Sample the m matrix
sample_m_jk <- function(idx, z_mat, m_mat, beta_vec, alpha){
  
  m_jk <- m_mat[idx[1],idx[2]]
  z_j <- z_mat[idx[1],]
  k <- idx[2]
  if(k %in% z_j){
    n_jdk <- sum(z_j == k, na.rm = T)
  }else{
    n_jdk <- 0
  }
  if(n_jdk > 1){
    m_vec_prob <- rep(0, n_jdk)
    gamma_frac <- gamma(alpha*beta_vec[k])/gamma(alpha*beta_vec[k] + n_jdk)
    #log_gamma_frac <- log(gamma(alpha*beta_vec[k])) - log(gamma(alpha*beta_vec[k] + n_jdk))
    ss_result <- stirling(n_jdk)
    for(m in 1:length(m_vec_prob)){
      #usig_ss <- ss_result$ss[m]*exp(ss_result$lmss)
      if(gamma_frac == 0){
        m_vec_prob[m] <- ss_result[m]*(alpha*beta_vec[k])^(m)
        #ss_result$ss[m]*(alpha*beta_vec[k])^(m)
      }else{
        m_vec_prob[m] <- gamma_frac*ss_result[m]*(alpha*beta_vec[k])^(m)
        #gamma_frac*ss_result$ss[m]*(alpha*beta_vec[k])^(m)
      }
      #m_vec_prob[m] <- log_gamma_frac + log(ss_result$ss[m]) + ss_result$lmss + m*log(alpha*beta_vec[k])
    }
    #print(n_jdk)
    #print(gamma_frac)
    #print(m_vec_prob)
    #samp_m <- sample(seq(1,n_jdk), size = 1, prob = m_vec_prob)
    samp_m <- sample(seq(1,n_jdk), size = 1, prob = m_vec_prob/sum(m_vec_prob))
  }else if(n_jdk == 1){
    samp_m <- 1
  }else{
    samp_m <- 0
  }
  return(samp_m)
  
}

# Sample the beta vector
sample_beta_vec <- function(m_mat, gamma){
  
  K <- length(m_mat[1,])
  
  k_list <- rep(NA, K+1)
  for(k in 1:K){
    m_k_vec <- m_mat[,k]
    #m_k_vec <- m_k_vec[!is.na(m_k_vec)]
    k_list[k] <- sum(m_k_vec)
  }
  k_list[K+1] <- gamma
  betas <- rdirichlet(1, k_list)
  return(betas)
  
}

# Sample alpha0
sample_alpha0 <- function(z_mat, m_mat, alpha0, a, b){
  J <- length(z_mat[,1])
  w_vec <- rep(NA, J)
  s_vec <- rep(NA, J)
  for(j in 1:J){
    n_j <- length(z_mat[j,!is.na(z_mat[j,])])
    w_vec[j] <- rbeta(1, alpha0+1, n_j)
    s_vec[j] <- runif(1)*(alpha0+n_j) < n_j
  }
  m_dd <- sum(as.vector(m_mat), na.rm = T)
  alpha0_upd <- rgamma(1, a+m_dd-sum(s_vec), b - sum(log(w_vec)))
  return(alpha0_upd)
}

sample_gamma <- function(m_mat, gamma0, a, b){
  m_dd <- sum(as.vector(m_mat), na.rm = T)
  K <- length(m_mat[1,])
  
  w <- rbeta(1, gamma0+1, m_dd)
  weights <- c((a+K-1)/(b-log(w)), m_dd)
  s <- rbinom(1, 1, prob = weights[1]/sum(weights))
  
  if(s == 1){
    gamma_upd <- rgamma(1, a+K, b-log(w))
  }else{
    gamma_upd <- rgamma(1, a+K-1, b-log(w))
  }
  
  return(gamma_upd)
}

init_mat <- function(data, gamma, num_c){
  d_mat <- data
  vec_d <- as.vector(d_mat)
  vec_d <- vec_d[!is.na(vec_d)]
  memship <- kmeans(vec_d, centers=num_c)$cluster
  init_z <- matrix(rep(NA, length(d_mat[1,])*length(d_mat[,1])), ncol = length(d_mat[1,]))
  for(j in 1:length(d_mat[,1])){
    for(i in 1:length(d_mat[1,])){
      if(!is.na(d_mat[j,i])){
        idx <- which(vec_d == d_mat[j,i])
        init_z[j,i] <- memship[idx]
      }
    }
  }
  init_m <- matrix(rep(0, num_c*length(d_mat[,1])), ncol = num_c)
  for(j in 1:length(d_mat[,1])){
    for(k in 1:num_c){
      max_table <- sum(which(init_z[j,] == k))
      if(max_table > 0){
        init_m[j,k] <- 1
      }else{
        init_m[j,k] <- 0
      }
    }
  }
  init_beta_vec_param <- rep(NA, num_c+1)
  for(k in 1:num_c){
    m_dk <- sum(init_m[,k])
    init_beta_vec_param[k] <- m_dk
  }
  init_beta_vec_param[num_c+1] <- gamma
  beta_vec_init <- rdirichlet(1, init_beta_vec_param)
  return(list(z_mat = init_z, m_mat = init_m, beta_vec = beta_vec_init))
}

init_mat_rand <- function(data, gamma, num_c){
  d_mat <- data
  init_z <- matrix(rep(NA, length(d_mat[1,])*length(d_mat[,1])), ncol = length(d_mat[1,]))
  for(j in 1:length(d_mat[,1])){
    for(i in 1:length(d_mat[1,])){
      if(!is.na(d_mat[j,i])){
        init_z[j,i] <- sample(seq(1,num_c), 1)
      }
    }
  }
  init_m <- matrix(rep(0, num_c*length(d_mat[,1])), ncol = num_c)
  for(j in 1:length(d_mat[,1])){
    for(k in 1:num_c){
      max_table <- sum(which(init_z[j,] == k))
      if(max_table > 0){
        init_m[j,k] <- 1
      }else{
        init_m[j,k] <- 0
      }
    }
  }
  init_beta_vec_param <- rep(NA, num_c+1)
  for(k in 1:num_c){
    m_dk <- sum(init_m[,k])
    init_beta_vec_param[k] <- m_dk
  }
  init_beta_vec_param[num_c+1] <- gamma
  beta_vec_init <- rdirichlet(1, init_beta_vec_param)
  return(list(z_mat = init_z, m_mat = init_m, beta_vec = beta_vec_init))
}


gibbs_sampler_iter <- function(data, z_mat, m_mat, beta_vec, alpha, hyperp, gamma, a1, a2, b1, b2){
  z_mat_upd <- z_mat
  m_mat_upd <- m_mat
  beta_vec_upd <- beta_vec
  alpha0_upd <- alpha
  gamma_upd <- gamma
  
  # sample for z
  for(j in 1:length(data[,1])){
    for(i in 1:length(data[1,])){
      if(!is.na(data[j,i])){
        idx <- c(j,i)
        z_mat_upd[j,i] <- sample_z_ji(idx, z_mat_upd, beta_vec, alpha0_upd, data, hyperp)
      }
    }
  }
  m_mat_upd <- add_col_m_mat(z_mat_upd, m_mat_upd)
  zm_mat_upd <- bookkeeping_z_m_mats(z_mat_upd, m_mat_upd)
  z_mat_upd <- zm_mat_upd$z_mat
  m_mat_upd <- zm_mat_upd$m_mat
  
  #print("samped z")
  #print(z_mat_upd)
  #print(m_mat_upd)
  
  # sample for m
  for(j in 1:length(m_mat_upd[,1])){
    for(k in 1:length(m_mat_upd[1,])){
      idx2 <- c(j,k)
      #print(idx2)
      m_mat_upd[j,k] <- sample_m_jk(idx2, z_mat_upd, m_mat_upd, beta_vec, alpha0_upd)
    }
  }
  
  #print("samped m")
  #print(m_mat_upd)
  
  # sample for beta
  beta_vec_upd <- sample_beta_vec(m_mat_upd, gamma_upd)
  
  #print("samped beta")
  
  alpha0_upd <- sample_alpha0(z_mat_upd, m_mat_upd, alpha0_upd, a1, b1)
  
  gamma_upd <- sample_gamma(m_mat_upd, gamma_upd, a2, b2)
  
  return(list(z_mat = z_mat_upd, m_mat = m_mat_upd, beta_vec = beta_vec_upd,
              alpha0 = alpha0_upd, gamma = gamma_upd))
}

mcmc_alg <- function(N_iter, burn_in, data, alpha, gamma, hyperp, a1, a2, b1, b2){
  # Initization
  init_mat_result <- init_mat(data, gamma, 3)
  z_mat <- init_mat_result$z_mat
  m_mat <- init_mat_result$m_mat
  beta_vec <- init_mat_result$beta_vec
  
  z_mat_list <- list()
  m_mat_list <- list()
  beta_vec_list <- list()
  alpha0_list <- list()
  gamma_list <- list()
  
  z_mat_list[[1]] <- z_mat
  m_mat_list[[1]] <- m_mat
  beta_vec_list[[1]] <- beta_vec
  alpha0_list[[1]] <- alpha
  gamma_list[[1]] <- gamma
  
  for(i in 1:(N_iter+burn_in)){
    gibbs_result <- gibbs_sampler_iter(data, z_mat_list[[i]], m_mat_list[[i]], beta_vec_list[[i]], 
                                       alpha0_list[[i]], hyperp, gamma_list[[i]],
                                       a1, a2, b1, b2)
    z_mat_list[[i+1]] <- gibbs_result$z_mat
    m_mat_list[[i+1]] <- gibbs_result$m_mat
    beta_vec_list[[i+1]] <- gibbs_result$beta_vec
    alpha0_list[[i+1]] <- gibbs_result$alpha0
    gamma_list[[i+1]] <- gibbs_result$gamma
    if(i %% 1000 == 0){
      print(i)
      #print(length(gibbs_result$m_mat[1,]))
    }
  }
  
  return(list(z_mat_list = z_mat_list, m_mat_list = m_mat_list, beta_vec_list = beta_vec_list, 
              alpha0 = alpha0_list, gamma = gamma_list))
  
}

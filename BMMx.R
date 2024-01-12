remotes::install_github("emilieodegaard/BayesMallows") # if not already installed
library(BayesMallows)
library(dplyr)
library(tidyr)

set.seed(20032023)

# Data dimensions
n = 20 # no. of items/variables
C = 3 # no. of clusters
N_cl = 30 # no. of people in each cluster
N = N_cl*C # total no. of assessors/people
K = 3 # no of covariates

# Parameters
alpha_true = 5
sigma <- 1 # for now, fix sigma
c1 <- 0.5 # default=0.5
c2 <- 10 # default=10
psi <- 1
cov_params = list(c=c(c1,c2), d=psi, theta=1) # tuning parameters
cov_type <- c(1,1,1) # "1" = continuous covariate, i.e. here  we simulate three 
                     # continuous covariates, and need to specify it to BMMx

# Function for generating the covariate cluster-wise means: mu
mu.eqdist.solver0 <- function(d, C, K) {
  mu <- matrix(data=NA, nrow=C, ncol=K)
  for(c in 1:C){
    mu[c,] <- rep(d*(c-1),K)
  }
  return(mu)
}

# Function for generating the cluster centroids (rho's) by swapping items
rho.swap <- function(n, W, d, C){
  rho = matrix(data=NA, nrow=C, ncol=n)
  for(c in 1:C){
    rho[c,] = c(1:n)
  }
  
  # Do swapping if d is not zero
  if(d != 0){ 
    # Different cluster centroids: sample unique swaps from rho1
    no_swaps = (C-1)*W # this is the total number of swaps needed to generate all the cluster centroids
    
    # Sample items to swap, test if they are valid, keep sampling them if not (maybe need an upper limit here?)
    my_test = T
    while(my_test){ # keep sampling the items to swap 
      
      # Need no_swaps unique swaps at distance d 
      items0 = sample.int(n = n, size = no_swaps) # original items to be swapped with another set of items
      add_sub = sample(x=c(-d,d), size=no_swaps, replace=T)
      
      # Randomly add or subtract d
      items_swapped = items0 + add_sub
      
      # Test if items set is valid
      my_test = isTRUE(max(items_swapped) > n | min(items_swapped) < 1 | any(duplicated(c(items0,items_swapped))))
      
    }
    
    
    # Swap original items in rho1 with the new set 
    for(c in 2:C){ # skip the first cluster centroid as we keep it fixed
      rho_swapped = c(1:n) # start here
      items_rm0 = array(NA, W)
      items_rm1 = array(NA, W)
      for(w in 1:W){
        # Swap back and forth
        rho_swapped[items0[w]] = items_swapped[w] 
        rho_swapped[items_swapped[w]] = items0[w] 
        
        # Save the items that were swapped so we can remove them and avoid any duplicates
        items_rm0[w] = items0[w] 
        items_rm1[w] = items_swapped[w] 
      }
      
      # Remove items that were swapped in the previous cluster
      items_swapped = items_swapped[! items_swapped %in% items_rm1]
      items0 = items0[! items0 %in% items_rm0]
      
      rho[c,] <- rho_swapped
    }
  }
  return(rho)  
}

# True labels
z_true <- c(rep("Cluster 1",N_cl), rep("Cluster 2",N_cl), rep("Cluster 3",N_cl)) 
shuffle <- sample.int(N, size=N)
z_true <- z_true[shuffle]

# MCMC params
M0 = 1e4 # no. of MCMC iterations for BMM
burn_in0 <- 0.1*M0
M1 = 1e4 # no. of MCMC iterations for BMMx
burn_in1 <- 0.1*M1
R <- 10 # no. of runs for d_rho
res <- data.frame()

# Distance parameters for data generation
W = 4 # number of swaps
dist_rs <- seq.int(0,10,2) # rho: cluster centroids
dist_cs <- c(0,1,3,6,9) # mu: covariates 


# Run BMM and BMMX with varying dist for cov and rho
for(dr in 1:length(dist_rs)){
  
  # Generate rho1, rho2, rho3 based on certain distance 
  dist_r <- dist_rs[dr]
  
  rho <- rho.swap(n, W, dist_r, C)
  rho_true <- rho
  
  for(r in 1:R){
    
    # Simulate the dataset based on the rhos and a certain alpha
    cl1 <- sample_mallows(rho0=rho[1,], alpha0=alpha_true, n_samples = N_cl, burnin = 5e3, thinning = 100)
    cl2 <- sample_mallows(rho0=rho[2,], alpha0=alpha_true, n_samples = N_cl, burnin = 5e3, thinning = 100)
    cl3 <- sample_mallows(rho0=rho[3,], alpha0=alpha_true, n_samples = N_cl, burnin = 5e3, thinning = 100)
    data <- rbind(cl1, cl2, cl3)
    data <- data[shuffle,]
    
    seed = sample(c(1:1000),1)
    
    # BMM
    t0 <- proc.time()[[3]] #Sys.time()
    test0 <- compute_mallows(rankings=data, n_clusters=C, save_clus=T, nmc=M0, seed=seed)
    t_end <- proc.time()[[3]] #Sys.time()
    comp_time_bmm <- t_end - t0
    
    for(dc in 1:length(dist_cs)){
        
      # Generate K covariates at a certain euclidian distance for each covariate over C clusters
      dist_x <- dist_cs[dc]
      mu <- mu.eqdist.solver0(dist_x, C, K)
      
      # Save results from BMM run
      res0 <- data.frame(index=numeric(1))
      model = "BMM"
      test0$burnin <- burn_in0
      
      # Confusion table for label assignments
      z.map <- assign_cluster(test0, soft=FALSE)
      z_est <- c(z.map$map_cluster)
      #conf_table <- table(z_true, z_est)
      #labels_correct_assigned <- sum(diag(conf_table))/N
      labels_correct_assigned <- sum(z_true==z_est)/N
      
      # Test for label switching: permute labels, and compute z_est for each, take the best one
      no_unique_z <- length(unique(z_est))
      z_perm <- gtools::permutations(no_unique_z, no_unique_z, unique(z_est))
      no_perm <- nrow(z_perm)
      z_perm_all <- list()
      z_perm_all[[1]] <- z_est
      z_res <- matrix(data=NA, nrow=1, ncol=no_perm)
      z_res[1] <- labels_correct_assigned
      label_switch=F
      if(no_perm > 1){
        for(i in 2:no_perm){
          z_est_tmp <- z_est
          for(j in 1:no_unique_z){
            z_est_tmp[z_est==z_perm[1,j]] <- z_perm[i,j]
          }
          #conf_table <- table(z_true, z_est_tmp)
          #labels_correct_assigned <- sum(diag(conf_table))/N
          labels_correct_assigned <- sum(z_true==z_est)/N
          z_res[i] <- labels_correct_assigned
          z_perm_all[[i]] <- z_est_tmp
          
        }
        labels_correct_assigned <- max(z_res)
        z_est <- z_perm_all[[which.max(z_res)]]
        if(which.max(z_res) != 1){label_switch = T}
        
      }
      
      # Re-label all assignments so that they correspond to the correct labels
      z.map["map_cluster"] <- z_est
      
      # MAP of correct cluster assignment for all assessors
      pp.map <- array(data=0, dim=c(N,1))
      ind.correct <- which(z_true == as.numeric(stringr::str_extract(z.map$map_cluster, "\\d")))
      pp.map[ind.correct,]  <- z.map[ind.correct,]$probability
      mean(pp.map)
      
      # Posterior probabilities of correct cluster assignment for all assessors
      z.pp <- assign_cluster(test0, soft=TRUE)
      if(label_switch == T){
        tmp <- z.pp
        for(c in 1:no_unique_z){
          tmp[which(z.pp["cluster"]== z_perm[1,c]),]$cluster <- z_perm[which.max(z_res),c]
        }
        z.pp <- tmp
      }
      pp.all <- array(data=0, dim= c(N,1))
      for(i in 1:N){
        tmp <- z.pp[z.pp$assessor==i, ]
        tmp <- tmp[as.numeric(stringr::str_extract(tmp$cluster, "\\d"))==z_true[i],]
        if(dim(tmp)[1]!=0){
          pp.all[i] <- tmp$probability
        }
      }
      mean(pp.all)
      
      # All results
      res0$rep <- r
      res0$labels_correct_assigned <- labels_correct_assigned  # total, i.e. the same for all c
      res0$model <- model
      res0$sim_fun_type <- "none" 
      res0$dist_rho <- dist_r
      res0$dist_cov <- dist_x
      res0$alpha <- alpha_true
      res0$comp_time <- comp_time_bmm
      res0$prob_map_mean <- mean(pp.map)
      res0$prob_all_mean <- mean(pp.all)
      res0$prob_map[[1]] <- list(pp.map)
      res0$prob_all[[1]] <- list(pp.all)
      #res0$labels[[c]] <- list(z_est)
      
      
      # Generate the covariates
      set.seed(dc) # we generate the same covariate dataset for each d_x
      Sigma = diag(K)*sigma
      cov_data <- c()
      for(c in 1:C){
        tmp <- MASS::mvrnorm(n = N_cl, mu=mu[c,], Sigma=Sigma)
        cov_data <- rbind(cov_data, tmp)
      }
      # Shuffle the dataset according to the correct labels
      cov_data <- cov_data[shuffle,]
      
      # BMMx: alternative model
      model = "BMMx"
      sim_fun_type = "augmented"
      t0 <- proc.time()[[3]] #Sys.time()
      test1 <- compute_mallows(rankings=data, n_clusters=C, save_clus=T, cov_mat=cov_data, cov_type=cov_type, cov_params=cov_params, prior_type=sim_fun_type, nmc=M1, seed=seed)
      t_end <- proc.time()[[3]] #Sys.time()
      comp_time_bmmx <- t_end - t0
      
      # Save results from run
      res1 <- data.frame(index=numeric(1))
      test1$burnin <- burn_in1
      
      # Confusion table for label assignments
      z.map <- assign_cluster(test1, soft=FALSE)
      z_est <- c(z.map$map_cluster)
      #conf_table <- table(z_true, z_est)
      #labels_correct_assigned <- sum(diag(conf_table))/N
      labels_correct_assigned <- sum(z_true==z_est)/N
      
      # Test for label switching: permute labels, and compute z_est for each, take the best one
      no_unique_z <- length(unique(z_est))
      z_perm <- gtools::permutations(no_unique_z, no_unique_z, unique(z_est))
      no_perm <- nrow(z_perm)
      z_perm_all <- list()
      z_perm_all[[1]] <- z_est
      z_res <- matrix(data=NA, nrow=1, ncol=no_perm)
      z_res[1] <- labels_correct_assigned
      label_switch=F
      if(no_perm > 1){
        for(i in 2:no_perm){
          z_est_tmp <- z_est
          for(j in 1:no_unique_z){
            z_est_tmp[z_est==z_perm[1,j]] <- z_perm[i,j]
          }
          labels_correct_assigned <- sum(z_true==z_est)/N
          z_res[i] <- labels_correct_assigned
          z_perm_all[[i]] <- z_est_tmp
        }
        labels_correct_assigned <- max(z_res)
        z_est <- z_perm_all[[which.max(z_res)]]
        if(which.max(z_res) != 1){label_switch = T}
      }
      z.map["map_cluster"] <- z_est
      
      # MAP of correct cluster assignment for all assessors
      pp.map <- array(data=0, dim=c(N,1))
      ind.correct <- which(z_true == as.numeric(stringr::str_extract(z.map$map_cluster, "\\d")))
      pp.map[ind.correct,]  <- z.map[ind.correct,]$probability
      mean(pp.map, na.rm=T)
      
      # Posterior probabilities of correct cluster assignment for all assessors
      z.pp <- assign_cluster(test1, soft=TRUE)
      if(label_switch == T){
        tmp <- z.pp
        for(c in 1:no_unique_z){
          tmp[which(z.pp["cluster"]== z_perm[1,c]),]$cluster <- z_perm[which.max(z_res),c]
        }
        z.pp <- tmp
      }
      pp.all <- array(data=0, dim= c(N,1))
      for(i in 1:N){
        tmp <- z.pp[z.pp$assessor==i, ]
        tmp <- tmp[as.numeric(stringr::str_extract(tmp$cluster, "\\d"))==z_true[i],]
        if(dim(tmp)[1]!=0){
          pp.all[i] <- tmp$probability
        }
      }
      
      # All results
      res1$rep <- r
      res1$labels_correct_assigned <- labels_correct_assigned  # total, i.e. the same for all c
      res1$model <- model
      res1$sim_fun_type <- sim_fun_type 
      res1$dist_rho <- dist_r
      res1$dist_cov <- dist_x
      res1$alpha <- alpha_true
      res1$comp_time<- comp_time_bmmx
      res1$prob_map_mean <- mean(pp.map)
      res1$prob_all_mean <- mean(pp.all)
      res1$prob_map[[1]] <- list(pp.map)
      res1$prob_all[[1]] <- list(pp.all)
      
      # BMMx: alternative model
      model = "BMMx"
      sim_fun_type = "alternative"
      t0 <- proc.time()[[3]] #Sys.time()
      test2 <- compute_mallows(rankings=data, n_clusters=C, save_clus=T, cov_mat=cov_data, cov_type=cov_type, cov_params=cov_params, prior_type=sim_fun_type, nmc=M1, seed=seed)
      t_end <- proc.time()[[3]] #Sys.time()
      comp_time_bmmx <- t_end - t0
      
      # Save results from run
      res2 <- data.frame(index=numeric(1))
      test2$burnin <- burn_in1
      
      # Confusion table for label assignments
      z.map <- assign_cluster(test2, soft=FALSE)
      z_est <- c(z.map$map_cluster)
      labels_correct_assigned <- sum(z_true==z_est)/N
      
      # Test for label switching: permute labels, and compute z_est for each, take the best one
      no_unique_z <- length(unique(z_est))
      z_perm <- gtools::permutations(no_unique_z, no_unique_z, unique(z_est))
      no_perm <- nrow(z_perm)
      z_perm_all <- list()
      z_perm_all[[1]] <- z_est
      z_res <- matrix(data=NA, nrow=1, ncol=no_perm)
      z_res[1] <- labels_correct_assigned
      label_switch=F
      if(no_perm > 1){
        for(i in 2:no_perm){
          z_est_tmp <- z_est
          for(j in 1:no_unique_z){
            z_est_tmp[z_est==z_perm[1,j]] <- z_perm[i,j]
          }
          #conf_table <- table(z_true, z_est_tmp)
          #labels_correct_assigned <- sum(diag(conf_table))/N
          labels_correct_assigned <- sum(z_true==z_est)/N
          z_res[i] <- labels_correct_assigned
          z_perm_all[[i]] <- z_est_tmp
        }
        labels_correct_assigned <- max(z_res)
        z_est <- z_perm_all[[which.max(z_res)]]
        if(which.max(z_res) != 1){label_switch = T}
      }
      z.map["map_cluster"] <- z_est
      
      # MAP of correct cluster assignment for all assessors
      pp.map <- array(data=0, dim=c(N,1))
      ind.correct <- which(z_true == as.numeric(stringr::str_extract(z.map$map_cluster, "\\d")))
      pp.map[ind.correct,]  <- z.map[ind.correct,]$probability
      mean(pp.map, na.rm=T)
      
      # Posterior probabilities of correct cluster assignment for all assessors
      z.pp <- assign_cluster(test2, soft=TRUE)
      if(label_switch == T){
        tmp <- z.pp
        for(c in 1:no_unique_z){
          tmp[which(z.pp["cluster"]== z_perm[1,c]),]$cluster <- z_perm[which.max(z_res),c]
        }
        z.pp <- tmp
      }
      pp.all <- array(data=0, dim= c(N,1))
      for(i in 1:N){
        tmp <- z.pp[z.pp$assessor==i, ]
        tmp <- tmp[as.numeric(stringr::str_extract(tmp$cluster, "\\d"))==z_true[i],]
        if(dim(tmp)[1]!=0){
          pp.all[i] <- tmp$probability
        }
      }
      
      # All results
      res2$rep <- r
      res2$labels_correct_assigned <- labels_correct_assigned  # total, i.e. the same for all c
      res2$model <- model
      res2$sim_fun_type <- sim_fun_type 
      res2$dist_rho <- dist_r
      res2$dist_cov <- dist_x
      res2$alpha <- alpha_true
      res2$comp_time<- comp_time_bmmx
      res2$prob_map_mean <- mean(pp.map)
      res2$prob_all_mean <- mean(pp.all)
      res2$prob_map[[1]] <- list(pp.map)
      res2$prob_all[[1]] <- list(pp.all)

      # Merge results
      res <- rbind(res, res0, res1, res2)
      
      
    }
  }
}



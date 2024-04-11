## code to reproduce table SM1

## BEGIN SETUP ##
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)

## first try for delta = -4 (original setting with illustrative example)
mu_temp <- -4

## set simulation settings (7 in total for each anticipated difference)
sig1_vec <- c(16.5, rep(18,3), rep(19.5,3))
sig2_vec <- c(16.5, rep(15,3), rep(13,3))
q_vec <- c(1, 1, 15/18, 18/15, 1, 13/19.5, 19.5/13)
q_title <- c("1", rep(c("1", "opt", "inv"), 2))
q_num <- c("1", rep("2", 3), rep("3", 3))
n1s <- seq(2,100)
n_sim <- 1000
n_sob <- 1024
setting <- 1

## simulation_settings
delta <- mu_temp
alpha <- 0.05
Eu <- 19.2
## implement process once for each simulation setting (1 through 7)
for (i in 1:length(sig1_vec)){
  print(i)
  n2s <- pmax(2,round(q_vec[i]*n1s))
  
  sigma_1 <- sig1_vec[i]
  sigma_2 <- sig2_vec[i]
  results <- NULL
  
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n_sim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  ## repeat 1000 times
  results <- foreach(j=1:n_sim, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts) %dopar% {
  
  ## generate a Sobol' sequence of length m
  sob <- sobol(n_sob, d = 3, randomize = "digital.shift", seed = (1000 + i*n_sim + j))
  check1 <- rep(1,n_sob)
  max_sd <- rep(0, n_sob)
  max_sd_samp <- rep(2, n_sob)
  issues_mat <- NULL
  
  ## compute se and lambda functions for each value of n
  for (k in 1:length(n1s)){
    n1 <- n1s[k]; n2 <- n2s[k]
    x <- qchisq(sob[,1], n1 - 1)
    y <- qchisq(sob[,2], n2 - 1)
    z <- qnorm(sob[,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2))
    
    sdv <- sqrt(sigma_1^2*x/((n1 - 1)*n1) + sigma_2^2*y/(n2*(n2 - 1)))
    
    dfw <- sdv^4/(sigma_1^4*x^2/((n1 - 1)^3*n1^2) + sigma_2^4*y^2/((n2 - 1)^3*n2^2))
    
    thres <- (Eu - abs(z))/qt(1-alpha, dfw)
    
    assign(paste0("check",n1),ifelse(sdv <= thres,1,0))
    
    ## this code checks whether points that departed from the rejection region
    ## have rentered to get the average duration column
    unresolved_inds <- which(issues_mat[,3] < 0)
    if (length(unresolved_inds) > 0){
      for (kk in unresolved_inds){
        if (get(paste0("check",n1))[issues_mat[kk,1]] == 1){
          issues_mat[kk,3] <- n1
        }
      }
    }
    ## check to see whether new points have left the rejection region
    check_temp <- ifelse(get(paste0("check",n1-1)), 1 - get(paste0("check",n1)),0)
    inds_issue_temp <- which(check_temp == 1)
    if (length(inds_issue_temp) > 0 & n1 > 2){
      ## negative 1 in third column means that the point has not returned to the rejection
      ## region; it will eventually be replaced with the sample size at which the point returns
      issues_mat <- rbind(issues_mat, cbind(inds_issue_temp, rep(n1,length(inds_issue_temp)), 
                          rep(-1,length(inds_issue_temp))))
    }
    
    ## determine the value of n at which the se function peaks
    max_sd <- pmax(max_sd, sdv)
    max_sd_samp <- ifelse(max_sd == sdv, n1, max_sd_samp)
    
  }
  ## print summary statistics for table SM1
  avg_samp <- mean(max_sd_samp); gt_5 <- mean(ifelse(max_sd_samp > 5, 1, 0)); 
  gt_10 <- mean(ifelse(max_sd_samp > 10, 1, 0)); gt_15 <- mean(ifelse(max_sd_samp > 15, 1, 0))
  gt_20 <- mean(ifelse(max_sd_samp > 20, 1, 0)); gt_25 <- mean(ifelse(max_sd_samp > 25, 1, 0))
  num_issues <- length(unique(issues_mat[,1]))
  
  ## output different results to indicate whether or not there were multiple intersections
  ## num_issues = 0 implies that there were no multiple intersections
  if (num_issues == 0){
    results_temp <- data.frame(num_issues = num_issues, avg_depart = NA, avg_duration = NA, avg_samp = avg_samp,
                  gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
  } else{
    results_temp <- data.frame(num_issues = num_issues, avg_depart = mean(issues_mat[,2]), 
                      avg_duration = mean(issues_mat[,3] - issues_mat[,2]), avg_samp = avg_samp,
                      gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
  }
  
  results_temp
  }
  
  write.csv(results, paste0("SM1_sim", q_num[i], "_q",q_title[i],"_mu", -1*mu_temp, ".csv"), row.names = FALSE)
}

## use this code with all anticipated mean differences to make that section of the table
results_tab <- NULL
for (i in 1:length(q_num)){
  full_res <- read.csv(paste0("SM1_sim", q_num[i], "_q",q_title[i],"_mu", -1*mu_temp, ".csv"))
  res1 <- colMeans(full_res)[c(1, 4, 5, 6)]
  full_res_expand <- full_res[,1]*full_res[,c(2,3)]
  res2 <- colSums(full_res_expand, na.rm = TRUE)/sum(full_res[,1])
  results_tab <- rbind(results_tab, c(res1[1], res2, res1[c(2,3,4)]))
}

write.csv(results_tab, paste0("v2table_SM1_mu", -1*mu_temp, ".csv"), row.names = FALSE)

## repeat simulation with anticipated difference of 0
## code is the same as the case for mu_temp = -4, so see that case for comments
mu_temp <- 0
n1s <- seq(2,100)

delta <- mu_temp
alpha <- 0.05
Eu <- 19.2
for (i in 1:length(sig1_vec)){
  print(i)
  n2s <- pmax(2,round(q_vec[i]*n1s))
  
  sigma_1 <- sig1_vec[i]
  sigma_2 <- sig2_vec[i]
  results <- NULL
  
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n_sim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  results <- foreach(j=1:n_sim, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts) %dopar% {

    sob <- sobol(n_sob, d = 3, randomize = "digital.shift", seed = (2000 + i*n_sim + j))
    check1 <- rep(1,n_sob)
    max_sd <- rep(0, n_sob)
    max_sd_samp <- rep(2, n_sob)
    issues_mat <- NULL
    
    for (k in 1:length(n1s)){
      n1 <- n1s[k]; n2 <- n2s[k]
      x <- qchisq(sob[,1], n1 - 1)
      y <- qchisq(sob[,2], n2 - 1)
      z <- qnorm(sob[,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2))
      
      sdv <- sqrt(sigma_1^2*x/((n1 - 1)*n1) + sigma_2^2*y/(n2*(n2 - 1)))
      
      dfw <- sdv^4/(sigma_1^4*x^2/((n1 - 1)^3*n1^2) + sigma_2^4*y^2/((n2 - 1)^3*n2^2))
      
      thres <- (Eu - abs(z))/qt(1-alpha, dfw)
      
      assign(paste0("check",n1),ifelse(sdv <= thres,1,0))
      
      unresolved_inds <- which(issues_mat[,3] < 0)
      if (length(unresolved_inds) > 0){
        for (kk in unresolved_inds){
          if (get(paste0("check",n1))[issues_mat[kk,1]] == 1){
            issues_mat[kk,3] <- n1
          }
        }
      }
      check_temp <- ifelse(get(paste0("check",n1-1)), 1 - get(paste0("check",n1)),0)
      inds_issue_temp <- which(check_temp == 1)
      if (length(inds_issue_temp) > 0 & n1 > 2){
        issues_mat <- rbind(issues_mat, cbind(inds_issue_temp, rep(n1,length(inds_issue_temp)), 
                                              rep(-1,length(inds_issue_temp))))
      }

      max_sd <- pmax(max_sd, sdv)
      max_sd_samp <- ifelse(max_sd == sdv, n1, max_sd_samp)
      
    }

    avg_samp <- mean(max_sd_samp); gt_5 <- mean(ifelse(max_sd_samp > 5, 1, 0)); 
    gt_10 <- mean(ifelse(max_sd_samp > 10, 1, 0)); gt_15 <- mean(ifelse(max_sd_samp > 15, 1, 0))
    gt_20 <- mean(ifelse(max_sd_samp > 20, 1, 0)); gt_25 <- mean(ifelse(max_sd_samp > 25, 1, 0))
    num_issues <- length(unique(issues_mat[,1]))
    
    if (num_issues == 0){
      results_temp <- data.frame(num_issues = num_issues, avg_depart = NA, avg_duration = NA, avg_samp = avg_samp,
                                 gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
    } else{
      results_temp <- data.frame(num_issues = num_issues, avg_depart = mean(issues_mat[,2]), 
                                 avg_duration = mean(issues_mat[,3] - issues_mat[,2]), avg_samp = avg_samp,
                                 gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
    }
    
    results_temp
  }
  
  write.csv(results, paste0("SM1_sim", q_num[i], "_q",q_title[i],"_mu", -1*mu_temp, ".csv"), row.names = FALSE)
}

## repeat simulation with anticipated difference of -8
## code is the same as the case for mu_temp = -4, so see that case for comments
mu_temp <- -8
n1s <- seq(2,200)

delta <- mu_temp
alpha <- 0.05
Eu <- 19.2
for (i in 1:length(sig1_vec)){
  print(i)
  n2s <- pmax(2,round(q_vec[i]*n1s))
  
  sigma_1 <- sig1_vec[i]
  sigma_2 <- sig2_vec[i]
  results <- NULL
  
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n_sim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  results <- foreach(j=1:n_sim, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts) %dopar% {

    sob <- sobol(n_sob, d = 3, randomize = "digital.shift", seed = (3000 + i*n_sim + j))
    check1 <- rep(1,n_sob)
    max_sd <- rep(0, n_sob)
    max_sd_samp <- rep(2, n_sob)
    issues_mat <- NULL
    
    for (k in 1:length(n1s)){
      n1 <- n1s[k]; n2 <- n2s[k]
      x <- qchisq(sob[,1], n1 - 1)
      y <- qchisq(sob[,2], n2 - 1)
      z <- qnorm(sob[,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2))
      
      sdv <- sqrt(sigma_1^2*x/((n1 - 1)*n1) + sigma_2^2*y/(n2*(n2 - 1)))
      
      dfw <- sdv^4/(sigma_1^4*x^2/((n1 - 1)^3*n1^2) + sigma_2^4*y^2/((n2 - 1)^3*n2^2))
      
      thres <- (Eu - abs(z))/qt(1-alpha, dfw)
      
      assign(paste0("check",n1),ifelse(sdv <= thres,1,0))
      
      unresolved_inds <- which(issues_mat[,3] < 0)
      if (length(unresolved_inds) > 0){
        for (kk in unresolved_inds){
          if (get(paste0("check",n1))[issues_mat[kk,1]] == 1){
            issues_mat[kk,3] <- n1
          }
        }
      }
      check_temp <- ifelse(get(paste0("check",n1-1)), 1 - get(paste0("check",n1)),0)
      inds_issue_temp <- which(check_temp == 1)
      if (length(inds_issue_temp) > 0 & n1 > 2){
        issues_mat <- rbind(issues_mat, cbind(inds_issue_temp, rep(n1,length(inds_issue_temp)), 
                                              rep(-1,length(inds_issue_temp))))
      }
      
      max_sd <- pmax(max_sd, sdv)
      max_sd_samp <- ifelse(max_sd == sdv, n1, max_sd_samp)
      
    }
    
    avg_samp <- mean(max_sd_samp); gt_5 <- mean(ifelse(max_sd_samp > 5, 1, 0)); 
    gt_10 <- mean(ifelse(max_sd_samp > 10, 1, 0)); gt_15 <- mean(ifelse(max_sd_samp > 15, 1, 0))
    gt_20 <- mean(ifelse(max_sd_samp > 20, 1, 0)); gt_25 <- mean(ifelse(max_sd_samp > 25, 1, 0))
    num_issues <- length(unique(issues_mat[,1]))
    
    if (num_issues == 0){
      results_temp <- data.frame(num_issues = num_issues, avg_depart = NA, avg_duration = NA, avg_samp = avg_samp,
                                 gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
    } else{
      results_temp <- data.frame(num_issues = num_issues, avg_depart = mean(issues_mat[,2]), 
                                 avg_duration = mean(issues_mat[,3] - issues_mat[,2]), avg_samp = avg_samp,
                                 gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
    }
    
    results_temp
  }
  
  write.csv(results, paste0("SM1_sim", q_num[i], "_q",q_title[i],"_mu", -1*mu_temp, ".csv"), row.names = FALSE)
}

## repeat simulation with anticipated difference of -12
## code is the same as the case for mu_temp = -4, so see that case for comments
mu_temp <- -12
n1s <- seq(2,500)

delta <- mu_temp
alpha <- 0.05
Eu <- 19.2
for (i in 1:length(sig1_vec)){
  print(i)
  n2s <- pmax(2,round(q_vec[i]*n1s))
  
  sigma_1 <- sig1_vec[i]
  sigma_2 <- sig2_vec[i]
  results <- NULL
  
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n_sim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  results <- foreach(j=1:n_sim, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts) %dopar% {
 
    sob <- sobol(n_sob, d = 3, randomize = "digital.shift", seed = (4000 + i*n_sim + j))
    check1 <- rep(1,n_sob)
    max_sd <- rep(0, n_sob)
    max_sd_samp <- rep(2, n_sob)
    issues_mat <- NULL
    
    for (k in 1:length(n1s)){
      n1 <- n1s[k]; n2 <- n2s[k]
      x <- qchisq(sob[,1], n1 - 1)
      y <- qchisq(sob[,2], n2 - 1)
      z <- qnorm(sob[,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2))
      
      sdv <- sqrt(sigma_1^2*x/((n1 - 1)*n1) + sigma_2^2*y/(n2*(n2 - 1)))
      
      dfw <- sdv^4/(sigma_1^4*x^2/((n1 - 1)^3*n1^2) + sigma_2^4*y^2/((n2 - 1)^3*n2^2))
      
      thres <- (Eu - abs(z))/qt(1-alpha, dfw)
      
      assign(paste0("check",n1),ifelse(sdv <= thres,1,0))
      
      unresolved_inds <- which(issues_mat[,3] < 0)
      if (length(unresolved_inds) > 0){
        for (kk in unresolved_inds){
          if (get(paste0("check",n1))[issues_mat[kk,1]] == 1){
            issues_mat[kk,3] <- n1
          }
        }
      }
      check_temp <- ifelse(get(paste0("check",n1-1)), 1 - get(paste0("check",n1)),0)
      inds_issue_temp <- which(check_temp == 1)
      if (length(inds_issue_temp) > 0 & n1 > 2){
        issues_mat <- rbind(issues_mat, cbind(inds_issue_temp, rep(n1,length(inds_issue_temp)), 
                                              rep(-1,length(inds_issue_temp))))
      }
      
      max_sd <- pmax(max_sd, sdv)
      max_sd_samp <- ifelse(max_sd == sdv, n1, max_sd_samp)
      
    }
    
    avg_samp <- mean(max_sd_samp); gt_5 <- mean(ifelse(max_sd_samp > 5, 1, 0)); 
    gt_10 <- mean(ifelse(max_sd_samp > 10, 1, 0)); gt_15 <- mean(ifelse(max_sd_samp > 15, 1, 0))
    gt_20 <- mean(ifelse(max_sd_samp > 20, 1, 0)); gt_25 <- mean(ifelse(max_sd_samp > 25, 1, 0))
    num_issues <- length(unique(issues_mat[,1]))
    
    if (num_issues == 0){
      results_temp <- data.frame(num_issues = num_issues, avg_depart = NA, avg_duration = NA, avg_samp = avg_samp,
                                 gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
    } else{
      results_temp <- data.frame(num_issues = num_issues, avg_depart = mean(issues_mat[,2]), 
                                 avg_duration = mean(issues_mat[,3] - issues_mat[,2]), avg_samp = avg_samp,
                                 gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
    }
    
    results_temp
  }
  
  write.csv(results, paste0("SM1_sim", q_num[i], "_q",q_title[i],"_mu", -1*mu_temp, ".csv"), row.names = FALSE)
}

## repeat simulation with anticipated difference of -16
## code is the same as the case for mu_temp = -4, so see that case for comments
mu_temp <- -16
n1s <- seq(2,2500)

delta <- mu_temp
alpha <- 0.05
Eu <- 19.2
for (i in 2:length(sig1_vec)){
  print(i)
  n2s <- pmax(2,round(q_vec[i]*n1s))
  
  sigma_1 <- sig1_vec[i]
  sigma_2 <- sig2_vec[i]
  results <- NULL
  
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n_sim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  results <- foreach(j=1:n_sim, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts) %dopar% {

                       sob <- sobol(n_sob, d = 3, randomize = "digital.shift", seed = (i*n_sim + j))
                       check1 <- rep(1,n_sob)
                       max_sd <- rep(0, n_sob)
                       max_sd_samp <- rep(2, n_sob)
                       issues_mat <- NULL
                       
                       for (k in 1:length(n1s)){
                         n1 <- n1s[k]; n2 <- n2s[k]
                         x <- qchisq(sob[,1], n1 - 1)
                         y <- qchisq(sob[,2], n2 - 1)
                         z <- qnorm(sob[,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2))
                         
                         sdv <- sqrt(sigma_1^2*x/((n1 - 1)*n1) + sigma_2^2*y/(n2*(n2 - 1)))
                         
                         dfw <- sdv^4/(sigma_1^4*x^2/((n1 - 1)^3*n1^2) + sigma_2^4*y^2/((n2 - 1)^3*n2^2))
                         
                         thres <- (Eu - abs(z))/qt(1-alpha, dfw)
                         
                         assign(paste0("check",n1),ifelse(sdv <= thres,1,0))
                         
                         ## check plot
                         unresolved_inds <- which(issues_mat[,3] < 0)
                         if (length(unresolved_inds) > 0){
                           for (kk in unresolved_inds){
                             if (get(paste0("check",n1))[issues_mat[kk,1]] == 1){
                               issues_mat[kk,3] <- n1
                             }
                           }
                         }
                         check_temp <- ifelse(get(paste0("check",n1-1)), 1 - get(paste0("check",n1)),0)
                         inds_issue_temp <- which(check_temp == 1)
                         if (length(inds_issue_temp) > 0 & n1 > 2){
                           issues_mat <- rbind(issues_mat, cbind(inds_issue_temp, rep(n1,length(inds_issue_temp)), 
                                                                 rep(-1,length(inds_issue_temp))))
                         }
                         
                         max_sd <- pmax(max_sd, sdv)
                         max_sd_samp <- ifelse(max_sd == sdv, n1, max_sd_samp)
                         
                       }
                       
                       avg_samp <- mean(max_sd_samp); gt_5 <- mean(ifelse(max_sd_samp > 5, 1, 0)); 
                       gt_10 <- mean(ifelse(max_sd_samp > 10, 1, 0)); gt_15 <- mean(ifelse(max_sd_samp > 15, 1, 0))
                       gt_20 <- mean(ifelse(max_sd_samp > 20, 1, 0)); gt_25 <- mean(ifelse(max_sd_samp > 25, 1, 0))
                       num_issues <- length(unique(issues_mat[,1]))
                       
                       if (num_issues == 0){
                         results_temp <- data.frame(num_issues = num_issues, avg_depart = NA, avg_duration = NA, avg_samp = avg_samp,
                                                    gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
                       } else{
                         results_temp <- data.frame(num_issues = num_issues, avg_depart = mean(issues_mat[,2]), 
                                                    avg_duration = mean(issues_mat[,3] - issues_mat[,2]), avg_samp = avg_samp,
                                                    gt_5 = gt_5, gt_10 = gt_10, gt_15 = gt_15, gt_20 = gt_20, gt_25 = gt_25)
                       }
                       
                       results_temp
                     }
  
  write.csv(results, paste0("SM1_sim", q_num[i], "_q",q_title[i],"_mu", -1*mu_temp, ".csv"), row.names = FALSE)
}
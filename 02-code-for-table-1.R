## code to reproduce Table 1

## BEGIN SETUP ##
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)

## code to implement Algorithm 1
## n1 and n2 are sample sizes
## delta is the anticipated difference mu_1 - mu_2
## sigma_1 and sigma_2 are the anticipated standard deviations
## delta_L and delta_U are the interval endpoints
## alpha is the significance level
## m is the number of simulation repetitions
## seed can be set to ensure reproducibility
alg1 <- function(n1, n2, delta, sigma_1, sigma_2, deltaL, deltaU, alpha, m, seed){
  sob <- sobol(m, d = 3, randomize = "digital.shift", seed = seed)
  x <- qchisq(sob[,1], n1 - 1)
  y <- qchisq(sob[,2], n2 - 1)
  z <- qnorm(sob[,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2))
  
  sdv <- sqrt(sigma_1^2*x/((n1 - 1)*n1) + sigma_2^2*y/(n2*(n2 - 1)))
  
  dfw <- sdv^4/(sigma_1^4*x^2/((n1 - 1)^3*n1^2) + sigma_2^4*y^2/((n2 - 1)^3*n2^2))
  
  if (deltaU == Inf){
    thres <- (z - deltaL)/stats::qt(1-alpha, dfw)
  }
  else if (deltaL == -Inf){
    thres <- (deltaU - z)/stats::qt(1-alpha, dfw)
  }
  else{
    thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    thres <- pmin(thresUp, thresLow)
  }
  
  return(mean(ifelse(sdv <= thres,1,0)))
}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 100, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## define parameters for illustrative example
delta <- -4; sigma_1 <- 18; sigma_2 <- 15 
deltaL <- -19.2; deltaU <- 19.2; alpha = 0.05
samps <- c(3,5,8,10,15,20,30,40,50,60)
pwr_sob <- NULL
for (j in 1:length(samps)){
  pwr_temp <- foreach(i=1:100, .combine='c', .packages = c("qrng"),
          .options.snow=opts, .errorhandling = "remove") %dopar% {
            alg1(samps[j], samps[j], delta, sigma_1, sigma_2, deltaL, deltaU, alpha, 65536, i)
          }
  pwr_sob <- rbind(pwr_sob, pwr_temp)
}

write.csv(pwr_sob, "pwr_sob.csv", row.names = FALSE)

## code to implement naive simulation
## inputs are the same as alg1()
alg1Naive <- function(n1, n2, delta, sigma_1, sigma_2, deltaL, deltaU, alpha, m, seed){
  res <- NULL
  for (i in 1:m){
    x1 <- rnorm(n1, 0, sigma_1)
    x2 <- rnorm(n2, delta, sigma_2)
    
    t1 <- t.test(x1, x2, "greater", mu = deltaL, conf.level = 1 - 2*alpha)$p.value
    t2 <- t.test(x1, x2, "less", mu = deltaU, conf.level = 1 - 2*alpha)$p.value
    
    res[i] <- ifelse(t1 < alpha, t2 < alpha,0)
  }
  return(mean(res))
}

pwr_naive <- NULL
for (j in 1:length(samps)){
  print(samps[j])
  pwr_temp <- foreach(i=1:100, .combine='c', .packages = c("qrng"),
                      .options.snow=opts, .errorhandling = "remove") %dopar% {
                        alg1Naive(samps[j], samps[j], delta, sigma_1, sigma_2, deltaL, deltaU, alpha, 65536, i + 100)
                      }
  pwr_naive <- rbind(pwr_naive, pwr_temp)
  print(mean(pwr_temp))
  print(sd(pwr_temp))
}

write.csv(pwr_naive, "pwr_naive.csv", row.names = FALSE)

## return data for Table 1

## Alg 1
as.numeric(round(rowMeans(pwr_sob),4))
as.numeric(apply(pwr_sob, 1, sd))

## Naive simulation
as.numeric(round(rowMeans(pwr_naive),4))
as.numeric(apply(pwr_naive, 1, sd))
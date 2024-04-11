## code to implement algorithm 2

## BEGIN SETUP ##
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)

#' Power Calculations for Two-Group Equivalence Tests with Unequal Variances
#'
#' Approximates the power of equivalence and one-sided hypothesis tests for two independent samples with unequal variances.
#'
#' @param diff the anticipated difference between the group means (\eqn{\mu_{1} - \mu_{2}}).
#' @param sigma1 the anticipated within-group standard deviation for group 1 (\eqn{\sigma_{1}}).
#' @param sigma2 the anticipated within-group standard deviation for group 2 (\eqn{\sigma_{2}}).
#' @param deltaL the lower bound for the interval of equivalence (can be `-Inf` for noninferiority test).
#' @param deltaU the upper bound for the interval of equivalence (can be `Inf` for noninferiority test).
#' @param alpha the significance level \eqn{\alpha \in (0,1)}.
#' @param targetPower the desired statistical power of the equivalence or noninferiority test, must be a single number between 0 and 1 (exclusive). Exactly one of the following pairs must be specified: `targetPower` and `q` or `n1` and `n2`. Specify both `targetPower` and `q` to find sample sizes `n1` and `n2` that yield desired power such that \eqn{n_{2} \approx q n_{1}}.
#' @param q the multiplicative constant for imbalanced two-group sample size determination (\eqn{n_{2} = q n_{1}}). The default value is 1.
#' @param plot a logical variable indicating whether to return a plot of the power curve. If `n1` and `n2` are specified instead of `q` and `targetPower`, this variable is automatically set to `FALSE`. If you wish to approximate many power curves, suppressing the plots will expedite this process.
#' @param seed if provided, a single positive integer is used to ensure reproducibility when randomizing the Sobol' sequence via `sobol()` in the `qrng` package.
#' @param sobol one of the following integers: \eqn{s \in \{0, 1, 2, 3, 4 \}}. When approximating the power curve using `targetPower` and `q`, \eqn{2^{s + 10}} points are generated from the Sobol' sequence. When estimating power for given samples sizes `n1` and `n2`, \eqn{2^{s + 16}} points are generated from the Sobol' sequence. The default setting is \eqn{s = 0}, which ensures that each function call should take less than two seconds. As \eqn{s} increases, the sample size calculation is less sensitive to simulation variability but takes longer to complete. However, all function calls should still take less than 30 seconds when \eqn{s = 4}.
#' @param check a logical variable indicating whether to return a check for incorrect inputs.
#' @examples
#' # specify targetPower and q to obtain sample sizes n1 and n2
#' alg2(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2,
#' alpha = 0.05, targetPower = 0.8, q = 1, plot = TRUE, seed = 1, sobol = 0)
#'
#'
#' @return The sample sizes `n1` and `n2` are returned as a list with supplementary information.
#' @export

## the sample sizes returned by the root-finding algorithm are the first item in the returned list.
## the sample size recommendation for group 1 is returned in the second item.
## the third item is a binary vector dictating how many times the root finding algorithm needed to be reimplemented
alg2 <- function(diff = NULL, sigma1 = NULL, sigma2 = NULL, deltaL = -Inf,
                 deltaU = Inf, alpha = NULL, targetPower = NULL, q = 1, 
                 plot = TRUE, seed = NULL, sobol = 0, check = FALSE){

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  if (check == TRUE){
    ## error checking
    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}
    if(!is.numeric(sigma1) | length(sigma1) != 1){
      stop("Please specify a valid number for sigma1.")}
    else if (sigma1 <= 0 | !is.finite(sigma1)){
      stop("Please specify a valid number for sigma1.")}
    if(!is.numeric(sigma2) | length(sigma2) != 1){
      stop("Please specify a valid number for sigma2.")}
    else if (sigma2 <= 0 | !is.finite(sigma2)){
      stop("Please specify a valid number for sigma2.")}
    if(!is.numeric(deltaL) | length(deltaL) != 1){
      stop("Please specify a valid number for deltaL.")}
    if(!is.numeric(deltaU) | length(deltaU) != 1){
      stop("Please specify a valid number for deltaU.")}
    if(deltaL == -Inf & deltaU == Inf){
      stop("Please specify valid interval endpoints deltaL and deltaU.")}
    if (deltaL >= deltaU){
      stop("Please specify valid interval endpoints deltaL and deltaU.")}
    if(!is.numeric(alpha) | length(alpha) != 1) {
      stop("Please specify a valid number for alpha.")}
    if (is.numeric(alpha)){
      if (alpha <= 0 | alpha >= 1){
        stop("Please specify a valid number for alpha.")}
    }
    if(!is.numeric(targetPower)) {
      stop("Please specify a valid number for targetPower.")}
    if (is.numeric(targetPower)){
      if (targetPower <= 0 | targetPower >= 1){
        stop("Please specify a valid number for targetPower.")}
      if(!is.numeric(q)) {
        stop("Please specify a valid number for q.")}
      else if (is.numeric(q)){
        if (q <= 0) {
          stop("Please specify a valid number for q.")}
      }
      if (diff >= deltaU | diff <= deltaL){
        stop("Please ensure diff is between deltaL and deltaU.")
      }
    }
    if(!is.null(seed) & (!is.numeric(seed) | length(seed) != 1)) {
      stop("Please specify a valid seed for random number generation.")}
    else if (!is.null(seed)){
      if (seed%%1 != 0 | seed < 0){
        stop("Please specify a valid seed for random number generation.")}
    }
    if(!is.numeric(sobol) | length(sobol) != 1) {
      stop("Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")}
    else if (sobol < 0){
      stop("Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")}
    else if (!(sobol %in% c(0,1,2,3,4))){
      stop("Please specify a valid integer between 0 and 4 to initialize the Sobol' sequence.")}
    if (length(plot) != 1 | !(plot %in% c(TRUE, FALSE))){
      stop("Please provide valid logical input for plot.")
    }
  }

  if (length(targetPower) == 1 & length(q) == 1){
    ## this function returns a positive result if the sample corresponds to the rejection region
    targetPowerfn <- function(u, deltaL, deltaU, diff, sigma1, sigma2, n_val, alpha, q){
      n1 <- n_val[1]; n2 <- max(2, q*n_val)
      x <- stats::qchisq(u[1], n1 - 1)
      y <- stats::qchisq(u[2], n2 - 1)
      z <- stats::qnorm(u[3], diff, sqrt(sigma1^2/n1 + sigma2^2/n2))

      sdv <- sqrt(sigma1^2*x/((n1 - 1)*n1) + sigma2^2*y/(n2*(n2 - 1)))

      dfw <- sdv^4/(sigma1^4*x^2/((n1 - 1)^3*n1^2) + sigma2^4*y^2/((n2 - 1)^3*n2^2))

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

      return(thres - sdv)
    }

    ## this is a leaner inplementation of uniroot in base R
    uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.0005, ...)
    {
      f <- function(x) fun(x, ...)
      x1 <- lower
      if (!is.null(f_lower)){
        f1 <- f_lower
      }
      else{
        f1 <- f(x1)
      }
      if (f1 > 0){return(x1)}
      x2 <- upper
      if (!is.null(f_upper)){
        f2 <- f_upper
      }
      else{
        f2 <- f(x2)
      }
      f2 <- f(x2)
      if (f2 < 0){return(x2)}
      x3 <- 0.5 * (lower + upper)
      niter <- 1
      while (niter <= maxiter) {
        f3 <- f(x3)
        if (abs(f3) < tol) {
          x0 <- x3
          return(x0)
        }
        if (f1 * f3 < 0) {
          upper <- x3}
        else {lower <- x3}
        if ((upper - lower) < tol2 * max(abs(upper), 1)) {
          x0 <- 0.5 * (lower + upper)
          return(x0)
        }
        denom <- (f2 - f1) * (f3 - f1) * (f2 - f3)
        numer <- x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 *
          (f2 - f3) + f1 * x2 * (f3 - f1)
        if (denom == 0) {
          dx <- upper - lower
        }
        else {
          dx <- f3 * numer/denom
        }
        x <- x3 + dx
        if ((upper - x) * (x - lower) < 0) {
          dx <- 0.5 * (upper - lower)
          x <- lower + dx
        }
        if (x1 < x3) {
          x2 <- x3
          f2 <- f3
        }
        else {
          x1 <- x3
          f1 <- f3
        }
        niter <- niter + 1
        if (abs(x - x3) < tol2) {
          x0 <- x
          return(x0)
        }
        x3 <- x
      }
      return(x0)
    }

    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    ## generate the Sobol' sequence
    sob <- qrng::sobol(2^(sobol + 10), d = 3, randomize = "digital.shift", seed = seed)

    ## find starting point for root-finding algorithm using 
    if (!is.finite(deltaU)){
      mid_val <- ((qnorm(targetPower) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(diff - deltaL))^2
      upper_val <- ((qnorm((3 +targetPower)/4) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(diff - deltaL))^2
    }
    else if (!is.finite(deltaL)){
      mid_val <- ((qnorm(targetPower) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - diff))^2
      upper_val <- ((qnorm((3 +targetPower)/4) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - diff))^2
    }
    else{
      ## find starting point for root-finding algorithm using normal approximations to
      ## the t-distribution
      a_cons <- (deltaU - diff)/sqrt(sigma1^2 + sigma2^2/q)
      b_cons <- (deltaL - diff)/sqrt(sigma1^2 + sigma2^2/q)
      c_cons <- qnorm(1-alpha)
      ## lower bound for root-finding algorithm
      lower_cons <- 2*c_cons*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - deltaL)
      upper_cons <- lower_cons
      
      fn_ci = function(n_sq, a, b, c, pwr){
        return(pnorm(a*n_sq - c) - pnorm(b*n_sq + c) - pwr)}
      
      upper_large <- FALSE
      while(upper_large == FALSE){
        upper_cons <- 10*upper_cons
        upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = targetPower)
        if (upper_check > 0){
          upper_large <- TRUE
        }
      }
      
      ## mid_val should be close to the final sample size
      mid_val <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = targetPower, 
                      lower = lower_cons, upper = upper_cons))^2
      
      upper_large <- FALSE
      upper_cons <- mid_val
      while(upper_large == FALSE){
        upper_cons <- 10*upper_cons
        upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = (3 + targetPower)/4)
        if (upper_check > 0){
          upper_large <- TRUE
        }
      }
      
      upper_val <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = (3 + targetPower)/4, 
                     lower = sqrt(mid_val), upper = sqrt(upper_cons)))^2
    }
    
    ## upper_val and lower_val will be the next sample sizes explored by the root-finding algorithm
    ## depending on whether the point corresponds to the rejection region at mid_val
    mid_val <- max(mid_val, 10)
    lower_val <- 0.5*mid_val
    
    print(c(ceiling(lower_val), ceiling(mid_val), ceiling(upper_val)))
    
    endpoints_vec <- rep(0, nrow(sob))
    samps <- NULL
    for (i in 1:nrow(sob)){
      power1 <- targetPowerfn(n_val = ceiling(mid_val),
                              deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                              u = sob[i,], alpha = alpha, q = q)
      if (power1 >= 0){
        power1b <- targetPowerfn(n_val = ceiling(lower_val),
                                 deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                 u = sob[i,], alpha = alpha, q = q)
        if (power1b >= 0){
          samps[i] <- uu(targetPowerfn, lower =2,
                             upper = ceiling(lower_val), f_upper = power1b,
                             deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                             u = sob[i,], alpha = alpha, q = q)
        }
        else{
          samps[i] <- uu(targetPowerfn, lower = ceiling(lower_val), f_lower = power1b,
                             upper = ceiling(mid_val), f_upper = power1,
                             deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                             u = sob[i,], alpha = alpha, q = q)
        }
      }
      else{
        power2 <- targetPowerfn(n_val = ceiling(upper_val),
                                deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                u = sob[i,], alpha = alpha, q = q)
        if (power2 >= 0){
          samps[i] <- uu(targetPowerfn, lower = ceiling(mid_val), f_lower = power1,
                             upper = ceiling(upper_val), f_upper = power2,
                             deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                             u = sob[i,], alpha = alpha, q = q)
        }
        else{
          endpoints_vec[i] <- 1 
        }
      }
    }
    
    last_group <- which(endpoints_vec == 1)
    if (length(last_group) == 0){
      upper_c <- 2
    } else{
      upper_c <- 1
      while(length(last_group) > 0){
        if (upper_c > 32){
          last_group <- NULL
        }
        upper_c <- 2*upper_c
        endpoints1_vec <- NULL
        for (i in 1:length(last_group)){
          endpoints1_vec[i] <- targetPowerfn(n_val = ceiling(upper_c*upper_val),
                                             deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                             u = sob[last_group[i],], alpha = alpha, q = q)
        }
        keep_vec <- ifelse(endpoints1_vec >= 0, FALSE, TRUE)
        ## only keep points that still do not satisfy power criterion after increasing
        ## the upper bound for the sample size
        last_group <- last_group[keep_vec]
      }
    }
    
    ## implement the root-finding algorithm for each point in the Sobol' sequence
    ## that required a large upper bound (i.e., those in last_group)
    for (i in 1:nrow(sob)){
      if (endpoints_vec[i] == 1){
        samps[i] <- uu(targetPowerfn, lower = ceiling(upper_val),
                       upper = ceiling(upper_c*upper_val),
                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                       u = sob[i,], alpha = alpha, q = q)
      }
    }
    
    n_rough <- quantile(samps, targetPower)

    ## explore sampling distribution at the preliminary sample size recommendation
    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)

    n1_temp <- n1; n2_temp <- n2
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))

    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))

    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))

    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }

    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }

    thres <- pmin(thresUp, thresLow)

    pwrs <- thres - sdv
    
    ## ensure that points correspond to the rejection region if their intersection
    ## is less than the recommended sample size
    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))
    
    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    samps_pre <- NULL
    samps_post <- NULL
    if (length(inconsistent) > 0){
      for (i in 1:length(inconsistent)){
        samps_pre[i] <- samps[inconsistent[i]]
        if (pwrs[inconsistent[i]] < 0){
          samps[inconsistent[i]] <- uu(targetPowerfn, lower = n_rough, upper = ceiling(upper_c*upper_val),
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfn, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        samps_post[i] <- samps[inconsistent[i]]
      }
    }
    funecdf <- stats::ecdf(samps)

    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)

    df_samps <- data.frame(n_plot = samps)

    ## plot power curve if desired
    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size=20,face="bold",
                                                          margin= ggplot2::margin(0,0,10,0))) +
        ggplot2::theme(axis.title.x = ggplot2::element_text(size = 16, margin= ggplot2::margin(10,0,0,0))) +
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = 16, margin= ggplot2::margin(0,10,0,0))) +
        ggplot2::geom_segment(ggplot2::aes(x = n_rough, y = 0, xend = n_rough, yend = targetPower), linetype="dashed", color = "black")

      plot_min <- ggplot2::ggplot_build(plot_pwr)$layout$panel_params[[1]]$x$breaks[1]
      if (is.na(plot_min)){
        plot_min <- floor(ggplot2::ggplot_build(plot_pwr)$layout$panel_params[[1]]$x$continuous_range[1])
      }

      plot_pwr <- plot_pwr +
        ggplot2::geom_segment(ggplot2::aes(x = plot_min, y = targetPower, xend = n_rough, yend = targetPower), linetype="dashed", color = "black")

      print(plot_pwr)
    }

    if (sigma1 == sigma2){
      message("This function works with equal population variances, but did you mean to use DesignParallelEqual()?")
    }
    return(list(samps, n_rough, consistency, ifelse(consistency == 0, pwrs, 0)))
  }
}

tic <- Sys.time()
temp = alg2(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2,
                         alpha = 0.05, targetPower = 0.8, q = 1, plot = FALSE, seed = 1, sobol = 0)
toc <- Sys.time()
toc-tic

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## define parameters for illustrative example
delta <- -4; sigma_1 <- 18; sigma_2 <- 15 
deltaL <- -19.2; deltaU <- 19.2; alpha = 0.05; targetP <- 0.8

## generate power curves for the left plot of Figure 2
## also used to assess whether or not the root-finding algorithm needed to be reimplemented
sob_80 <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                    .options.snow=opts, .errorhandling = "remove") %dopar% {
                      unlist(alg2(diff = delta, sigma1 = sigma_1, sigma2 = sigma_2, deltaL = deltaL, deltaU = deltaU,
                          alpha = alpha, targetPower = targetP, q = 1, plot = FALSE, seed = i + 200, sobol = 0))
                    }

write.csv(sob_80[, 1:1024], "samps_SOB_1024.csv", row.names = FALSE)
write.csv(sob_80[, 1026:2049], "consistency_80.csv", row.names = FALSE)

## assess whether or not the root-finding algorithm needed to be reimplemented for the other values of power
targetPs <- c(20, 30, 40, 50, 60, 70, 90)
for (k in 1:length(targetPs)){
  tempP <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                    .options.snow=opts, .errorhandling = "remove") %dopar% {
                      unlist(alg2(diff = delta, sigma1 = sigma_1, sigma2 = sigma_2, deltaL = deltaL, deltaU = deltaU,
                                  alpha = alpha, targetPower = targetPs[k]/100, q = 1, plot = FALSE, seed = 200 + k*1000 + i, sobol = 0))
                    }
  print(c(targetPs[k], sum(tempP[, 1026:2049] == 0)))
  write.csv(tempP[, 1026:2049], paste0("consistency_", targetPs[k],".csv"), row.names = FALSE)
  # write.csv(tempP[, 2050:3073], paste0("pwrs_", targetPs[k],".csv"), row.names = FALSE)
}

## this function provides an implementation of algorithm 2 with pseudorandom sequences
alg2PRNG <- function(diff = NULL, sigma1 = NULL, sigma2 = NULL, deltaL = -Inf,
                 deltaU = Inf, alpha = NULL, targetPower = NULL, q = 1, 
                 plot = TRUE, seed = NULL, m = 1024, check = FALSE){
  
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  if (check == TRUE){
    ## error checking
    if(!is.numeric(diff) | length(diff) != 1) {
      stop("Please specify a valid number for diff.")}
    if(!is.numeric(sigma1) | length(sigma1) != 1){
      stop("Please specify a valid number for sigma1.")}
    else if (sigma1 <= 0 | !is.finite(sigma1)){
      stop("Please specify a valid number for sigma1.")}
    if(!is.numeric(sigma2) | length(sigma2) != 1){
      stop("Please specify a valid number for sigma2.")}
    else if (sigma2 <= 0 | !is.finite(sigma2)){
      stop("Please specify a valid number for sigma2.")}
    if(!is.numeric(deltaL) | length(deltaL) != 1){
      stop("Please specify a valid number for deltaL.")}
    if(!is.numeric(deltaU) | length(deltaU) != 1){
      stop("Please specify a valid number for deltaU.")}
    if(deltaL == -Inf & deltaU == Inf){
      stop("Please specify valid interval endpoints deltaL and deltaU.")}
    if (deltaL >= deltaU){
      stop("Please specify valid interval endpoints deltaL and deltaU.")}
    if(!is.numeric(alpha) | length(alpha) != 1) {
      stop("Please specify a valid number for alpha.")}
    if (is.numeric(alpha)){
      if (alpha <= 0 | alpha >= 1){
        stop("Please specify a valid number for alpha.")}
    }
    if(!is.numeric(targetPower)) {
      stop("Please specify a valid number for targetPower.")}
    if (is.numeric(targetPower)){
      if (targetPower <= 0 | targetPower >= 1){
        stop("Please specify a valid number for targetPower.")}
      if(!is.numeric(q)) {
        stop("Please specify a valid number for q.")}
      else if (is.numeric(q)){
        if (q <= 0) {
          stop("Please specify a valid number for q.")}
      }
      if (diff >= deltaU | diff <= deltaL){
        stop("Please ensure diff is between deltaL and deltaU.")
      }
    }
    if(!is.null(seed) & (!is.numeric(seed) | length(seed) != 1)) {
      stop("Please specify a valid seed for random number generation.")}
    else if (!is.null(seed)){
      if (seed%%1 != 0 | seed < 0){
        stop("Please specify a valid seed for random number generation.")}
    }
    if (length(plot) != 1 | !(plot %in% c(TRUE, FALSE))){
      stop("Please provide valid logical input for plot.")
    }
  }
  
  if (length(targetPower) == 1 & length(q) == 1){
    ## this function returns a positive result if the sample corresponds to the rejection region
    targetPowerfn <- function(u, deltaL, deltaU, diff, sigma1, sigma2, n_val, alpha, q){
      n1 <- n_val[1]; n2 <- max(2, q*n_val)
      x <- stats::qchisq(u[1], n1 - 1)
      y <- stats::qchisq(u[2], n2 - 1)
      z <- stats::qnorm(u[3], diff, sqrt(sigma1^2/n1 + sigma2^2/n2))
      
      sdv <- sqrt(sigma1^2*x/((n1 - 1)*n1) + sigma2^2*y/(n2*(n2 - 1)))
      
      dfw <- sdv^4/(sigma1^4*x^2/((n1 - 1)^3*n1^2) + sigma2^4*y^2/((n2 - 1)^3*n2^2))
      
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
      
      return(thres - sdv)
    }
    
    ## this is a leaner inplementation of uniroot in base R
    uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.0005, ...)
    {
      f <- function(x) fun(x, ...)
      x1 <- lower
      if (!is.null(f_lower)){
        f1 <- f_lower
      }
      else{
        f1 <- f(x1)
      }
      if (f1 > 0){return(x1)}
      x2 <- upper
      if (!is.null(f_upper)){
        f2 <- f_upper
      }
      else{
        f2 <- f(x2)
      }
      f2 <- f(x2)
      if (f2 < 0){return(x2)}
      x3 <- 0.5 * (lower + upper)
      niter <- 1
      while (niter <= maxiter) {
        f3 <- f(x3)
        if (abs(f3) < tol) {
          x0 <- x3
          return(x0)
        }
        if (f1 * f3 < 0) {
          upper <- x3}
        else {lower <- x3}
        if ((upper - lower) < tol2 * max(abs(upper), 1)) {
          x0 <- 0.5 * (lower + upper)
          return(x0)
        }
        denom <- (f2 - f1) * (f3 - f1) * (f2 - f3)
        numer <- x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 *
          (f2 - f3) + f1 * x2 * (f3 - f1)
        if (denom == 0) {
          dx <- upper - lower
        }
        else {
          dx <- f3 * numer/denom
        }
        x <- x3 + dx
        if ((upper - x) * (x - lower) < 0) {
          dx <- 0.5 * (upper - lower)
          x <- lower + dx
        }
        if (x1 < x3) {
          x2 <- x3
          f2 <- f3
        }
        else {
          x1 <- x3
          f1 <- f3
        }
        niter <- niter + 1
        if (abs(x - x3) < tol2) {
          x0 <- x
          return(x0)
        }
        x3 <- x
      }
      return(x0)
    }
    
    if (is.null(seed)){
      seed <- ceiling(1000*stats::runif(1))
    }
    ## generate the pseudrandom sequence
    # sob <- qrng::sobol(2^(sobol + 10), d = 3, randomize = "digital.shift", seed = seed)
    sob <- matrix(runif(3*m), ncol = 3)
    
    ## find starting point for root-finding algorithm using 
    if (!is.finite(deltaU)){
      mid_val <- ((qnorm(targetPower) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(diff - deltaL))^2
      upper_val <- ((qnorm((3 +targetPower)/4) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(diff - deltaL))^2
    }
    else if (!is.finite(deltaL)){
      mid_val <- ((qnorm(targetPower) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - diff))^2
      upper_val <- ((qnorm((3 +targetPower)/4) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - diff))^2
    }
    else{
      ## find starting point for root-finding algorithm using normal approximations to
      ## the t-distribution
      a_cons <- (deltaU - diff)/sqrt(sigma1^2 + sigma2^2/q)
      b_cons <- (deltaL - diff)/sqrt(sigma1^2 + sigma2^2/q)
      c_cons <- qnorm(1-alpha)
      ## lower bound for root-finding algorithm
      lower_cons <- 2*c_cons*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - deltaL)
      upper_cons <- lower_cons
      
      fn_ci = function(n_sq, a, b, c, pwr){
        return(pnorm(a*n_sq - c) - pnorm(b*n_sq + c) - pwr)}
      
      upper_large <- FALSE
      while(upper_large == FALSE){
        upper_cons <- 10*upper_cons
        upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = targetPower)
        if (upper_check > 0){
          upper_large <- TRUE
        }
      }
      
      ## mid_val should be close to the final sample size
      mid_val <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = targetPower, 
                     lower = lower_cons, upper = upper_cons))^2
      
      upper_large <- FALSE
      upper_cons <- mid_val
      while(upper_large == FALSE){
        upper_cons <- 10*upper_cons
        upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = (3 + targetPower)/4)
        if (upper_check > 0){
          upper_large <- TRUE
        }
      }
      
      upper_val <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = (3 + targetPower)/4, 
                       lower = sqrt(mid_val), upper = sqrt(upper_cons)))^2
    }
    
    ## upper_val and lower_val will be the next sample sizes explored by the root-finding algorithm
    ## depending on whether the point corresponds to the rejection region at mid_val
    mid_val <- max(mid_val, 10)
    lower_val <- 0.5*mid_val
    
    print(c(ceiling(lower_val), ceiling(mid_val), ceiling(upper_val)))
    
    endpoints_vec <- rep(0, nrow(sob))
    samps <- NULL
    for (i in 1:nrow(sob)){
      power1 <- targetPowerfn(n_val = ceiling(mid_val),
                              deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                              u = sob[i,], alpha = alpha, q = q)
      if (power1 >= 0){
        power1b <- targetPowerfn(n_val = ceiling(lower_val),
                                 deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                 u = sob[i,], alpha = alpha, q = q)
        if (power1b >= 0){
          samps[i] <- uu(targetPowerfn, lower =2,
                         upper = ceiling(lower_val), f_upper = power1b,
                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                         u = sob[i,], alpha = alpha, q = q)
        }
        else{
          samps[i] <- uu(targetPowerfn, lower = ceiling(lower_val), f_lower = power1b,
                         upper = ceiling(mid_val), f_upper = power1,
                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                         u = sob[i,], alpha = alpha, q = q)
        }
      }
      else{
        power2 <- targetPowerfn(n_val = ceiling(upper_val),
                                deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                u = sob[i,], alpha = alpha, q = q)
        if (power2 >= 0){
          samps[i] <- uu(targetPowerfn, lower = ceiling(mid_val), f_lower = power1,
                         upper = ceiling(upper_val), f_upper = power2,
                         deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                         u = sob[i,], alpha = alpha, q = q)
        }
        else{
          endpoints_vec[i] <- 1 
        }
      }
    }
    
    last_group <- which(endpoints_vec == 1)
    if (length(last_group) == 0){
      upper_c <- 2
    } else{
      upper_c <- 1
      while(length(last_group) > 0){
        if (upper_c > 32){
          last_group <- NULL
        }
        upper_c <- 2*upper_c
        endpoints1_vec <- NULL
        for (i in 1:length(last_group)){
          endpoints1_vec[i] <- targetPowerfn(n_val = ceiling(upper_c*upper_val),
                                             deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                             u = sob[last_group[i],], alpha = alpha, q = q)
        }
        keep_vec <- ifelse(endpoints1_vec >= 0, FALSE, TRUE)
        ## only keep points that still do not satisfy power criterion after increasing
        ## the upper bound for the sample size
        last_group <- last_group[keep_vec]
      }
    }
    
    ## implement the root-finding algorithm for each point in the Sobol' sequence
    ## that required a large upper bound (i.e., those in last_group)
    for (i in 1:nrow(sob)){
      if (endpoints_vec[i] == 1){
        samps[i] <- uu(targetPowerfn, lower = ceiling(upper_val),
                       upper = ceiling(upper_c*upper_val),
                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                       u = sob[i,], alpha = alpha, q = q)
      }
    }
    
    n_rough <- quantile(samps, targetPower)
    
    ## explore sampling distribution at the preliminary sample size recommendation
    pwrs <- NULL
    n1 <- n_rough
    n2 <- pmax(2, q*n1)
    
    n1_temp <- n1; n2_temp <- n2
    ## generate two sds and one mean
    x <- stats::qchisq(sob[,1], n1_temp - 1)
    y <- stats::qchisq(sob[,2], n2_temp - 1)
    z <- stats::qnorm(sob[,3], diff, sqrt(sigma1^2/n1_temp + sigma2^2/n2_temp))
    
    sdv <- sqrt(sigma1^2*x/((n1_temp - 1)*n1_temp) + sigma2^2*y/(n2_temp*(n2_temp - 1)))
    
    dfw <- sdv^4/(sigma1^4*x^2/((n1_temp - 1)^3*n1_temp^2) + sigma2^4*y^2/((n2_temp - 1)^3*n2_temp^2))
    
    if (deltaU == Inf){
      thresUp <- rep(Inf, length(z))
    }
    else{
      thresUp <- (deltaU - z)/stats::qt(1-alpha, dfw)
    }
    
    if (deltaL == -Inf){
      thresLow <- rep(Inf, length(z))
    }
    else{
      thresLow <- (z - deltaL)/stats::qt(1-alpha, dfw)
    }
    
    thres <- pmin(thresUp, thresLow)
    
    pwrs <- thres - sdv
    
    ## ensure that points correspond to the rejection region if their intersection
    ## is less than the recommended sample size
    consistency <- ifelse(samps < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.02))
    
    ## for any points where the root-finding algorithm has caused issues,
    ## re-run the root-finding algorithm starting at n_*
    inconsistent <- which(consistency == 0)
    samps_pre <- NULL
    samps_post <- NULL
    if (length(inconsistent) > 0){
      for (i in 1:length(inconsistent)){
        samps_pre[i] <- samps[inconsistent[i]]
        if (pwrs[inconsistent[i]] < 0){
          samps[inconsistent[i]] <- uu(targetPowerfn, lower = n_rough, upper = ceiling(upper_c*upper_val),
                                       f_lower = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        else{
          samps[inconsistent[i]] <- uu(targetPowerfn, lower = 2, upper = n_rough,
                                       f_upper = pwrs[inconsistent[i]],
                                       deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                       u = sob[inconsistent[i],], alpha = alpha, q = q)
        }
        samps_post[i] <- samps[inconsistent[i]]
      }
    }
    funecdf <- stats::ecdf(samps)
    
    ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
    n_rough <- uu(ecdf_root, lower = stats::quantile(samps, targetPower*0.5), upper = stats::quantile(samps, targetPower + 0.5*(1 - targetPower)),
                  pwr = targetPower)
    
    df_samps <- data.frame(n_plot = samps)
    
    ## plot power curve if desired
    n_plot <- NULL
    if (plot == TRUE){
      plot_pwr <- ggplot2::ggplot(df_samps, ggplot2::aes(x = n_plot)) + ggplot2::theme_bw() +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, colour = cbbPalette[6], linewidth = 2) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13)) +
        ggplot2::theme(axis.text.x =  ggplot2::element_text(size = 13)) +
        ggplot2::labs(title = "Approximated Power Curve") +
        ggplot2::labs(y = "Power", x=bquote(italic(n)[1]*"  ("*italic(n)[2]*" = "*.(round(q,3))*italic(n)[1]*")")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size=20,face="bold",
                                                          margin= ggplot2::margin(0,0,10,0))) +
        ggplot2::theme(axis.title.x = ggplot2::element_text(size = 16, margin= ggplot2::margin(10,0,0,0))) +
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = 16, margin= ggplot2::margin(0,10,0,0))) +
        ggplot2::geom_segment(ggplot2::aes(x = n_rough, y = 0, xend = n_rough, yend = targetPower), linetype="dashed", color = "black")
      
      plot_min <- ggplot2::ggplot_build(plot_pwr)$layout$panel_params[[1]]$x$breaks[1]
      if (is.na(plot_min)){
        plot_min <- floor(ggplot2::ggplot_build(plot_pwr)$layout$panel_params[[1]]$x$continuous_range[1])
      }
      
      plot_pwr <- plot_pwr +
        ggplot2::geom_segment(ggplot2::aes(x = plot_min, y = targetPower, xend = n_rough, yend = targetPower), linetype="dashed", color = "black")
      
      print(plot_pwr)
    }
    
    if (sigma1 == sigma2){
      message("This function works with equal population variances, but did you mean to use DesignParallelEqual()?")
    }
    return(list(samps, n_rough, consistency))
  }
}

## define parameters for illustrative example
delta <- -4; sigma_1 <- 18; sigma_2 <- 15 
deltaL <- -19.2; deltaU <- 19.2; alpha = 0.05; targetP <- 0.8

## get 1000 power curves to create bootstrap confidence intervals with pseudorandom sequences (m = 1024)
prng_1024 <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                  .options.snow=opts, .errorhandling = "remove") %dopar% {
                    unlist(alg2PRNG(diff = delta, sigma1 = sigma_1, sigma2 = sigma_2, deltaL = deltaL, deltaU = deltaU,
                                alpha = alpha, targetPower = targetP, q = 1, plot = FALSE, seed = i + 8200, m = 1024))
                  }

write.csv(prng_1024[, 1:1024], "samps_PRNG_1024.csv", row.names = FALSE)

## get 1000 power curves to create bootstrap confidence intervals with pseudorandom sequences (m = 10000)
prng_10k <- foreach(i=1:1000, .combine='rbind', .packages = c("qrng"),
                     .options.snow=opts, .errorhandling = "remove") %dopar% {
                       unlist(alg2PRNG(diff = delta, sigma1 = sigma_1, sigma2 = sigma_2, deltaL = deltaL, deltaU = deltaU,
                                       alpha = alpha, targetPower = targetP, q = 1, plot = FALSE, seed = i + 9200, m = 10000))
                     }

write.csv(prng_10k[, 1:10000], "samps_PRNG_10k.csv", row.names = FALSE)

## now compute the bootstrap confidence intervals by taking the
## 2.5th percentile and 97.5th percentile of power estimates
data_sim_x <- c(3,5,8,10,15,20,30,40,50,60)
upperSOB <- NULL
lowerSOB <- NULL
alpha <- 0.025
for (i in 1:length(data_sim_x)){
  powers <- rowMeans(sob_80[, 1:1024] <= data_sim_x[i])
  upperSOB <- c(upperSOB, quantile(powers, 1 - alpha))
  lowerSOB <- c(lowerSOB, quantile(powers, alpha))
}

write.csv(upperSOB, "upperSOB.csv", row.names = FALSE)
write.csv(lowerSOB, "lowerSOB.csv", row.names = FALSE)

## repeat for PRNG (m = 1024)
data_sim_x <- c(3,5,8,10,15,20,30,40,50,60)
upperPRNG <- NULL
lowerPRNG <- NULL
alpha <- 0.025
for (i in 1:length(data_sim_x)){
  powers <- rowMeans(prng_1024[, 1:1024] <= data_sim_x[i])
  upperPRNG <- c(upperPRNG, quantile(powers, 1 - alpha))
  lowerPRNG <- c(lowerPRNG, quantile(powers, alpha))
}

write.csv(upperPRNG, "upperPRNG.csv", row.names = FALSE)
write.csv(lowerPRNG, "lowerPRNG.csv", row.names = FALSE)

## repeat for PRNG (m = 10000)
data_sim_x <- c(3,5,8,10,15,20,30,40,50,60)
upperPRNG10k <- NULL
lowerPRNG10k <- NULL
alpha <- 0.025
for (i in 1:length(data_sim_x)){
  powers <- rowMeans(prng_10k[, 1:10000] <= data_sim_x[i])
  upperPRNG10k <- c(upperPRNG10k, quantile(powers, 1 - alpha))
  lowerPRNG10k <- c(lowerPRNG10k, quantile(powers, alpha))
}

write.csv(upperPRNG10k, "upperPRNG_10k.csv", row.names = FALSE)
write.csv(lowerPRNG10k, "lowerPRNG_10k.csv", row.names = FALSE)

## generate figure 2

## start with left plot
## the PDF file of the plot becomes too large when all 1000 curves are imposed on plot.
## We visualize the range between the minimum and minimum power estimates for a variety of sample sizes.
fun_lt <- function(x, cutoff){mean(x <= cutoff)}
samps_mat_1024 <- read.csv("samps_SOB_1024.csv")
samps_cut <- seq(2, 65, 0.25)
res_min <- NULL
res_max <- NULL
for (i in 1:length(samps_cut)){
  temp <- apply(samps_mat_1024, 1, fun_lt, cutoff = samps_cut[i])
  res_min <- c(res_min, min(temp))
  res_max <- c(res_max, max(temp))
}
data_full <- data.frame(n = c(samps_cut, rev(samps_cut)), power = c(res_min, rev(res_max)),
                        curve = c(rep(1, length(samps_cut)), rep(2, length(samps_cut))))

## get power estimates from Alg 1
data_sim_y <- c(0.0414, 0.1283, 0.3801, 0.5366, 0.7699, 0.8815, 0.9687, 0.9923, 0.9982, 0.9996)
data_sim <- data.frame(n = data_sim_x, power = data_sim_y)

## combine polygon representing 1000 power curves with the unbiased power estimates
## from Algorithm 1
assign(paste0("plot2"), ggplot(get(paste0("data_full")), aes(x=n)) + theme_bw() +
         geom_polygon(aes(x = n, y = power), alpha = 0.25) +
         labs(x= bquote(italic(n)), y= bquote('Power')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.position="none") +
         scale_color_manual(name = " ",
                            values = c("firebrick", "grey")) +
         theme(legend.text=element_text(size=18)) +
         ylim(0,1) + 
         xlim(0, 65) +
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + 
         geom_point(data = get(paste0("data_sim")), aes(x = n, y = power, color="firebrick"), size =2)
)
plot2

## create right plot for confidence interval lines
upperSOB <- read.csv("upperSOB.csv")$x
upperPRNG <- read.csv("upperPRNG.csv")$x
upperPRNG10k <- read.csv("upperPRNG_10k.csv")$x
df_upper <- data.frame(Sequence = c(rep("1Sobol'",10), rep("2PRNG",10), rep("3PRNG10k",10)), 
                       Difference = c(upperSOB - data_sim_y, upperPRNG - data_sim_y, upperPRNG10k - data_sim_y), 
                       n = rep(c(3,5,8,10,15,20,30,40,50,60),3))

lowerSOB <- read.csv("lowerSOB.csv")$x
lowerPRNG <- read.csv("lowerPRNG.csv")$x
lowerPRNG10k <- read.csv("lowerPRNG_10k.csv")$x
df_lower <- data.frame(Sequence = c(rep("4Sobol'",10), rep("5PRNG",10), rep("6PRNG10k",10)), 
                       Difference = c(lowerSOB - data_sim_y, lowerPRNG - data_sim_y, lowerPRNG10k - data_sim_y),
                       n = rep(c(3,5,8,10,15,20,30,40,50,60),3))
df_full <- rbind(df_upper, df_lower)

## use this first plot to get the legend (otherwise each label shows up twice for upper and lower limit)
plot_temp <- ggplot(df_upper, aes(x=n, y=Difference, color = as.factor(Sequence), linetype = as.factor(Sequence))) + theme_bw() +
  geom_line() + ylim(-0.04, 0.04) +
  geom_point() +
  scale_color_manual(name = "", labels = rep(c("Sobol'  ", "PRNG  ", "Long PRNG"),1),
                     values = rep(c("firebrick", "steelblue", "gold4"),1)) +
  scale_linetype_manual(name = "", labels = rep(c("Sobol'  ", "PRNG  ", "Long PRNG"),1),
                        values = rep(c("solid", "dashed", "dotted"),1)) +
  labs(color  = "", linetype = "") +
  labs(x= bquote(italic(n)), y= bquote('Centered 95% Confidence Interval')) +
  theme(plot.title = element_text(size=20,face="bold",
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.text=element_text(size=18)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.position="bottom") +
  theme(legend.key.size = unit(0.85, "cm"))

## extract legend
mylegend<-get_legend(plot_temp)

## combine the following plot with the previous legend
plot3a <- ggplot(df_full, aes(x=n, y=Difference, color = as.factor(Sequence), linetype = as.factor(Sequence))) + theme_bw() +
  geom_line() + ylim(-0.04, 0.04) +
  geom_point() +
  scale_color_manual(name = "", labels = rep(c("Sobol'  ", "PRNG  ", "Long PRNG"),2),
                     values = rep(c("firebrick", "steelblue", "gold4"),2)) +
  scale_linetype_manual(name = "", labels = rep(c("Sobol'  ", "PRNG  ", "Long PRNG"),2),
                        values = rep(c("solid", "dashed", "dotted"),2)) +
  labs(color  = "", linetype = "") +
  labs(x= bquote(italic(n)), y= bquote('Centered 95% Confidence Interval')) +
  theme(plot.title = element_text(size=20,face="bold",
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.text=element_text(size=18)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.position="none")

plot3 <- plot_grid(plot3a, mylegend, nrow=2, rel_heights = c(10,1))

SimCurve <- plot_grid(plot2 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                      plot3 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), rel_widths = c(1,1))

# output as .pdf file for the article
pdf(file = "Figure2.pdf",   # The directory you want to save the file in
    width = 11.025, # The width of the plot in inches (12.41)
    height = 4.9875) # The height of the plot in inches (10.7)

SimCurve

dev.off()
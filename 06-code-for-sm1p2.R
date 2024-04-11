## code to reproduce numerical study in Section SM1.2 of the supplement

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

## we require algorithm 2 again for this study (code is taken from file 03)
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
    uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.0000001, ...)
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

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 100, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## define inputs for 280 simulation settings
sig1_vec <- c(16.5, rep(18,3), rep(19.5,3))
sig2_vec <- c(16.5, rep(15,3), rep(13,3))
q_vec <- c(1, 1, 15/18, 18/15, 1, 13/19.5, 19.5/13)
q_title <- c("1", rep(c("1", "opt", "inv"), 2))
q_num <- c("1", rep("2", 3), rep("3", 3))

targetPs <- seq(0.2, 0.9, 0.1)
mu_temps <- c(0, -4, -8, -12, -16)
## track2 tracks the seed for each setting for further investigation
## problem_set tracks the settings where the root-finding algorithm needed to be 
## reimplemented (and how many times it needed to be reimplemented)
track_seed <- 0
problem_set <- NULL
track2 <- NULL
for (jj in 1:length(mu_temps)){
  for (ll in 1:length(sig1_vec)){
    for (kk in 1:length(targetPs)){
      delta <- mu_temps[jj]; sigma_1 <- sig1_vec[ll]; sigma_2 <- sig2_vec[ll]; q <- q_vec[ll] 
      deltaL <- -19.2; deltaU <- 19.2; alpha = 0.05; targetP <- targetPs[kk]
      
      track2 <- rbind(track2, c(jj, ll, kk, track_seed))
      temp <- foreach(i=1:100, .combine='rbind', .packages = c("qrng"),
                       .options.snow=opts, .errorhandling = "remove") %dopar% {
                         unlist(alg2(diff = delta, sigma1 = sigma_1, sigma2 = sigma_2, deltaL = deltaL, deltaU = deltaU,
                                     alpha = alpha, targetPower = targetP, q = q, plot = FALSE, seed = track_seed + i, sobol = 0))
                       }
      track_seed <- track_seed + 100
      ## the entries in this matrix are 0 if we need to reimplement the root-finding algorithm
      if (sum(temp[, 1026:2049] == 0) > 0){
        problem_set <- rbind(problem_set, c(mu_temps[jj], ll, targetPs[kk], sum(temp[, 1026:2049] == 0)))
      }
      print(c(mu_temps[jj], ll, targetPs[kk], sum(temp[, 1026:2049] == 0)))
      write.csv(temp[, 1026:2049], paste0("consistency_mu_",-1*delta,"_set_", ll, "_pwr_", 100*targetP,".csv"), row.names = FALSE)
    }
  }
}
write.csv(track2, "track2.csv", row.names = FALSE)
write.csv(problem_set, "problem_set.csv", row.names = FALSE)
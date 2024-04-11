## code to reproduce figure 3

## BEGIN SETUP ##
require(qrng)
require(foreach)
require(doParallel)
require(doSNOW)
require(ggplot2)
require(cowplot)
require(ggpubr)
require(imager)
require(magick)

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

## this version of the algorithm returns the sample sizes that were explored by the root-finding algorithm
## (i.e., the intermediate steps, not just the final ones)
alg2Root <- function(diff = NULL, sigma1 = NULL, sigma2 = NULL, deltaL = -Inf,
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
    
    uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.005, ...)
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
    sob <- qrng::sobol(2^(sobol + 10), d = 3, randomize = "digital.shift", seed = seed)
    
    ## find starting point for root-finding algorithm
    if (!is.finite(deltaU)){
      mid_val <- ((qnorm(targetPower) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(diff - deltaL))^2
      upper_val <- ((qnorm((3 +targetPower)/4) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(diff - deltaL))^2
    }
    else if (!is.finite(deltaL)){
      mid_val <- ((qnorm(targetPower) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - diff))^2
      upper_val <- ((qnorm((3 +targetPower)/4) + qnorm(1 - alpha))*sqrt(sigma1^2 + sigma2^2/q)/(deltaU - diff))^2
    }
    else{
      ## find more conservative upper bound using criteria for credible intervals
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
    
    mid_val <- max(mid_val, 10)
    lower_val <- 0.5*mid_val
    
    print(c(ceiling(lower_val), ceiling(mid_val), ceiling(upper_val)))
    
    uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.005, ...)
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
      save <- c(x1, x2)
      niter <- 1
      while (niter <= maxiter) {
        save <- c(save, x3)
        f3 <- f(x3)
        if (abs(f3) < tol) {
          x0 <- x3
          return(c(x0, save))
        }
        if (f1 * f3 < 0) {
          upper <- x3}
        else {lower <- x3}
        if ((upper - lower) < tol2 * max(abs(upper), 1)) {
          x0 <- 0.5 * (lower + upper)
          return(c(x0, save))
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
          return(c(x0, save))
        }
        x3 <- x
      }
      return(c(x0, save))
    }
    
    endpoints_vec <- rep(0, nrow(sob))
    samps <- matrix(0, nrow = nrow(sob), ncol = 100)
    for (i in 1:nrow(sob)){
      power1 <- targetPowerfn(n_val = ceiling(mid_val),
                              deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                              u = sob[i,], alpha = alpha, q = q)
      if (power1 >= 0){
        power1b <- targetPowerfn(n_val = ceiling(lower_val),
                                 deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                 u = sob[i,], alpha = alpha, q = q)
        if (power1b >= 0){
          temp <- uu(targetPowerfn, lower =2,
                     upper = ceiling(lower_val), f_upper = power1b,
                     deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                     u = sob[i,], alpha = alpha, q = q)
          samps[i, 1: length(temp)] <- temp
        }
        else{
          temp <- uu(targetPowerfn, lower = ceiling(lower_val), f_lower = power1b,
                     upper = ceiling(mid_val), f_upper = power1,
                     deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                     u = sob[i,], alpha = alpha, q = q)
          samps[i, 1: length(temp)] <- temp
        }
      }
      else{
        power2 <- targetPowerfn(n_val = ceiling(upper_val),
                                deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                                u = sob[i,], alpha = alpha, q = q)
        if (power2 >= 0){
          temp <- uu(targetPowerfn, lower = ceiling(mid_val), f_lower = power1,
                     upper = ceiling(upper_val), f_upper = power2,
                     deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                     u = sob[i,], alpha = alpha, q = q)
          samps[i, 1: length(temp)] <- temp
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
        temp <- uu(targetPowerfn, lower = ceiling(upper_val),
                   upper = ceiling(upper_c*upper_val),
                   deltaL = deltaL, deltaU = deltaU, diff = diff, sigma1 = sigma1, sigma2 = sigma2,
                   u = sob[i,], alpha = alpha, q = q)
        samps[i, 1: length(temp)] <- temp
      }
    }
    
    samps_final <- samps[,1]
    samps_pwr <- samps[,-1]
    
    n_rough <- quantile(samps_final, targetPower)
    
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
    
    ## check consistency
    consistency <- ifelse(samps_final < n_rough, round(pwrs,2) >= 0, round(pwrs,2) <= 0)
    consistency <- ifelse(consistency == 1, 1, as.numeric(abs(pwrs) < 0.01))
    
    return(list(samps_final, samps_pwr, n_rough, consistency, sob))
  }
}

## output the results from the root-finding algorithm for one power curve approximation with
## the illustrative example for a target power of 0.8
temp = alg2Root(diff = -4, sigma1 = 18, sigma2 = 15, deltaL = -19.2, deltaU = 19.2,
                               alpha = 0.05, targetPower = 0.8, q = 1, plot = FALSE, seed = 1, sobol = 0)

## separate the final sample sizes from the intermediary ones  
finals <- temp[[1]]
prelims <- temp[[2]]

## this code block determines whether or not a point was used to explore at least
## one sample size in a given range
how_many <- apply(prelims, 1, function(x){sum(x > 0)})
prelims <- prelims[, seq(1, range(how_many)[2])]
prelims <- cbind(rep(16, length(finals)), prelims)
  
fun <- function(x, xs){
  return(sum(ifelse(x >= xs[1], x <= xs[2], 0))>0)
}
 
## our cut points for these intervals are 2, 8, 16, and 26 
end1 <- c(2.05, 8.05, 16.05, 26.1)
end2 <- c(7.95, 15.95, 25.95, 60)
df.set <- data.frame(u1 = -1, u2 = -1, u3 = -1, type = 0, index = 0)
for (i in 1:length(end1)){
  results_temp <- apply(prelims, 1, fun, xs = c(end1[i],end2[i]))
  df.set <-  rbind(df.set,
                    data.frame(u1 = temp[[5]][results_temp, 1],
                                u2 = temp[[5]][results_temp, 2], 
                                u3 = temp[[5]][results_temp, 3], type = i,
                                index = which(results_temp)))
}
df.set <- df.set[-1,]
## add the points that corresponded to the rejection region at n = 2
twos <- which(!(1:1024 %in% df.set$index))
df.set <-  rbind(
  data.frame(u1 = temp[[5]][twos, 1],
               u2 = temp[[5]][twos, 2], 
               u3 = temp[[5]][twos, 3], type = 0,
               index =twos), df.set)
df.set$type <- df.set$type + 1

## there are five coloured categorizations; visualize all of them on the same plot  
plot3d(x = subset(df.set[,1], df.set$type == 1),
        y = subset(df.set[,2], df.set$type == 1),
        z = subset(df.set[,3], df.set$type == 1), col = "#5E4FA2", xlim = c(0,1), ylim = c(0,1), zlim = c(0,1), size =10,
        xlab = "", ylab = "", zlab = "",
        axes = FALSE)
axis3d('y--',nticks =0,labels = FALSE)
axis3d('y-+',nticks=0,labels = FALSE)
axis3d('y+-',nticks=1,labels = TRUE, cex = 2)
mtext3d(expression(italic(u)[1]), "y+-", line = 1.5, cex = 2)
axis3d('y++',nticks=0,labels = FALSE)
axis3d('x--',nticks=1,labels = TRUE, cex = 2)
mtext3d(expression(italic(u)[2]), "x--", line = 1.5, cex = 2)
axis3d('x-+',nticks=0,labels = FALSE)
axis3d('x+-',nticks=0,labels = FALSE)
axis3d('x++',nticks=0,labels = FALSE)
axis3d('z--',nticks=1,labels = TRUE, cex = 2)
mtext3d(expression(italic(u)[3]), "z--", line = 1.5, cex = 2)
axis3d('z-+',nticks=0,labels = FALSE)
axis3d('z+-',nticks=0,labels = FALSE)
axis3d('z++',nticks=0,labels = FALSE)
plot3d(x = subset(df.set[,1], df.set$type == 2),
      y = subset(df.set[,2], df.set$type == 2),
      z = subset(df.set[,3], df.set$type == 2), col = "#3288BD", size =10, add = TRUE)
plot3d(x = subset(df.set[,1], df.set$type == 3),
       y = subset(df.set[,2], df.set$type == 3),
       z = subset(df.set[,3], df.set$type == 3), col = "forestgreen", size =10, add = TRUE)
plot3d(x = subset(df.set[,1], df.set$type == 4),
       y = subset(df.set[,2], df.set$type == 4),
       z = subset(df.set[,3], df.set$type == 4), col = "#FDAE61", size =10, add = TRUE)
plot3d(x = subset(df.set[,1], df.set$type == 5),
       y = subset(df.set[,2], df.set$type == 5),
       z = subset(df.set[,3], df.set$type == 5), col = "#D53E4F", size =10, add = TRUE)

rgl.snapshot('targeted.png', fmt = 'png')

## this code is used to create a legend for the previous .png file in ggplot()
k = 1
assign(paste0("df.set", k, "a2"), data.frame( u1 = rep(0.5,5), u2 = rep(0.5, 5), type = seq(1,5,1)))
## extract legend
for (k in 1){
  assign(paste0("p", k, "2"), ggplot() + theme_bw() + geom_point(
    data = get(paste0("df.set", k, "a2"))
    ,aes(
      x=u1
      ,y=u2
      ,group=factor(type), colour = factor(type), fill = factor(type)
    ),
    alpha = 0.75, size = 2.5
  ) + xlim(0,1) + ylim(0,1) +
    labs(title=paste0("Setting ", k, "a")) +
    theme(plot.title = element_text(size=20,face="bold",
                                    margin = margin(t = 0, 0, 5, 0))) +
    theme(axis.text.x = element_text(size = 15)) +
    theme(axis.text.y = element_text(size = 15)) +
    labs(x= expression(u[1])) +
    labs(y= expression(u[2])) +
    scale_colour_manual(
      name = expression(italic(n)),
      labels = c("2", "(2, 8)", "(8, 16)", "(16, 26)", "(26, 60)"),
      values = c("#5E4FA2", "#3288BD", "forestgreen","#FDAE61","#D53E4F"),
      aesthetics = c("colour", "fill")
    ) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_text(size = 16, margin = margin(t = 5, r = 0, b = 0, l = 0))) +
    theme(axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 5, b = 0, l = 0))))
}

## combine the png file and legend
p12_legend <- p12 + theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=18)) +
  theme(legend.title=element_text(size=18)) +
  guides(col=guide_legend(nrow=2,byrow=TRUE))

img1 <- ggplot()
image <- load.image(paste0("targeted.png"))

img1 <- ggdraw(img1) + draw_image(image)

fig_left <- plot_grid(plot_grid(NULL, img1, NULL, ncol =3, rel_widths = c(0.55, 4.65, 0.55)), 
                 get_legend(p12_legend), nrow = 2, rel_heights = c(10, 2.5))

## now get the p-values at n = 8 for right plot
delta <- -4; sigma_1 <- 18; sigma_2 <- 15 
deltaL <- -19.2; deltaU <- 19.2; alpha = 0.05
n1 <- 8; n2 <- 8

## get summary statistics
x <- qchisq(sob[,1], n1 - 1)
y <- qchisq(sob[,2], n2 - 1)
z <- qnorm(sob[,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2))
  
sdv <- sqrt(sigma_1^2*x/((n1 - 1)*n1) + sigma_2^2*y/(n2*(n2 - 1)))
dfw <- sdv^4/(sigma_1^4*x^2/((n1 - 1)^3*n1^2) + sigma_2^4*y^2/((n2 - 1)^3*n2^2))

## p-value is the maximum of the p-values for each test  
pL <- pt((deltaL - z)/sdv, df = dfw)
pU <- pt((z - deltaU)/sdv, df = dfw)
pBoth <- pmax(pL, pU)

## add p-values to the previous data frame
df.set_ord <- df.set[order(df.set$index),]
df.set_ord$p <- pBoth

## get violin plots of p-values conditional on categorizations (should show separation)
violin_right <- ggplot(df.set_ord, aes(x=factor(type), y=p, 
                                       group=factor(type), colour = factor(type), fill = factor(type))) + 
  geom_violin() + theme_bw() +
  scale_colour_manual(
    name = expression(italic(n)),
    values = c("#5E4FA2", "#3288BD", "forestgreen","#FDAE61","#D53E4F"),
    aesthetics = c("colour", "fill")
  ) + theme(legend.position = "none") + ylim(0,1) +
  geom_hline(yintercept=0.05, lty = 2) + coord_flip() +
  theme(plot.title = element_text(size=20,face="bold",
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  labs(y= bquote(italic(p)*'-value ('*italic(n)*' = 8)')) +
  labs(x= "") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 5, b = 0, l = 0))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## combine left and right plots
Fig3 <- plot_grid(fig_left + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                      violin_right + theme(plot.margin=unit(c(0.75,0.5,0.75,0.5),"cm")), rel_widths = c(1,1))

# output as .pdf file for the article
pdf(file = "Figure3.pdf",   # The directory you want to save the file in
    width = 11.025, # The width of the plot in inches (12.41)
    height = 4.9875) # The height of the plot in inches (10.7)

Fig3

dev.off()
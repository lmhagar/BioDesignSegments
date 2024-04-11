## code to reproduce Figure SM2 in the online supplement

## BEGIN SETUP ##
require(pracma)

## code adapted from Shieh et al. (2022) to confirm power estimates in Table 1
## use sample sizes from Table 1
for (ii in c(3, 5, 8, 10, 15, 20, 30, 40, 50, 60)){
n1=ii #sample sizes
n2=ii
alpha=0.05 #type I error rate
del=19.2 #equivalence bound
mu1=96 #group means
mu2=92
sigma1=15 #group standard deviations
sigma2=18
#USER SPECIFICATIONS PORTION
mud=mu1-mu2
sigsq1=sigma1^2
sigsq2=sigma2^2
## default in their code is 100 here
## implementation of Simpson's rule
numintb=100
lb=numintb+1
dd=1e-10
coevecb<-c(1,rep(c(4,2),numintb/2-1),4,1)
intlb=(1-dd-dd)/numintb
bvec=dd+intlb*(0:numintb)
numintc=5000
lc=numintc+1
cl=1e-10
coevecc<-c(1,rep(c(4,2),numintc/2-1),4,1)
sigsqd=sigsq1/n1+sigsq2/n2
sigmad=sqrt(sigsqd)
df1=n1-1
df2=n2-1
wbpdf<-(intlb/3)*coevecb*dbeta(bvec,df1/2,df2/2)
dft=df1+df2
sn1=sigsq1/n1
sn2=sigsq2/n2
g=(sn1/df1)*bvec+(sn2/df2)*(1-bvec)
q1=(sn1/df1)*bvec/g
q2=1-q1
dfvvec=1/(q1*q1/df1+q2*q2/df2)
wstavec=qt(1-alpha,dfvvec)
epvec=rep(0,lb)
for (i in seq(lb)) {
  cu=(del^2)/(g[i]*wstavec[i]^2)
  intc=cu-cl
  intlc=intc/numintc
  cvec=cl+intlc*(0:numintc)
  wcpdf<-(intlc/3)*coevecc*dchisq(cvec,dft)
  st=sqrt(cvec*g[i])*wstavec[i]
  epvec[i]<-sum(wcpdf*(pmax(0,pnorm((-st+del-mud)/sigmad)-pnorm((st-del
                                                          -mud)/sigmad))))
}
wsepower=sum(wbpdf*epvec)
ntotal=n1+n2
## print results
print(c(n1,round(wsepower,4)))
}

## this code block illustrates the consistency issues for n = 2
for (ii in c(2)){
  n1=ii #sample sizes
  n2=ii
  alpha=0.05 #type I error rate
  del=19.2 #equivalence bound
  mu1=96 #group means
  mu2=92
  sigma1=15 #group standard deviations
  sigma2=18
  #USER SPECIFICATIONS PORTION
  tic <- Sys.time()
  mud=mu1-mu2
  sigsq1=sigma1^2
  sigsq2=sigma2^2
  ## now we try multiple values here to make a more fine grid
  ## for Simpson's rule
  numintbs=c(100, 1000, 10000, 50000)
  for (kk in 1:length(numintbs)){
  numintb <- numintbs[kk]
  lb=numintb+1
  dd=1e-10
  coevecb<-c(1,rep(c(4,2),numintb/2-1),4,1)
  intlb=(1-dd-dd)/numintb
  bvec=dd+intlb*(0:numintb)
  numintc=5000
  lc=numintc+1
  cl=1e-10
  coevecc<-c(1,rep(c(4,2),numintc/2-1),4,1)
  sigsqd=sigsq1/n1+sigsq2/n2
  sigmad=sqrt(sigsqd)
  df1=n1-1
  df2=n2-1
  wbpdf<-(intlb/3)*coevecb*dbeta(bvec,df1/2,df2/2)
  dft=df1+df2
  sn1=sigsq1/n1
  sn2=sigsq2/n2
  g=(sn1/df1)*bvec+(sn2/df2)*(1-bvec)
  q1=(sn1/df1)*bvec/g
  q2=1-q1
  dfvvec=1/(q1*q1/df1+q2*q2/df2)
  wstavec=qt(1-alpha,dfvvec)
  epvec=rep(0,lb)
  for (i in seq(lb)) {
    cu=(del^2)/(g[i]*wstavec[i]^2)
    intc=cu-cl
    intlc=intc/numintc
    cvec=cl+intlc*(0:numintc)
    wcpdf<-(intlc/3)*coevecc*dchisq(cvec,dft)
    st=sqrt(cvec*g[i])*wstavec[i]
    epvec[i]<-sum(wcpdf*(pmax(0,pnorm((-st+del-mud)/sigmad)-pnorm((st-del
                                                                   -mud)/sigmad))))
  }
  wsepower=sum(wbpdf*epvec)
  ntotal=n1+n2
  
  print(c(numintb, n1,round(wsepower,4)))
  toc <- Sys.time()
  toc - tic
  }
}

## this function illustrates integration with the pracma package
## to create Figure SM2

## modify Shieh et al.'s code to use numerical integration functions in R instead
## of Simpson's rule
intpower <- function(x, y, n1, n2, alpha, del, mu1, mu2, sigma1, sigma2){
  b <- x; v <- y
  mud=mu1-mu2
  sigsq1=sigma1^2
  sigsq2=sigma2^2
  sigsqd=sigsq1/n1+sigsq2/n2
  sigmad=sqrt(sigsqd)
  df1=n1-1
  df2=n2-1
  
  dft=df1+df2
  sn1=sigsq1/n1
  sn2=sigsq2/n2
  g=(sn1/df1)*b+(sn2/df2)*(1-b)
  q1=(sn1/df1)*b/g
  q2=1-q1
  dfvvec=1/(q1*q1/df1+q2*q2/df2)
  wstavec=qt(1-alpha,dfvvec)
  
  st <- sqrt(v*g)*wstavec
  
  t1 <- dbeta(b,df1/2,df2/2)
  t2 <- dchisq(v,dft)
  t3 <- pmax(0,pnorm((-st+del-mud)/sigmad)-pnorm((st-del -mud)/sigmad))
  return(t1*t2*t3)
}

## explore different sample sizes
nn <- c(5, 8, 10, 15)
for (j in 1:length(nn)){
  fun <- function(x,y) {intpower(x,y, nn[j],nn[j],0.05,19.2,92,96,15,18)}
  xmin <- 0; xmax <- 1
  ymin <- 0; ymax <- seq(35,1000,5)
  
  results <- NULL
  for (i in 1:length(ymax)){
    results[i] <- integral2(fun, xmin, xmax, ymin, ymax[i])$Q
  }
  
  write.csv(cbind(ymax, results), paste0("pracma", nn[j], ".csv"), row.names = FALSE)
}

pracma5 <- read.csv("pracma5.csv")
pracma8 <- read.csv("pracma8.csv")
pracma10 <- read.csv("pracma10.csv")
pracma15 <- read.csv("pracma15.csv")

## create a subplot for each sample size
t1 <- ggplot(data=pracma5, aes(x=ymax, y=results)) + theme_bw() + 
  geom_line() +
  labs(x= expression("Upper Bound ("*italic(n)*" = "*5*")"), y= "Estimated Power") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept = 0.1283, linetype = 2, color = "firebrick") + xlim(0,1000)

t2 <- ggplot(data=pracma8, aes(x=ymax, y=results)) + theme_bw() + 
  geom_line() +
  labs(x= expression("Upper Bound ("*italic(n)*" = "*8*")"), y= "Estimated Power") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept = 0.3801, linetype = 2, color = "firebrick") + xlim(0,1000)

t3 <- ggplot(data=pracma10, aes(x=ymax, y=results)) + theme_bw() + 
  geom_line() +
  labs(x= expression("Upper Bound ("*italic(n)*" = "*10*")"), y= "Estimated Power") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept = 0.5366, linetype = 2, color = "firebrick") + xlim(0,1000)

t4 <- ggplot(data=pracma15, aes(x=ymax, y=results)) + theme_bw() + 
  geom_line() +
  labs(x= expression("Upper Bound ("*italic(n)*" = "*15*")"), y= "Estimated Power") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept = 0.7699, linetype = 2, color = "firebrick") + xlim(0,1000)

## aggregate plots and output
row1 <- plot_grid(t1 + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                        t2 + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                        nrow = 1, rel_widths = c(1,1))

row2 <- plot_grid(t3 + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                  t4 + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                  nrow = 1, rel_widths = c(1,1))

body <- plot_grid(row1, row2,
                  rel_heights = c(1,1), nrow = 2)

pdf(file = "FigPracma.pdf",   # The directory you want to save the file in
    width = 11.75, # The width of the plot in inches (12.41)
    height = 6) # The height of the plot in inches (10.7)

body

dev.off()
## code to reproduce Figure 1 and Figure SM1 in the online supplement

## BEGIN SETUP ##
require(ggplot2)
require(cowplot)
require(ggpubr)
require(rgl)
require(qrng)
require(data.table)
require(imager)
require(magick)

## generate Sobol' sequence
sob <- sobol(n = 1024, d = 3, seed = 1, randomize = "digital.shift")
## load color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## now show the unit hypercube
plot3d(x = sob[c(159),1],
       y = sob[c(159),2],
       z = sob[c(159),3], 
       col = c(cbbPalette[2]), xlim = c(0,1), ylim = c(0,1), zlim = c(0,1), size =10,
       xlab = "", ylab = "", zlab = "",
       axes = FALSE)
axis3d('y--',nticks =0,labels = FALSE)
axis3d('y-+',nticks=0,labels = FALSE)
axis3d('y+-',nticks=1,labels = TRUE, cex = 2)
mtext3d(expression(italic(u)[2]), "y+-", line = 1.5, cex = 2)
axis3d('y++',nticks=0,labels = FALSE)
axis3d('x--',nticks=1,labels = TRUE, cex = 2)
mtext3d(expression(italic(u)[1]), "x--", line = 1.5, cex = 2)
axis3d('x-+',nticks=0,labels = FALSE)
axis3d('x+-',nticks=0,labels = FALSE)
axis3d('x++',nticks=0,labels = FALSE)
axis3d('z--',nticks=1,labels = TRUE, cex = 2)
mtext3d(expression(italic(u)[3]), "z--", line = 1.5, cex = 2)
axis3d('z-+',nticks=0,labels = FALSE)
axis3d('z+-',nticks=0,labels = FALSE)
axis3d('z++',nticks=0,labels = FALSE)

rgl.snapshot('fig1_pt.png', fmt = 'png')

## use points in Sobol' sequence to get z_mat (which holds d bar value for each point),
## sdv_mat (which holds standard error for each point), and dfw_mat (which holds degrees of
## freedom for each point); this will be used later for plotting
## define equivalence margin and significance level
dU <- 19.2; alpha <- 0.05
z_mat <- NULL
sdv_mat <- NULL
dfw_mat <- NULL
thres_mat <- NULL
for (j in 1:nrow(sob)){
  n1 <- seq(2,100); n2 <- seq(2,100)
  x <- qchisq(sob[j,1], n1 - 1)
  y <- qchisq(sob[j,2], n2 - 1)
  z <- qnorm(sob[j,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2))
  
  sdv <- sqrt(sigma_1^2*x/((n1 - 1)*n1) + sigma_2^2*y/(n2*(n2 - 1)))
  
  dfw <- sdv^4/(sigma_1^4*x^2/((n1 - 1)^3*n1^2) + sigma_2^4*y^2/((n2 - 1)^3*n2^2))
  
  thres <- (dU - abs(z))/qt(1-alpha, dfw)
  
  z_mat <- rbind(z_mat, z)
  sdv_mat <- rbind(sdv_mat, sdv)
  dfw_mat <- rbind(dfw_mat, dfw)
  thres_mat <- rbind(thres_mat, thres)
}

## create data frame with normal distribution for d bar
## true mu and sigma parameters are taken from illustrative example
## sample sizes of 8 in each group are used to create the figure
delta <- -4; n1 <- 8; n2 <- 8; sigma_1 <- 18; sigma_2 <- 15
df.bard <- data.frame(x = seq(-28, 20, by = 0.01),
                      y = dnorm(seq(-28, 20, by = 0.01), delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2)))

plot.bard1 <- ggplot(df.bard, aes(x=x, y=y)) + theme_bw() +
  geom_line(color="black", size = 1) +
  labs(x= expression(bar(italic(d))), y= " ") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  geom_segment(aes(x = qnorm(sob[159,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2)), y = 0, 
                   xend = qnorm(sob[159,3], delta, sqrt(sigma_1^2/n1 + sigma_2^2/n2)), yend = Inf), color = cbbPalette[2], size = 1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_x_continuous(breaks = seq(-28,28,length.out = 5)) +
  scale_y_continuous(breaks = seq(0, 0.05, length.out = 3))
plot.bard1

## create data frame with the true distribution of s1
df.s1 <- data.frame(x = seq(0,40, by = 0.01),
                    y = dchisq(seq(0,40, by = 0.01)^2*(n1-1)/sigma_1^2, n1-1)*(2*seq(0,40, by = 0.01)*(n1-1)/sigma_1^2))

plot.s1 <- ggplot(df.s1, aes(x=x, y=y)) + theme_bw() +
  geom_line(color="black", size = 1) +
  labs(x= expression(italic(s)[1]), y= " ") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  geom_segment(aes(x = sqrt(sigma_1^2*qchisq(sob[159,1], n1-1)/(n1-1)), y = 0, 
                   xend = sqrt(sigma_1^2*qchisq(sob[159,1], n1-1)/(n1-1)), yend = Inf), color = cbbPalette[2], size = 1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_x_continuous(breaks = seq(0,40,20)) +
  scale_y_continuous(breaks = seq(0, 0.085, length.out = 3))
plot.s1

## repeat the process for s2
df.s2 <- data.frame(x = seq(0,40, by = 0.01),
                    y = dchisq(seq(0,40, by = 0.01)^2*(n2-1)/sigma_2^2, n2-1)*(2*seq(0,40, by = 0.01)*(n2-1)/sigma_2^2))

plot.s2 <- ggplot(df.s2, aes(x=x, y=y)) + theme_bw() +
  geom_line(color="black", size = 1) +
  labs(x= expression(italic(s)[2]), y= " ") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  geom_segment(aes(x = sqrt(sigma_2^2*qchisq(sob[159,2], n2-1)/(n2-1)), y = 0, 
                   xend = sqrt(sigma_2^2*qchisq(sob[159,2], n2-1)/(n2-1)), yend = Inf), color = cbbPalette[2], size = 1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_x_continuous(breaks = seq(0,40,20)) +
  scale_y_continuous(breaks = seq(0, 0.1, 0.05))
plot.s2

## plot the triangular rejection region corresponding to this sample
dt.triangle <- data.table(group = c(1,1,1), polygon.x = c(-19.2,0,19.2), polygon.y = c(0,19.2/qt(0.95, dfw_mat[159,7]),0))

df.extend = data.frame(x = z_mat[159,7], y = sdv_mat[159,7]) 

p <- ggplot() + theme_bw() + geom_polygon(
  data = dt.triangle
  ,aes(
    x=polygon.x
    ,y=polygon.y
    ,group=group
  ),
  alpha = 0, col = cbbPalette[2], size = 1, linetype = "longdash"
) +
  geom_point(data = df.extend, aes(x = x, y = y), col = cbbPalette[2], size = 2) +
  geom_point(aes(x = 0.25, y = 11.3), col = "white", size = 2) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  labs(x= expression(bar(italic(d)))) +
  labs(y= expression(italic(se))) +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 5, r = 0, b = 0, l = 0))) +
  theme(axis.title.y = element_text(size = 18, margin = margin(t = 0, r = -10, b = 0, l = 0))) +
  scale_x_continuous(breaks = seq(-19.2,19.2,length.out = 3),
                     labels = c(expression(delta[L]), rep("", 1), expression(delta[U]))) +
  scale_y_continuous(breaks = seq(0, (3/2)*19.2/qt(0.95, dfw_mat[159,7]),length.out = 4),
                     labels = c(0, "", expression(frac(delta[U] - delta[L],2*t[1-alpha](nu))), ""))
p

## create composite figure
img1 <- ggplot()
image <- load.image("fig1_pt.png")

img1 <- ggdraw(img1) + draw_image(image)

fig.col1 <- plot_grid(plot.s1 + theme(plot.margin=unit(c(0.2,0.5,0.2,0.1),"cm")),
                      plot.s2 + theme(plot.margin=unit(c(0.2,0.5,0.2,0.1),"cm")), 
                      plot.bard1 + theme(plot.margin=unit(c(0.2,0.5,0.2,0.1),"cm")),
                      nrow = 3)
fig.col2 <- plot_grid(p)
fig <- plot_grid(plot_grid(NULL, img1, NULL, nrow =3, rel_heights = c(0.79, 2.15, 0.79)), fig.col1, 
                 plot_grid(NULL, fig.col2, NULL, nrow =3, rel_heights = c(0.2, 1, 0.05)), ncol = 3, rel_widths = c(0.75, 0.85, 1.25))

pdf(file = "Figure1.pdf",   # The directory you want to save the file in
    width = 9.75, # The width of the plot in inches (12.41)
    height = 4.125) # The height of the plot in inches (10.7)

fig

dev.off()

## code to reproduce Figure SM1

## this code indicates that points 549 and 621 have multiple intersections
reject <- thres_mat - sdv_mat >= 0
diffReject <- t(apply(reject, 1, diff))
fn_lt0 <- function(x){return(sum(x < 0))}
indices <- apply(diffReject, 1, fn_lt0)

## create data frame with se and lambda functions for point 549 to 
## illustrate multiple intersections
df_549 <- data.frame(n = rep(seq(2,100),2), y = c(thres_mat[549,], sdv_mat[549,]), 
                     type = c(rep("1Thres",99), rep("2Sdv",99)))

## expanded view
plotA1a <- ggplot(df_549, aes(x=n, y=y, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
  geom_line() +
  # geom_point() +
  scale_color_manual(name = "", rep(c(expression("  "*Lambda[italic(r)]^(italic(n)*","*italic(q))*"  "), 
                                      expression(paste(italic(se)[italic(r)]^(italic(n)*","*italic(q))))),1),
                     values = rep(c("firebrick", "steelblue"),1)) +
  scale_linetype_manual(name = "", rep(c(expression("  "*Lambda[italic(r)]^(italic(n)*","*italic(q))*"  "), 
                                         expression(paste(italic(se)[italic(r)]^(italic(n)*","*italic(q))))),1),
                        values = rep(c("solid", "dashed"),1)) +
  labs(color  = "", linetype = "") +
  labs(x= bquote(italic(n)), y= " ") +
  theme(plot.title = element_text(size=20,face="bold",
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.text=element_text(size=18)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.position="bottom")

## close up view to illustrate intersections
plotA1b <- ggplot(df_549[c(1:9, 100:108),], aes(x=n, y=y, color = as.factor(type), linetype = as.factor(type))) + theme_bw() +
  geom_line() +
  geom_point() +
  scale_color_manual(name = "", labels = rep(c(expression("  "*Lambda[italic(r)]^(italic(n)*","*italic(q))*"  "), 
                                               expression(paste(italic(se)[italic(r)]^(italic(n)*","*italic(q))))),1),
                     values = rep(c("firebrick", "steelblue"),1)) +
  scale_linetype_manual(name = "", labels = rep(c(expression("  "*Lambda[italic(r)]^(italic(n)*","*italic(q))*"  "), 
                                                  expression(paste(italic(se)[italic(r)]^(italic(n)*","*italic(q))))),1),
                        values = rep(c("solid", "dashed"),1)) +
  labs(color  = "", linetype = "") +
  labs(x= bquote(italic(n)), y= " ") +
  theme(plot.title = element_text(size=20,face="bold",
                                  margin = margin(t = 0, 0, 5, 0))) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18)) +
  theme(legend.text=element_text(size=18)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(legend.position="bottom") +
  ylim(2.198109, 9.014143)

FigA1 <- plot_grid(plotA1a + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                   plotA1b + theme(legend.position="none") + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), rel_widths = c(1,1))

## extract common legend and combine left and right plots
mylegend<-get_legend(plotA1a)
FigA11 <- plot_grid(FigA1, mylegend, nrow=2, rel_heights = c(10,1))

# output as .pdf file for the article
pdf(file = "FigureSM1.pdf",   # The directory you want to save the file in
    width = 11.025, # The width of the plot in inches (12.41)
    height = 4.9) # The height of the plot in inches (10.7)

FigA11

dev.off()
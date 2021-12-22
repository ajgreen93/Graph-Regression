#------------------------------#
# Script to create warmup example
#------------------------------#
setwd("~/Statistics/graph_regression/code")
source("graph.R")
source("sample.R")
library(Matrix)
library(magrittr)
library(RANN)
library(reshape2)
library(RSpectra)
setwd("~/Statistics/graph_regression/talks/laplacian_smoothing/MURI/figures/")

# Data
set.seed(4444)
n <- 200
d <- 2
sampler <- function(n)
{
  sample_x <- make_sample_gaussian_mixture(d = 1, mu = c(-1,1), pi = c(.5,.5),sigma = c(.5,.5),truncate = T)
  x <- sample_x(n)
  y <- sapply(x, FUN = function(x_i){runif(1,- (sqrt(x_i^4*(1 - x_i^2)) + .2), sqrt(x_i^4*(1 - x_i^2)) + .2)}) # dumbbell
  X <- matrix(c(x,y),ncol = 2)
}
X <- sampler(n)
f0 <- function(x){2*x[1]}
Y <- apply(X,1,f0) + rnorm(n,0,.5)

# Compute graph, keep largest connected component.
radius <- c(.15)
kernel_fxn <- function(s){exp(-2*s^2)}
G <- neighborhood_graph(X,radius,kernel = kernel_fxn)
X <- X[which(rowSums(G) > 0),]
Y <- Y[which(rowSums(G) > 0)]
G <- neighborhood_graph(X,radius,kernel = kernel_fxn) # one connected component
adj_list <- which(G > 0, arr.ind = T)
adj_list <- cbind(adj_list, G[adj_list])

# 1. Plot X and Y
coul <- colorRampPalette(c("blue","red"))
coul <- coul(25)
colors <- coul[cut(Y,25)]
pdf("warmup_example_1.pdf")
par(pty="s",mar = c(1,1,1,1))
plot(X[,1],X[,2], 
     ylim = c(-1.1,1.1),xlim = c(-1.1,1.1),
     axes = F,xlab = "",ylab = "",type = "n")
points(X[,1],X[,2],pch = 21, bg = colors, col = "black",cex = 2)
dev.off()

# 2. Plot neighborhood graph.

# Add more colors to this palette :
coul <- colorRampPalette(c("lightblue", "blue"))
coul <- coul(25)
colors <- coul[cut(log(adj_list[,3]),25)]

pdf("warmup_example_2.pdf")
par(pty="s",mar = c(1,1,1,1))
plot(X[,1],X[,2], 
     ylim = c(-1.1,1.1),xlim = c(-1.1,1.1),
     axes = F,xlab = "",ylab = "",type = "n")
adj_indx <- which(1:n %in% adj_list[,-3])
points(X[adj_indx,1],X[adj_indx,2],pch = 21, bg = "grey62", col = "black",cex = 2)
for(ii in 1:nrow(adj_list))
{
  segments(x0 = X[adj_list[ii,1],1],y0 = X[adj_list[ii,1],2],
           x1 = X[adj_list[ii,2],1],y1 = X[adj_list[ii,2],2],
           col = colors[ii],
           lwd = 2)
}
dev.off()

# 3. Plot estimate

# Estimate
n <- nrow(X)
L <- Laplacian(G)
lambda <- 1/rev(eigen(L)$values)[10]
f <- solve(diag(n) + lambda*L,Y) %>% as.numeric()
coul <- colorRampPalette(c("blue","red"))
coul <- coul(25)
colors <- coul[cut(f,25)]
pdf("warmup_example_3.pdf")
par(pty="s",mar = c(1,1,1,1))
plot(X[,1],X[,2], 
     ylim = c(-1.1,1.1),xlim = c(-1.1,1.1),
     axes = F,xlab = "",ylab = "",type = "n")
points(X[,1],X[,2],pch = 21, bg = colors, col = "black",cex = 2)
dev.off()

# 4. Plot Nystrom extension

# Nystrom
X_new <- expand.grid(seq(-1.25,1.25,length.out = 200),
                     seq(-.75,.75,length.out = 200))
Y_blur <- lambda * G  %*% solve(diag(n) + lambda*diag(rowSums(G))) %*% Y
f <- solve(diag(n) + lambda*L,Y_blur) %>% as.numeric()

G_new <- neighborhood_graph(X,radius,X_new = X_new,kernel = kernel_fxn) %>% t()
D_new <- rowSums(G_new)
theta_new <- lambda/(1 + lambda*D_new) * G_new %*% (f  + Y/(1 + lambda*rowSums(G)))

# Special plotting
u1 <- unique(X_new[,1])
u2 <- unique(X_new[,2])
indx <- matrix(nrow = nrow(X_new),ncol = 2)
indx[,1] <- sapply(X_new[,1], FUN = function(x){which(u1 == x)})
indx[,2] <- sapply(X_new[,2], FUN = function(x){which(u2 == x)})
z <- matrix(NA,ncol = length(unique(X_new[,1])),
            nrow = length(unique(X_new[,2])))
z[indx] <- theta_new
z[which(z == 0)] <- NA
pdf("warmup_example_4.pdf")
image(u1,u2,z = z,
      ylim = c(-1.1,1.1),
      col = colorRampPalette(c("blue", "red"))(25),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "",
      axes = F)
dev.off()

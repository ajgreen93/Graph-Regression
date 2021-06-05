#--------------------------------------------------#
# Script to generate plots analyzing impact of LE tuning parameters.
#--------------------------------------------------#
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("plot_methods.R")
source("misc.R")

# User entered information.
data_directory <- "data/thesis/cluster_assumption/stepfunction"
plot_directory <- "plots/thesis/cluster_assumption/tuning/stepfunction" # Please change this to whichever directory you prefer.
if(!exists(plot_directory)){dir.create(plot_directory)}

plot_n <- 1111                                 # Which sample size would you like to use?
plot_methods <- c("laplacian_eigenmaps",
                  "laplacian_smoothing",
                  "laplacian_eigenmaps_plus_kernel_smoothing",
                  "spectral_projection",
                  "least_squares",
                  "kernel_smoothing")

# Load data
load(file.path(data_directory,"configs.R"))
d <- configs$d; s <- configs$s; ns <- configs$ns; methods <- configs$methods
load(file.path(data_directory,"best_fits_by_method.R"))
load(file.path(data_directory,"thetas.R"))
load(file.path(data_directory,"mse.R"))
if("test_mse.R" %in% list.files(data_directory)) load(file.path(data_directory,"test_mse.R"))
load(file.path(data_directory,"Xs.R"))
load(file.path(data_directory,"f0s.R"))
load(file.path(data_directory,"Ys.R"))

stopifnot(plot_n %in% ns)
mse <- mse[ns %in% plot_n][[1]][names(methods) %in% plot_methods]
if(exists("test_mse")) test_mse <- test_mse[ns %in% plot_n][[1]][names(methods) %in% plot_methods]
thetas <- thetas[ns %in% plot_n][[1]][names(methods) %in% plot_methods]

# Plotting parameters for all plots
title <- "Stepfunction."

## Plot 1: Mean squared error as a function of K.

# Find best parameters
best_parameters <- find_best_parameters(list(mse),list(thetas))[[1]]
alg_indx <- sapply(thetas,FUN = function(theta){"K" %in% names(theta)}) %>% which()
Ks <- sapply(thetas[alg_indx],FUN = function(theta){unique(theta$K)})

# Mse as a function of K, for best other parameters.
plot_mse <- matrix(nrow = nrow(Ks),ncol = ncol(Ks))
colnames(plot_mse) <- names(methods)[alg_indx]
for(jj in 1:length(alg_indx))
{
  ii <- alg_indx[jj]
  best_parameters_ii <- best_parameters[[ii]]
  mse_ii <- mse[[ii]]
  thetas_ii <- thetas[[ii]]
  if(exists("test_mse")) test_mse_ii <- test_mse[[ii]]
  if(is.null(names(best_parameters_ii))){
    # Spectral projection or least squares
    plot_mse[,jj] <- rowMeans(mse_ii)
  } else if(all(names(best_parameters_ii) == c("r","K"))){
    # Laplacian eigenmaps
    plot_mse[,jj] <- rowMeans(mse_ii[thetas_ii$r == best_parameters_ii$r,])
  } else if(all(names(best_parameters_ii) == c("r","K","h"))){
    # Laplacian eigenmaps plus kernel smoothing
    plot_mse[,jj] <- rowMeans(test_mse_ii[thetas_ii$r == best_parameters_ii$r & thetas_ii$h == best_parameters_ii$h,])
  }
}

# Plotting parameters
xlims <- c(min(Ks),max(Ks))
ylims <- c(min(plot_mse),max(plot_mse))
cols <- c("red","blue","green")
stopifnot(ncol(plot_mse) <= 3)

plot_name <- paste0("mse_by_number_of_eigenvectors_",plot_n,"n_",d,"d_",s,"s.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
# shadow plot
plot(x = Ks[,1], xlim = xlims, ylim = ylims, xlab = "Number of eigenvectors", ylab = "Mean Squared Error", 
     main = title,
     cex.main = 2.5, cex.lab = 2.5, cex.axis = 2, type = "n")

# add points and lines
for(jj in 1:ncol(plot_mse))
{
  points(x = Ks[,jj], y =  plot_mse[,jj], col = cols[jj], pch = 20)
  lines(x = Ks[,jj], y =  plot_mse[,jj], col = cols[jj], lwd = 1.5)
}
grid(lwd = 2)
dev.off()

## Plot 2: Mean squared error as a function of radius.
r_methods <- c("kernel_smoothing")
alg_indx_r <- which(names(methods) %in% r_methods)
rs <- sapply(thetas[alg_indx_r],FUN = function(theta){unique(theta$r)})

# Mse as a function of r, for best other parameters.
plot_mse_r <- matrix(nrow = nrow(rs),ncol = ncol(rs))
for(jj in 1:length(alg_indx_r))
{
  ii <- alg_indx_r[jj]
  best_parameters_ii <- best_parameters[[ii]]
  mse_ii <- mse[[ii]]
  thetas_ii <- thetas[[ii]]
  if(exists("test_mse")) test_mse_ii <- test_mse[[ii]]
  if(is.null(names(best_parameters_ii))){
    # Kernel smoothing
    plot_mse_r[,jj] <- rowMeans(mse_ii)
  } else if(all(names(best_parameters_ii) == c("r","K"))){
    # Laplacian eigenmaps
    plot_mse_r[,jj] <- unique(rowMeans(mse_ii[thetas_ii$K == best_parameters_ii$K,]))
  } else if(all(names(best_parameters_ii) == c("r","K","h"))){
    # Laplacian eigenmaps plus kernel smoothing
    plot_mse_r[,jj] <- rowMeans(test_mse_ii[thetas_ii$K == best_parameters_ii$K & thetas_ii$h == best_parameters_ii$h,])
  } else if(all(names(best_parameters_ii) == c("r","rho"))){
    # Laplacian smoothing
    plot_mse_r[,jj] <- unique(rowMeans(mse_ii[thetas_ii$rho == best_parameters_ii$rho,]))
  }
}

# Mse as a function of h, for best other parameters.
alg_indx_h <- sapply(thetas,FUN = function(theta){"h" %in% names(theta)}) %>% which()
if(length(alg_indx_h) > 0)
{
  hs <- sapply(thetas[alg_indx_h],FUN = function(theta){unique(theta$h)})
  plot_mse_h <- matrix(nrow = nrow(hs),ncol = ncol(hs))
  for(jj in 1:length(alg_indx_h))
  {
    ii <- alg_indx_h[jj]
    best_parameters_ii <- best_parameters[[ii]]
    mse_ii <- mse[[ii]]
    thetas_ii <- thetas[[ii]]
    test_mse_ii <- test_mse[[ii]]
    if(all(names(best_parameters_ii) == c("r","K","h"))){
      # Laplacian eigenmaps plus kernel smoothing
      plot_mse_h[,jj] <- rowMeans(test_mse_ii[thetas_ii$K == best_parameters_ii$K & thetas_ii$r == best_parameters_ii$r,])
    }
  }
  plot_mse <- cbind(plot_mse_r,plot_mse_h)
  rs <- cbind(rs,hs)
  colnames(plot_mse) <- c(paste("r",names(methods)[alg_indx_r],sep = "_"),
                          paste("h",names(methods)[alg_indx_h],sep = "_"))
} else{
  plot_mse <- plot_mse_r
  colnames(plot_mse) <- c(paste("r",names(methods)[alg_indx_r],sep = "_"))
}

# Plotting parameters
xlims <- c(min(rs),max(rs))
ylims <- c(min(plot_mse),max(plot_mse))
cols <- c("purple")
stopifnot(ncol(plot_mse) <= 3)

plot_name <- paste0("mse_by_radius_",plot_n,"n_",d,"d_",s,"s.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
# shadow plot
plot(x = rs[,1], xlim = xlims, ylim = ylims, xlab = "Radius", ylab = "Mean Squared Error", 
     main = title,
     cex.main = 2.5, cex.lab = 2.5, cex.axis = 2, type = "n")

# add points and lines
for(jj in 1:ncol(plot_mse))
{
  if(startsWith(colnames(plot_mse)[jj],"h_"))
  {
    points(x = rs[,jj], y =  plot_mse[,jj], col = cols[jj], pch = 17, cex = .9)
    lines(x = rs[,jj], y =  plot_mse[,jj], col = cols[jj], lwd = 1.5, lty = 2)
  } else{
    points(x = rs[,jj], y =  plot_mse[,jj], col = cols[jj], pch = 20)
    lines(x = rs[,jj], y =  plot_mse[,jj], col = cols[jj], lwd = 1.5)
  }
}
grid(lwd = 2)
dev.off()
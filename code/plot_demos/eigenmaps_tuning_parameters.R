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
data_directory <- "data/laplacian_eigenmaps/tuning/20210525094744"
plot_directory <- "plots/laplacian_eigenmaps/tuning" # Please change this to whichever directory you prefer.
if(!exists(plot_directory)){dir.create(plot_directory)}
plot_n <- 1000                                       # Which sample size would you like to use?

# Load data
load(file.path(data_directory,"configs.R"))
d <- configs$d; s <- configs$s; ns <- configs$ns; methods <- configs$methods
load(file.path(data_directory,"best_fits_by_method.R"))
load(file.path(data_directory,"thetas.R"))
load(file.path(data_directory,"mse.R"))
load(file.path(data_directory,"Xs.R"))
load(file.path(data_directory,"f0s.R"))
load(file.path(data_directory,"Ys.R"))

stopifnot(plot_n %in% ns)
mse <- mse[plot_n %in% ns]
methods <- methods[plot_n %in% ns]
thetas <- thetas[plot_n %in% ns]

# Plotting parameters for all plots
title <- paste0("d = ",d,", s = ",s,".")

## Plot 1: Mean squared error as a function of K.

# Find best radius
best_r <- find_best_radius(mse,thetas)
alg_indx <- sapply(thetas[[1]],FUN = function(theta){"K" %in% names(theta)}) %>% which()
Ks <- sapply(thetas[[1]][alg_indx],FUN = function(theta){unique(theta$K)})

# Parameters
plot_mse <- matrix(nrow = nrow(Ks),ncol = ncol(Ks))
for(ii in alg_indx)
{
  if(is.na(best_r[ii]))
  {
    plot_mse[,ii] <- rowMeans(mse[[1]][[ii]])
  } else{
    best_r_row <- which(thetas[[1]][[ii]]$r %in% best_r)
    plot_mse[,ii] <- rowMeans(mse[[1]][[ii]][best_r_row,])
  }
}

# Plotting parameters
xlims <- c(min(Ks),max(Ks))
ylims <- c(min(plot_mse),max(plot_mse))
cols <- c("red","blue","green")
stopifnot(ncol(plot_mse) <= 3)

pdf(file.path(plot_directory,"mse_by_number_of_eigenvectors.pdf"))
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

# Parameters
alg_indx <- which(names(configs$methods) == "laplacian_eigenmaps")
n_indx <- length(ns)
rs <- unique(thetas[[n_indx]][[alg_indx]]$r)
best_K <- find_best_K(list(mse[[n_indx]]),list(thetas[[n_indx]]))[1]
Ks <- round(seq(5,40,length.out = 9))

# Data for plotting
plot_mse <- mse[[n_indx]][[alg_indx]][which(thetas[[n_indx]][[alg_indx]]$K %in% Ks),] %>%
  rowMeans()
plot_rs <- thetas[[n_indx]][[alg_indx]][which(thetas[[n_indx]][[alg_indx]]$K %in% Ks),"r"]
plot_Ks <- thetas[[n_indx]][[alg_indx]][which(thetas[[n_indx]][[alg_indx]]$K %in% Ks),"K"]
plot_mat <- cbind(plot_mse,plot_rs,plot_Ks); colnames(plot_mat) <- c("mse","r","K")

# Plotting parameters
xlims <- c(min(rs),max(rs))
ylims <- c(min(plot_mse),max(plot_mse))
cols <- viridis(length(Ks))

pdf(file.path(plot_directory,"mse_by_radius.pdf"))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.

# shadow plot
title <- paste0("d = ",d,", s = ",s,".")
plot(x = plot_mat[,"r"], xlim = xlims, ylim = ylims, xlab = "Radius", ylab = "Mean squared error", 
     main = title, cex.main = 2.5, cex.lab = 2.5, cex.axis = 2, type = "n")

# add points and lines
for(jj in 1:length(Ks))
{
  plot_indx <- which(plot_mat[,"K"] == Ks[jj])
  points(x = plot_mat[plot_indx,"r"], y =  plot_mat[plot_indx,"mse"], col = cols[jj], pch = 20)
  lines(x = plot_mat[plot_indx,"r"], y =  plot_mat[plot_indx,"mse"], col = cols[jj], lwd = 1.5)
}
grid(lwd = 2)
dev.off()
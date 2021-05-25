#--------------------------------------------------#
# Script to generate plots analyzing impact of LE tuning parameters.
#--------------------------------------------------#
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("plot_methods.R")

# User entered information.
data_directory <- "data/laplacian_eigenmaps/eigenmaps_parameters_2s"
plot_directory <- "plots/laplacian_eigenmaps/eigenmaps_parameters_2s" # Please change this to whichever directory you prefer.
if(!exists(plot_directory)){dir.create(plot_directory)}

# Load data
load(file.path(data_directory,"configs.R"))
d <- configs$d; ns <- configs$ns; methods <- configs$methods
load(file.path(data_directory,"best_fits_by_method.R"))
load(file.path(data_directory,"thetas.R"))
load(file.path(data_directory,"mse.R"))
load(file.path(data_directory,"Xs.R"))
load(file.path(data_directory,"f0s.R"))
load(file.path(data_directory,"Ys.R"))

# Plotting parameters for all plots
title <- paste0("d = ",d,", s = ",s,".")

## Plot 1: Mean squared error as a function of K.

# Parameters
alg_indx <- which(names(configs$methods) == "laplacian_eigenmaps")
n_indx <- length(ns)
all_rs <- unique(thetas[[n_indx]][[alg_indx]]$r)
rs <- all_rs[round(seq(1,length(all_rs),length.out = 9))]
Ks <- unique(thetas[[n_indx]][[alg_indx]][,"K"])

# Data for plotting
plot_mse <- mse[[n_indx]][[alg_indx]][which(thetas[[n_indx]][[alg_indx]]$r %in% rs),] %>%
  rowMeans()
plot_rs <- thetas[[n_indx]][[alg_indx]][which(thetas[[n_indx]][[alg_indx]]$r %in% rs),"r"]
plot_Ks <- thetas[[n_indx]][[alg_indx]][which(thetas[[n_indx]][[alg_indx]]$r %in% rs),"K"]
plot_mat <- cbind(plot_mse,plot_rs,plot_Ks); colnames(plot_mat) <- c("mse","r","K")

# Plotting parameters
xlims <- c(min(Ks),max(Ks))
ylims <- c(min(plot_mse),max(plot_mse))
cols <- viridis(length(rs))

pdf(file.path(plot_directory,"mse_by_number_of_eigenvectors.pdf"))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
# shadow plot
plot(x = plot_mat[,"K"], xlim = xlims, ylim = ylims, xlab = "Number of eigenvectors", ylab = "Mean Squared Error", 
     main = title,
     cex.main = 2.5, cex.lab = 2.5, cex.axis = 2, type = "n")

# add points and lines
for(jj in 1:length(rs))
{
  plot_indx <- which(plot_mat[,"r"] == rs[jj])
  points(x = plot_mat[plot_indx,"K"], y =  plot_mat[plot_indx,"mse"], col = cols[jj], pch = 20)
  lines(x = plot_mat[plot_indx,"K"], y =  plot_mat[plot_indx,"mse"], col = cols[jj], lwd = 1.5)
}
grid(lwd = 2)
dev.off()

## Plot 1: Mean squared error as a function of K.

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
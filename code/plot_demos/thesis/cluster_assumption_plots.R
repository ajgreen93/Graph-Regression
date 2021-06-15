#--------------------------------------------------#
# Script to generate plots for cluster assumption regression fxn in thesis.
#--------------------------------------------------#
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("misc.R")
source("plot_methods.R")
source("sample.R")

# User entered information.
data_directory <- "data/thesis/cluster_assumption/20210608233004"
plot_directory <- "plots/thesis/cluster_assumption/gaussian_mixture_stepfunction" # Please change this to whichever directory you prefer.
if(!exists(plot_directory)){dir.create(plot_directory,recursive = T)}

# Load data
load(file.path(data_directory,"configs.R"))
d <- configs$d; ns <- configs$ns; s <- configs$s; M <- configs$M;
methods <- configs$methods
load(file.path(data_directory,"best_fits_by_method.R"))
load(file.path(data_directory,"thetas.R"))
load(file.path(data_directory,"mse.R"))
load(file.path(data_directory,"Xs.R"))
load(file.path(data_directory,"f0s.R"))
load(file.path(data_directory,"Ys.R"))

# Subset data to methods you actually want to plot.
plot_methods <- c("laplacian_eigenmaps",
                  "laplacian_smoothing",
                  "least_squares",
                  "kernel_smoothing")
mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})
methods <- methods[names(methods) %in% plot_methods]

## Plots 1-3: Actual and estimated regression functions.
ii <- 1
n <- ns[ii]
# sep <- (log(n))^2/(2 * n) # TODO: Should have a way to grab this from data!
# domains <- list(c(0,1/2 - sep),c(1/2 + sep,1))
domains <- c(-1,1)

# Plot 1: Regression function
plot_name <- paste0("actual_function.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
plot_fxn(x = Xs[[ii]], 
         y = Ys[[ii]],
         f = apply(Xs[[ii]],1,f0s[[ii]]),
         domains = domains,
         title = "True function.")
dev.off()

# Plot 2: fits of Laplacian smoothing and eigenmaps
plot_name <- paste0("laplacian_fits.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
indx <- names(methods) %in% c("laplacian_eigenmaps","laplacian_smoothing")
plot_fxn(x = Xs[[ii]], 
         y = Ys[[ii]],
         f = (best_fits_by_method[[ii]][indx]),
         domains = domains,
         title = "Laplacian methods.",
         cols = c("red","darkgreen"))
dev.off()

# Plot 3: fits of spectral projection and kernel smoothing
plot_name <- paste0("smooth_fits.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
indx <- names(methods) %in% c("least_squares","kernel_smoothing")
plot_fxn(x = Xs[[ii]], 
         y = Ys[[ii]],
         f = (best_fits_by_method[[ii]][indx]), 
         domains = domains,
         title = "Kernel smoothing, least squares.",
         cols = c("blue","purple"))
dev.off()

## Plot 4: Mean squared error as a function of n.
plot_name <- paste0("mse_by_sample_size.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
plot_best_mse(methods,mse,sd = T, cols = c("red","darkgreen","blue","purple"),
              legend = F, title = "Mean squared error.", rate = F)
dev.off()
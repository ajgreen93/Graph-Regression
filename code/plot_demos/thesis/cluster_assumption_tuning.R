#--------------------------------------------------#
# Script to generate plots analyzing impact of tuning parameters.
#--------------------------------------------------#
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("plot_methods.R")
source("misc.R")

# User entered information.
data_directory <- "data/thesis/cluster_assumption/20210607144055"
plot_directory <- "plots/thesis/cluster_assumption/tuning/gaussian_mixture_signed_rt" # Please change this to whichever directory you prefer.
if(!exists(plot_directory)){dir.create(plot_directory)}

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

# Subset data to methods you actually want to plot.
plot_methods <- c("laplacian_eigenmaps",
                  "laplacian_smoothing",
                  "least_squares",
                  "kernel_smoothing")
mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})
thetas <- lapply(thetas,FUN = function(t){t[names(methods) %in% plot_methods]})

# Plot 1: Error as a function of K
plot_methods <- c("laplacian_eigenmaps","least_squares")
plot_mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})[1:9]
plot_thetas <- lapply(thetas,FUN = function(t){t[names(methods) %in% plot_methods]})[1:9]

cols <- c("red","blue")
title <- paste0("n = ",ns,".")

plot_name <- paste0("mse_by_K.pdf")
pdf(file.path(plot_directory,plot_name))
par(mfrow = c(2,2))
plot_mse_by_tuning_parameter(parameter_name = "K",
                             plot_methods,plot_thetas,plot_mse,
                             cols,title,
                             ylims = c(0,.2))
dev.off()

# Plot 2: Error as a function of r
plot_methods <- c("laplacian_eigenmaps","laplacian_smoothing","kernel_smoothing")
plot_mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})[1:9]
plot_thetas <- lapply(thetas,FUN = function(t){t[names(methods) %in% plot_methods]})[1:9]

cols <- c("red","darkgreen","purple")
title <- paste0("n = ",ns,".")

plot_name <- paste0("mse_by_r.pdf")
pdf(file.path(plot_directory,plot_name))
par(mfrow = c(2,2))
plot_mse_by_tuning_parameter(parameter_name = "r",
                             plot_methods,plot_thetas,plot_mse,
                             cols,title,
                             ylims = c(0,.2))
dev.off()




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
data_directory <- "data/thesis/mse/eigenfunction_1s_1d"
plot_directory <- "plots/thesis/mse/eigenfunction/1s_1d" # Please change this to whichever directory you prefer.
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
                  "spectral_projection")
mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})
thetas <- lapply(thetas,FUN = function(t){t[names(methods) %in% plot_methods]})

# Plot 1: Error as a function of K
plot_methods <- c("laplacian_eigenmaps","spectral_projection")
plot_methods <- plot_methods[plot_methods %in% names(methods)]
plot_mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})
plot_thetas <- lapply(thetas,FUN = function(t){t[names(methods) %in% plot_methods]})

cols <- c("red","blue")
title <- paste0("n = ",ns,".")

plot_name <- paste0("mse_by_K_",d,"d_",s,"s",".pdf")
pdf(file.path(plot_directory,plot_name))
par(mfrow = c(2,2))
plot_mse_by_tuning_parameter(parameter_name = "K",
                             plot_methods,plot_thetas,plot_mse,
                             cols,title, 
                             ylims = NULL)
dev.off()

# Plot 2: Error as a function of r
plot_methods <- c("laplacian_eigenmaps","laplacian_smoothing","kernel_smoothing")
plot_methods <- plot_methods[plot_methods %in% names(methods)]
plot_mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})
plot_thetas <- lapply(thetas,FUN = function(t){t[names(methods) %in% plot_methods]})

cols <- c("red","darkgreen","purple")
title <- paste0("n = ",ns,".")

plot_name <- paste0("mse_by_r_",d,"d_",s,"s",".pdf")
pdf(file.path(plot_directory,plot_name))
par(mfrow = c(2,2))
plot_mse_by_tuning_parameter(parameter_name = "r",
                             plot_methods,plot_thetas,plot_mse,
                             cols,title,
                             ylims = NULL)
dev.off()

# Plot 3: Error as a function of rho
plot_methods <- c("laplacian_smoothing")
plot_methods <- plot_methods[plot_methods %in% names(methods)]
plot_mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})
plot_thetas <- lapply(thetas,FUN = function(t){t[names(methods) %in% plot_methods]})

cols <- c("darkgreen")
title <- paste0("n = ",ns,".")

plot_name <- paste0("mse_by_rho_",d,"d_",s,"s",".pdf")
pdf(file.path(plot_directory,plot_name))
par(mfrow = c(2,2))
plot_mse_by_tuning_parameter(parameter_name = "rho",
                             plot_methods,plot_thetas,plot_mse,
                             cols,title,
                             xlims = c(quantile(unique(bind_rows(plot_thetas)$rho),0),
                                       quantile(unique(bind_rows(plot_thetas)$rho),.99)),
                             log_axis = "x")
dev.off()

# Plot 4: Error for eigenmaps by r and K
plot_methods <- c("laplacian_eigenmaps")
plot_methods <- plot_methods[plot_methods %in% names(methods)]
plot_mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})
plot_thetas <- lapply(thetas,FUN = function(t){t[names(methods) %in% plot_methods]})

parameter_names = c("r","K")

title <- paste0("n = ",ns,".")
plot_name <- paste0("mse_by_K_and_r_",d,"d_",s,"s",".pdf")
pdf(file.path(plot_directory,plot_name))
par(mfrow = c(2,2))
plot_mse_by_tuning_parameter_2d(parameter_names = c("r","K"),
                             plot_thetas,plot_mse,
                             title,
                             ylims = NULL,
                             log_axis = "")
dev.off()

# Plot 5: Error for smoothing by r and rho
plot_methods <- c("laplacian_smoothing")
plot_methods <- plot_methods[plot_methods %in% names(methods)]
plot_mse <- lapply(mse,FUN = function(m){m[names(methods) %in% plot_methods]})
plot_thetas <- lapply(thetas,FUN = function(t){t[names(methods) %in% plot_methods]})

parameter_names = c("r","rho")

title <- paste0("n = ",ns,".")
plot_name <- paste0("mse_by_rho_and_r_",d,"d_",s,"s",".pdf")
pdf(file.path(plot_directory,plot_name))
par(mfrow = c(2,2))
plot_mse_by_tuning_parameter_2d(parameter_names = c("r","rho"),
                                plot_thetas,plot_mse,
                                title,
                                ylims = c(quantile(unique(bind_rows(plot_thetas)$rho),0),
                                          quantile(unique(bind_rows(plot_thetas)$rho),.99)),
                                log_axis = "y")
dev.off()


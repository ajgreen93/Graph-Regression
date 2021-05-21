#--------------------------------------------------#
# Script to generate plots used for Figure 1 of Laplacian Eigenmaps paper.
#--------------------------------------------------#
library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("misc.R")
source("plot_methods.R")

# User entered information.
data_directory <- "data/laplacian_eigenmaps/mse/eigenmaps_mse_2d_2s"
plot_directory <- "plots/laplacian_eigenmaps/eigenmaps_mse" # Please change this to whichever directory you prefer.
if(!exists(plot_directory)){dir.create(plot_directory)}

load(file.path(data_directory,"configs.R"))
d <- configs$d; ns <- configs$ns; s <- configs$s; methods <- configs$methods
load(file.path(data_directory,"best_fits_by_method.R"))
load(file.path(data_directory,"thetas.R"))
load(file.path(data_directory,"mse.R"))
load(file.path(data_directory,"Xs.R"))
load(file.path(data_directory,"f0s.R"))
load(file.path(data_directory,"Ys.R"))

## Plot 1: Mean squared error as a function of n.
plot_name <- paste0("mse_by_sample_size_",d,"d_",s,"s.pdf")
pdf(file.path(plot_directory,plot_name))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.

## ALDEN GET RID OF
best_K <- find_best_K(mse,thetas)
best_radius <- find_best_radius(mse,thetas) 
##

plot_best_mse(methods,mse)
dev.off()
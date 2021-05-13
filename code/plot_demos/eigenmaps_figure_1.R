library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("plot_methods.R")

# User entered information.
data_directory <- "data/laplacian_eigenmaps/20210511125310"
plot_directory <- "plots/laplacian_eigenmaps/eigenmaps_approximates_spectral_projection" # Please change this to whichever directory you prefer.
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

## Plot 1: Regression fxn and fitted values.
plot1_indx <- 5
plot1_K <- 10
rs <- unique(thetas[[plot1_indx]][[1]]$r)
plot1_r <- rs[which.min(abs(rs - .02))]    # arbitrary choice

# Regression fxn
pdf(file.path(plot_directory,"f0.pdf"))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
plot_fxn(Xs[[plot1_indx]],Ys[[plot1_indx]],f = f0s[[plot1_indx]])
dev.off()

# Eigenmaps fitted values
pdf(file.path(plot_directory,"eigenmaps.pdf"))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
f <- subset(fits[[plot1_indx]][[1]], 
            thetas[[plot1_indx]][[1]]$r == plot1_r & 
            thetas[[plot1_indx]][[1]]$K == plot1_K)[1,]
            
plot_fxn(Xs[[plot1_indx]],Ys[[plot1_indx]],f,title = "Eigenmaps, K = 10",col = "red")
dev.off()

# Eigenmaps fitted values
pdf(file.path(plot_directory,"spectral_projection.pdf"))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
f <- subset(fits[[plot1_indx]][[2]],thetas[[plot1_indx]][[2]]$K == plot1_K)[1,]
plot_fxn(Xs[[plot1_indx]],Ys[[plot1_indx]],f,title = "Spectral Projection, K = 10",col = "blue")
dev.off()


## Plot 2: Mean squared error as a function of n.
pdf(file.path(plot_directory,"mse_by_sample_size.pdf"))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
plot_best_mse(methods,mse)
dev.off()

## Plot 3: Mean squared difference between estimators.
plot3K <- 10
plot_mat <- matrix(nrow = length(ns),ncol = 1)
for(ii in 1:length(ns))
{
  rs <- unique(thetas[[ii]][[1]]$r)
  r <-  rs[round(length(rs)/2)]
  fit_sp <- best_fits_by_method[[ii]][[2]]
  fit_le <- best_fits_by_method[[ii]][[1]]
  # fit_sp <- subset(fits[[ii]][[2]],thetas[[ii]][[2]]$K == plot3K)
  # fit_le <- subset(fits[[ii]][[1]],
  #                 thetas[[ii]][[1]]$K == plot3K & thetas[[ii]][[1]]$r == r)[1,]
  plot_mat[ii,] <-  mean( (fit_sp - fit_le)^2) 
}

# shadow plot
pdf(file.path(plot_directory,"squared_difference_by_sample_size.pdf"))
par(mar = c(5.1,6,4.1,4.1))
title <- "Difference."
plot(x = ns,y = plot_mat,
     log = "xy", xlab = "Sample size", ylab = "Mean squared difference", main = title,
     pch = 20,
     cex.main = 2.5, cex.lab = 2.5, cex.axis = 2)
lines(x = ns, y = plot_mat, lwd = 1.5)
dev.off()


## Plot 4: Mean squared error as a function of K.

# Parameters
plot4n <- ns[length(ns)]
plot4rs <- unique(thetas[[length(ns)]][[1]]$r)
plot4r <- plot4rs[round(length(plot4rs)/2)]
plot4Ks <- unique(thetas[[length(ns)]][[1]]$K)

# Data for plotting
plot4_mat<- matrix(nrow = length(plot4Ks),ncol = 2)
colnames(plot4_mat) <- c("le_mean_mse","sp_mean_mse")
plot4_mat[,"le_mean_mse"] <- rowMeans(subset(mse[[which(ns == plot4n)]][[1]],thetas[[length(ns)]][[1]]$r == plot4r))
plot4_mat[,"sp_mean_mse"] <- rowMeans(mse[[which(ns == plot4n)]][[2]])

# Plotting parameters
xlims <- c(min(plot4Ks),max(plot4Ks))
ylims <- c(min(plot4_mat),max(plot4_mat))
cols <- c("red","blue")

pdf(file.path(plot_directory,"mse_by_number_of_eigenvectors.pdf"))
par(mar = c(5.1,6,4.1,4.1)) # Prevent the left hand side from getting cut off.
# shadow plot
plot(x = plot4Ks, xlim = xlims, ylim = ylims, xlab = "K", ylab = "mse", 
     cex.main = 2.5, cex.lab = 2.5, cex.axis = 2, type = "n")

# add points and lines
for(jj in 1:ncol(plot4_mat))
{
  points(x = plot4Ks, y = plot4_mat[,jj], col = cols[jj], pch = 20)
  lines(x = plot4Ks, y = plot4_mat[,jj], col = cols[jj], lwd = 1.5)
}
dev.off()

## Plot 5: Mean squared error as a function of r.

# Parameters
plot5_thetas <- thetas[[length(ns)]][[1]]
plot5n <- ns[length(ns)]
plot5rs <- plot5_thetas$r[plot5_thetas$K == 1]
plot5Ks <- c(3,5,7,9,11,13,15)

# Plot data
plot5_data <- list()
for(K in plot5Ks)
{
  plot5_data[[as.character(K)]] <- data.frame(y = rowMeans(subset(mse[[which(ns == plot5n)]][[1]],
                                                              thetas[[length(ns)]][[1]]$K %in% K &
                                                              thetas[[length(ns)]][[1]]$r %in% plot5rs)))
}
plot5_df <- bind_rows(plot5_data,.id = "K")

# Plotting parameters
xlims <- c(min(plot5rs),max(plot5rs))
ylims <- c(min(plot5_df$y),max(plot5_df$y))
cols <- viridis(length(plot5Ks))

# shadow plot
pdf(file.path(plot_directory,"mse_by_radius.pdf"))
plot(x = plot5rs, xlim = xlims, ylim = ylims, xlab = "radius", ylab = "mse", 
     cex.main = 2.5, cex.lab = 2.5, cex.axis = 2, type = "n")
for(ii in 1:length(plot5Ks)){
  K <- plot5Ks[ii]
  points(x = plot5rs, y = plot5_df$y[plot5_df$K == K], col = cols[ii], pch = 20)
  lines(x = plot5rs, y = plot5_df$y[plot5_df$K == K], col = cols[ii], lwd = 1.5)
}
legend_text <- unique(plot5Ks)
legend("topright", legend = legend_text, col = cols, pch = 20,
       bg = "white", inset = .01, cex = 1.75)
dev.off()
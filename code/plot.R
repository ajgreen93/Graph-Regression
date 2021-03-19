library(dplyr)
library(ggplot2)
library(gridExtra)

### Plot. ###
save_directory <- "data/20210318210219" # Please change this to whichever directory you prefer.
plot_directory <- file.path(save_directory,"plots")
dir.create(plot_directory)

# Load data.
load(file.path(save_directory,"configs.R"))
d <- configs$d; ns <- configs$ns; methods <- configs$methods
load(file.path(save_directory,"best_fits_by_method.R"))
load(file.path(save_directory,"thetas.R"))
load(file.path(save_directory,"mse.R"))
load(file.path(save_directory,"Xs.R"))
load(file.path(save_directory,"f0s.R"))
load(file.path(save_directory,"Ys.R"))

# Plot of mse by tuning parameter.
pdf(file.path(plot_directory,"mse.pdf"))
for(jj in 1:length(methods))
{
  tuning_plots <- vector(mode = "list",length = length(ns))
  for(ii in 1:length(ns))
  {
    mse_ii_jj <- mse[[ii]][[jj]]
    thetas_ii_jj <- thetas[[ii]][[jj]]
    mean_mse <- rowMeans(mse_ii_jj)
    sd_mse <- apply(mse_ii_jj,1,sd)
    
    if(names(methods)[jj] == "laplacian_smoothing" || is.null(methods))
    {
      plot_df <- data.frame(x = log(thetas_ii_jj[,2]),y = mean_mse,
                            col = as.factor(round(thetas_ii_jj[,1],3)))
      tuning_plots[[ii]] <- 
        ggplot(data = plot_df,aes(x = x,y = y, color = col)) + geom_line() + 
        labs(x = "log(rho)",
             y = "mse",
             title = paste0("n = ",ns[ii],"")) +
        theme_bw()
    } else if(names(methods)[jj] == "knn")
    {
      plot_df <- data.frame(x = thetas_ii_jj[,1], y = mean_mse)
      tuning_plots[[ii]] <- 
        ggplot(data = plot_df,aes(x = x,y = y)) + geom_line() + 
        labs(x = "k",
             y = "mse",
             title = paste0("n = ",ns[ii],"")) +
        theme_bw()
    }
   
  }
  ncol_plots <- floor(sqrt(length(ns)))
  do.call("grid.arrange", c(tuning_plots, ncol=ncol_plots))
}
dev.off()

# Plot of true regression function

# Parameters
xlim = c(-1,1)
ylim = c(min(unlist(Ys)),max(unlist(Ys)))
col = "grey34"
cex = .6
lwd = 1.1
cex.main = 2.5
cex.lab = 2.5
cex.axis = 2

if(d == 1)
{
  for(ii in 1:length(ns))
  {
    pdf(file.path(plot_directory,paste0("regression_function_",ii,".pdf")))
    plot_df <- data.frame(x = Xs[[ii]],f0 = f0s[[ii]], y = Ys[[ii]])
    plot(x = plot_df$x,y = plot_df$y, 
         col = col, 
         cex = cex, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
         xlim = xlim, ylim = ylim, lwd = lwd,
         xlab = "",ylab = "", main = "True function")
    lines(x = plot_df$x[order(plot_df$x)],y = plot_df$f0[order(plot_df$x)],lwd = 1.5)
    dev.off()
  }
}

# Plot of fitted values for best choice of tuning parameter.
if(d == 1)
{
  for(jj in 1:length(methods))
  {
    for(ii in 1:length(ns))
    {
      plot_name <- paste0(names(methods)[jj],"_estimate_",ii,".pdf")
      pdf(file.path(plot_directory,plot_name))
      plot_df <- data.frame(x = Xs[[ii]], y = Ys[[ii]],fhat = best_fits_by_method[[ii]][[jj]])
      plot(x = plot_df$x,y = plot_df$y, 
           col = col,
           cex = cex, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
           xlim = xlim, ylim = ylim, lwd = lwd,
           xlab = "",ylab = "", main = "Laplacian smoothing")
      lines(x = plot_df$x[order(plot_df$x)],y = plot_df$fhat[order(plot_df$x)],lwd = 1.5,
            col = "blue")
      dev.off()
    }
  }
}

# Plot of mse---for best choice of tuning parameter---by sample size
pdf(file.path(plot_directory,"mse_by_sample_size.pdf"))
par(mar = c(5.1,4.5,4.1,4.1)) # The left hand side was getting cut off a tad.
plot_dfs_best_mse <- vector(mode = "list", length = length(methods))
names(plot_dfs_best_mse) <- names(methods)
fitted_slopes <- numeric(length(methods))
names(fitted_slopes) <- names(methods)
for(jj in 1:length(methods))
{
  if(names(methods)[jj] == "knn") next
  best_mse <- numeric()
  minimax_mse <- numeric()
  for(ii in 1:length(ns))
  {
    mse_ii_jj <- mse[[ii]][[jj]]
    best_mse[ii] <- min(rowMeans(mse_ii_jj))
    minimax_mse[ii] <- ns[ii]^{-2/(2 + d)}
  }
  
  # Rescale minimax mse to match intercept with best_mse
  minimax_mse <- minimax_mse * (best_mse[1]/minimax_mse[1])
  plot_dfs_best_mse[[jj]] <- data.frame(x = ns, y = best_mse, z = minimax_mse) 
  
  # fitted slope
  log_best_mse <- log(best_mse)
  log_ns <- log(ns)
  fitted_slopes[jj] <- round( lm(log_best_mse ~ log_ns)$coefficients[2], 2)
  
  # hack to change names for plotting
  if(names(plot_dfs_best_mse)[[jj]] == "laplacian_smoothing") names(plot_dfs_best_mse)[[jj]] <- "LS"
  names(plot_dfs_best_mse)[[jj]] <- paste0(names(plot_dfs_best_mse)[[jj]],
                                           " [Slope = ", fitted_slopes[jj],"].")
}

title <- paste0("d = ", d,". Minimax slope = ", -2,"/", d+2, ".")
plot_df_best_mse <- bind_rows(plot_dfs_best_mse, .id = "method")
legend_text <- unique(plot_df_best_mse$method)


# Plotting parameters
xlims <- c(min(ns),max(ns))
ylims <- c(.025, 1.2)

# only two methods
cols <- c("red","blue") 
pchs <- c(21,23)
ltys <- c(2,3)

# shadow plot
plot(x = ns, xlim = xlims, ylim = ylims,
     log = "xy", xlab = "Sample size", ylab = "Mean squared error", main = title, cex.main = 2, cex.lab = 2, cex.axis = 1.5)

# add points and lines
for(jj in 1:length(methods))
{
  if(names(methods)[jj] == "knn") next
  points(x = ns, y = plot_dfs_best_mse[[jj]]$y, bg = cols[jj], pch = pchs[jj], cex = 1.5)
  lines(x = ns, y = plot_dfs_best_mse[[jj]]$y, col = cols[jj], lty = ltys[jj], lwd = 1.5)
}
lines(x = ns, y = minimax_mse)
grid(equilogs = F, lwd = 2)
legend("bottomleft", legend = legend_text, col = cols[2], pch = pchs[2], lty = ltys[2], pt.bg = cols[2],
       bg = "white", inset = .01, cex = 1.75)
dev.off()

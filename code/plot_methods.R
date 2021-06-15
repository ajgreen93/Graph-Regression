plot_fxn <- function(x,y,f,title = NULL,cols = "black",domains = c(0,1)){
    # Observed values.
    if(!is.list(domains)) domains <- list(domains)
    domain <- c(min(unlist(domains)), max(unlist(domains)))
    plot(x = x,y = y, 
         col = "grey34",cex = .6,
         cex.main = 2.5, cex.lab = 2.5, cex.axis = 2,
         lwd = 1.1,xlab = "",ylab = "", main = title)
  
    # Plot functions.
    if(!is.list(f))
    {
      f <- list(f)
    }
    for(jj in 1:length(f))
    {
      for(interval in domains)
      {
        interval_indx <- c(x >= interval[1] & x <= interval[2])
        plot_x <- x[interval_indx]
        plot_f <- f[[jj]][interval_indx]
        lines(x = plot_x[order(plot_x)], y = plot_f[order(plot_x)],lwd = 2,
              col = cols[jj])
      }
    }
}


# Plot of mse---for best choice of tuning parameter---by sample size
plot_best_mse <- function(methods,mse,sd = T,validate = validate_mse,cols = NULL,
                          legend = T, title = NULL, rate = T){
  if(is.null(cols)){
    stopifnot(length(methods) <= 3)
    cols <- c("red","blue","green")[1:length(methods)]     
  }
  
  plot_dfs_best_mse <- vector(mode = "list", length = length(methods))
  names(plot_dfs_best_mse) <- names(methods)
  fitted_slopes <- numeric(length(methods))
  names(fitted_slopes) <- names(methods)
  for(jj in 1:length(methods))
  {
    method <- methods[[jj]]
    best_mse <- numeric()
    minimax_mse <- numeric()
    sd_best_mse <- numeric()
    for(ii in 1:length(ns))
    {
      mse_ii_jj <- mse[[ii]][[jj]]
      best_mse[ii] <- min(rowMeans(mse_ii_jj),na.rm = T)
      sd_best_mse[ii] <- apply(mse_ii_jj,1,sd)[which.min(rowMeans(mse_ii_jj))]/sqrt(ncol(mse_ii_jj))
      minimax_mse[ii] <- ns[ii]^{-2*s/(2*s + d)}
    }
    # Rescale minimax mse to match intercept with best_mse
    minimax_mse <- minimax_mse * (best_mse[1]/minimax_mse[1])
    plot_dfs_best_mse[[jj]] <- data.frame(x = ns, y = best_mse,sd = sd_best_mse) 
    
    # fitted slope
    log_best_mse <- log(best_mse)
    log_ns <- log(ns)
    fitted_slopes[jj] <- round( lm(log_best_mse ~ log_ns)$coefficients[2], 2)
    

    # hack to change names for plotting
    name <- names(plot_dfs_best_mse)[[jj]]
    names(plot_dfs_best_mse)[[jj]] <- case_when(
      name == "laplacian_smoothing" ~ "LS",
      name == "laplacian_eigenmaps" ~ "LE",
      name == "spectral_projection" ~ "SP",
      name == "least_squares"       ~ "LS",
      name == "laplacian_eigenmaps_plus_kernel_smoothing" ~ "LE+KS"
    ) 
    names(plot_dfs_best_mse)[[jj]] <- paste0(names(plot_dfs_best_mse)[[jj]],
                                             " [Slope = ", fitted_slopes[jj],"].")
  }
  
  if(is.null(title)) title <- paste0("d = ", d,", s = ",s,".", "Minimax slope = ", -2*s,"/", d+2*s, ".")
  plot_df_best_mse <- bind_rows(plot_dfs_best_mse, .id = "method") 
  legend_text <- unique(plot_df_best_mse$method)
  
  # Plotting parameters
  xlims <- c(min(ns),max(ns))
  ylims <- c(min(plot_df_best_mse$y - plot_df_best_mse$sd,minimax_mse)     , 
             max(plot_df_best_mse$y + plot_df_best_mse$sd,minimax_mse))
  
  # shadow plot
  plot(x = ns, xlim = xlims, ylim = ylims,
       log = "xy", xlab = "Sample size", ylab = "Mean squared error", 
       main = title,
       cex.main = 2.5, cex.lab = 2.5, cex.axis = 2)
  
  # add points and lines
  for(jj in 1:length(methods))
  {
    points(x = ns, y = plot_dfs_best_mse[[jj]]$y, col = cols[jj], pch = 20)
    lines(x = ns, y = plot_dfs_best_mse[[jj]]$y, col = cols[jj],lwd = 1.5)
    if(sd)
    {
      lines(x = ns, y = plot_dfs_best_mse[[jj]]$y + plot_dfs_best_mse[[jj]]$sd, 
            col = cols[jj],
            lwd = 1.5,
            lty = 2)
      lines(x = ns, y = plot_dfs_best_mse[[jj]]$y - plot_dfs_best_mse[[jj]]$sd, 
            col = cols[jj],
            lwd = 1.5,
            lty = 2)
    }
  }
  
  # Complete the plot
  if(rate) lines(x = ns, y = minimax_mse)
  grid(equilogs = F, lwd = 2)
  if(legend){
    legend("bottomleft", legend = legend_text, col = cols, pch = 20,
         bg = "white", inset = .01, cex = 1.75)
  }
}

# Plot the critical radius of tests.
# Inputs:
# -- methods: list, as in configs/testing/eigenfunction_1d_1s.R
# -- power: list, as computed by testing_error_pipeline.R
# -- type_II_error: the maximum tolerated level of type II error.
# -- sd: plot standard error bars?
# -- cols: a vector giving colors for different methods
plot_testing_critical_radius <- function(methods,err,type_II_error = .5,
                                         sd = T, cols = NULL, verbose = T){
  # Calculate critical radius.
  critical_radius <- find_critical_radius(err,methods,make_f0s,type_II_error,sd = T,
                                          verbose)
  
  # Calculate slope
  slope <- numeric(length(critical_radius))
  names(slope) <- names(methods)
  for(jj in 1:length(slope))
  {
    slope[jj] <- mse_slope(critical_radius[[jj]][,"crit_radius"],ns)
  }
  
  # Calculate minimax rate
  minimax_crit_radius <- numeric(length(ns))
  for(ii in 1:length(ns))
  {
    minimax_crit_radius[ii] <- ns[ii]^(-4*s/(4*s + d))
  }
  minimax_crit_radius <- minimax_crit_radius * (critical_radius[[1]][1,1] / minimax_crit_radius[1])
  
  # Data frame for plotting
  plot_df_list <- list()
  for(jj in 1:length(critical_radius))
  {
    plot_df_list[[jj]] <- data.frame(critical_radius[[jj]])
    name <- names(methods)[jj]
    abbr <- case_when(
      name == "laplacian_smoothing" ~ "LS",
      name == "laplacian_eigenmaps" ~ "LE",
      name == "spectral_projection" ~ "SP",
      name == "least_squares"       ~ "LS",
      name == "laplacian_eigenmaps_plus_kernel_smoothing" ~ "LE+KS"
    ) 
    names(plot_df_list)[jj] <- paste0(abbr," [Slope = ", round(slope[jj],2),"].")
  }
  plot_df <- bind_rows(plot_df_list,.id = "method")
  
  title <- paste0("d = ", d,", s = ",s,".", "Minimax slope = ", -4*s,"/", d+4*s, ".")
  legend_text <- unique(plot_df$method)
  
  # Plotting parameters
  xlims <- c(min(ns),max(ns))
  ylims <- c(min(plot_df[,-1],minimax_crit_radius), max(plot_df[-1],minimax_crit_radius))
  
  # shadow plot
  plot(x = ns, xlim = xlims, ylim = ylims,
       log = "xy", xlab = "Sample size", ylab = "Critical radius", 
       main = title,
       cex.main = 2.5, cex.lab = 2.5, cex.axis = 2)
  
  # add points and lines
  for(jj in 1:length(methods))
  {
    points(x = ns, y = plot_df_list[[jj]]$crit_radius, col = cols[jj], pch = 20)
    lines(x = ns, y = plot_df_list[[jj]]$crit_radius, col = cols[jj],lwd = 1.5)
    if(sd)
    {
      lines(x = ns, y = plot_df_list[[jj]]$crit_radius_lb, 
            col = cols[jj],
            lwd = 1.5,
            lty = 2)
      lines(x = ns, y = plot_df_list[[jj]]$crit_radius_ub, 
            col = cols[jj],
            lwd = 1.5,
            lty = 2)
    }
  }
  
  # Complete the plot
  lines(x = ns, y = minimax_crit_radius)
  grid(equilogs = F, lwd = 2)
  legend("bottomleft", legend = legend_text, col = cols, pch = 20,
         bg = "white", inset = .01, cex = 1.75)
}

# Plot the best mse for each value of a tuning parameter.
# Inputs:
# -- parameter_name: name of the parameter
plot_mse_by_tuning_parameter <- function(parameter_name,
                                         methods,thetas,mse,
                                         cols = rep("black",length(methods)),
                                         title = "",
                                         ylims = NULL)
{
  for(ii in 1:length(mse))
  {
    # Collect plotting data.
    plot_data <- list(length(methods))
    for(jj in 1:length(methods))
    {
      parameter_values <- unique(thetas[[ii]][[jj]][[parameter_name]])
      plot_data[[jj]] <- numeric(length(parameter_values))
      names(plot_data[[jj]]) <- parameter_values
      for(kk in 1:length(parameter_values))
      {
        plot_data[[jj]][kk] <- mse[[ii]][[jj]] %>% 
          subset(thetas[[ii]][[jj]][[parameter_name]] == parameter_values[kk]) %>% rowMeans() %>% min()
      }
    }
    
    # Plot
    xlims <- c(unlist(plot_data) %>% names() %>% as.numeric() %>% min(), 
               unlist(plot_data) %>% names() %>% as.numeric() %>% max())
    if(is.null(ylims)) ylims <- c(unlist(plot_data) %>% min(), unlist(plot_data) %>% max())
    
    plot(x = NULL, xlim = xlims, ylim = ylims, xlab = parameter_name, ylab = "mse", 
         main = title[ii], type = "n")
    for(jj in 1:length(plot_data))
    {
      points(x = as.numeric(names(plot_data[[jj]])), y =  plot_data[[jj]], col = cols[jj], pch = 20)
      lines(x = as.numeric(names(plot_data[[jj]])), y =  plot_data[[jj]], col = cols[jj], lwd = 1.5)
    }
    grid(lwd = 2)
  }
}

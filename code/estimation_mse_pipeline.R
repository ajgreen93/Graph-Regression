library(RANN)
#----------------------------------------------------#
# This is a pipeline for running regression estimation experiments.
# It reads in:
# -- a function sample_X, for sampling design points
# -- a function f_0, for sampling responses
#----------------------------------------------------#

source("estimation_mse_config.R")

mse <- vector(mode = "list",length = length(ns))
for(ii in 1:length(ns))
{
  n <- ns[ii]
  
  # Initialize tuning parameters (thetas).
  # Note: I have granted methods the right to access the distribution P from which X is sampled
  # in initalizing hyper-parameters.
  thetas <- vector(mode = "list", length(methods))
  for(jj in 1:length(methods))
  {
    method <- methods[[jj]]
    thetas[[jj]] <- initialize_thetas[[jj]](sample_X, n)
  }
  
  for(iter in iters)
  {
    ### Data. ###
    f0_evaluations <- apply(X,1,f_0)
    Y <- f0_evaluations + rnorm(n,0,1)
    
    ### Analysis. ##
    for(method in methods)
    {
      for(theta in thetas[[method]]){
        estimator <- method(theta)
        fhat[[theta]] <- estimator(Y,X)
      }
    }
    
    ### Evaluation. ###
    for(theta in theta_choices){
      mse[[n]][[iter]][[method]][theta] <- mean( (Y - fhat[[theta]])^2 )
    }
  }
}


### Save. ###
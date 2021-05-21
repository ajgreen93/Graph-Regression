library(RANN)
library(dplyr)
library(Matrix)
library(reshape2)
library(RSpectra)
source("sample.R")
source("graph.R")
source("estimators.R")
source("manifold.R")
source("misc.R")

# Please change this line to whichever config you wish to run.
source("configs/laplacian_eigenmaps/mse/sobolev_2d_2s.R")

# Options for running the pipeline
verbose <- F
save_fits <- F

#----------------------------------------------------#
# This is the pipeline for running regression estimation experiments.
#----------------------------------------------------#

# For output.
best_fits_by_method <- vector(mode = "list",length = length(ns))
mse <- vector(mode = "list",length = length(ns))
Xs <- vector(mode = "list", length = length(ns))
f0s <- vector(mode = "list", length = length(ns))
Ys <- vector(mode = "list", length = length(ns))
thetas <- vector(mode = "list", length = length(ns))
if(save_fits) fits <- vector(mode = "list", length = length(ns))

for(ii in 1:length(ns))
{
  n <- ns[ii]
  f0 <- make_f0(d,n)
  mse[[ii]] <- vector(mode = "list", length = length(methods))
  if(save_fits)
  {
    fits[[ii]] <- vector(mode = "list",length = length(methods))
  }

  # Initialize tuning parameters (thetas).
  # Note: I have granted methods the right to access the distribution P from which X is sampled
  # in initalizing hyper-parameters.
  thetas[[ii]] <- vector(mode = "list", length(methods))
  for(jj in 1:length(methods))
  {
    thetas[[ii]][[jj]] <- initialize_thetas[[jj]](sample_X, n)
    mse[[ii]][[jj]] <- matrix(nrow = nrow(thetas[[ii]][[jj]]), ncol = iters)
    if(save_fits) fits[[ii]][[jj]] <- array(dim = c(nrow(thetas[[ii]][[jj]]),n,iters))
  }
  for(iter in 1:iters)
  {
    fits_by_method <- vector(mode = "list",length = length(methods))
    
    ### Data. ###
    X <- sample_X(n)
    f0_evaluations <- apply(X,1,f0)
    Y <- f0_evaluations + rnorm(n,0,1)
    
    ### Analysis. ##
    for(jj in 1:length(methods))
    {
      method <- methods[[jj]]
      fits_method <- matrix(nrow = nrow(thetas[[ii]][[jj]]), ncol = n)
      
      # Do any computations here that we would have to repeat for each hyperparameter (theta)
      # E.g. compute eigenvectors.
      precompute_data <- attributes(method)$precompute_data
      if(!is.null(precompute_data)){
        precomputed_data <- precompute_data(Y,X, thetas[[ii]][[jj]])
      } else{
        precomputed_data <- NULL
      }
      for(kk in 1:nrow(thetas[[ii]][[jj]])){
        theta <- slice(thetas[[ii]][[jj]],kk)
        estimator <- method(theta)
        fits_method[kk,] <- estimator(Y,X)
        if(verbose){
          logger::log_info("n: ", n, ".",
                         "Iter: ", iter," out of ", iters, ".",
                         "Method: ", jj, " out of ", length(methods), ".",
                         "Parameter: ", kk, " out of ", nrow(thetas[[ii]][[jj]]))
        }
      }
      fits_by_method[[jj]] <- fits_method
      logger::log_info("n: ", n, ".",
                       "Iter: ", iter," out of ", iters, ".",
                       "Method: ", jj, " out of ", length(methods), ".")
    }
    
    ### Evaluation. ###
    for(jj in 1:length(methods)){
      for(kk in 1:nrow(thetas[[ii]][[jj]])){
        mse[[ii]][[jj]][kk,iter] <- mean( (f0_evaluations - fits_by_method[[jj]][kk,])^2 )
      }
    }
    
    ### Save fits.
    if(save_fits) fits[[ii]][[jj]][,,iter] <- fits_by_method[[jj]]
  }
  
  ### Save best fits. ###
  for(jj in 1:length(methods))
  {
    mse_jj <- mse[[ii]][[jj]]
    mean_mse <- rowMeans(mse_jj)
    which_min_mse <- which.min(mean_mse)
    
    # Taking the last iteration fit, it doesn't matter.
    best_fits_by_method[[ii]][[jj]] <- fits_by_method[[jj]][which_min_mse,]
  }
  
  # Taking the last iteration X, it doesn't matter. 
  Xs[[ii]] <- X
  f0s[[ii]] <- f0_evaluations
  Ys[[ii]] <- Y
}

### Save. ###
save_directory <- file.path("data",gsub("[^[:alnum:]]", "", Sys.time()))
for(directory in c(save_directory))
{
  dir.create(directory)
}
configs <- list(d = d,
                s = s,
                ns = ns,
                methods = methods,
                sample_X = sample_X,
                make_f0 = make_f0)
save(configs, file = paste0(save_directory, '/configs.R'))
save(best_fits_by_method, file = paste0(save_directory, '/best_fits_by_method.R'))
save(thetas, file = paste0(save_directory, '/thetas.R'))
save(mse, file = paste0(save_directory, '/mse.R'))
save(Xs, file = paste0(save_directory, '/Xs.R'))
save(Ys, file = paste0(save_directory, '/Ys.R'))
save(f0s, file = paste0(save_directory, '/f0s.R'))
if(save_fits) save(fits, file = paste0(save_directory, '/fits.R'))

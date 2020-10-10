#------------------------------------------------#
# Configs for first estimation MSE experiment.
#------------------------------------------------#

### General configs.
d <- 1
ns <- round( 4 * 10^{seq(2,3,length.out = 10)} )
iters <- 5

### Configs for sampling. ###
sample_X <- function(n){
  X <- matrix(runif(n,-1,1),ncol = 1)
  return(X)
}

make_f_0 <- function(n)
{
  # Choose M as small as possible such that 
  # minimax rate of f_0 is n^{-2/(2 + d)}
  m <- round(1/2 * n^{1/d * (2 - (4 + d)/(2 + d))})
  f0 <- function(x){
    return(1/m * cos(pi * m * x))
  }
}


### Methods. ###
methods <- list(make_laplacian_smoothing)
initialize_thetas <- list(initialize_laplacian_smoothing_thetas)
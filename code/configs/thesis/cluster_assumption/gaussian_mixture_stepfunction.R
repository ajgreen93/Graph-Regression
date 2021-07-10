#------------------------------------------------#
# Configs for Cluster Assumption plot in thesis, when 
# design density is a mixture of two Gaussians.
#------------------------------------------------#

### General configs.
d <- 1
ns <- round(seq(400,2000,length.out = 10))
iters <- 1

### Configs for sampling. ###
sample_X <- function(n){
  # Parameters for Gaussians
  mu <- c(-1,1)
  sigma <- c(.4,.4)
  pi <- c(.5,.5)
  
  # Sample
  sample_fxn <- make_sample_gaussian_mixture(d, mu = mu, sigma = sigma, pi = pi)
  sample_fxn(n)
}
make_f0 <- function(d,n){
  height <- 1
  f0 <- make_stepfunction(height)
}


### Methods. ###
methods <- list(
  laplacian_eigenmaps = make_laplacian_eigenmaps,
  laplacian_smoothing = make_laplacian_smoothing,
  least_squares = make_least_squares,
  kernel_smoothing = make_kernel_smoothing
)

### Tuning parameters for methods. ###
# NOTE: to minimize computation time, we choose max_K_graph smaller than max_K_least_squares
#       This value of max_K_graph will be more than big enough for the stepfunction with cluster assumption.
max_K_graph <- 10
max_K_least_squares <- 30
max_degree_kernel_smoothing <- function(n){round(log(n)^2)}
dist_to_r <- function(dists){quantile(dists,.75)}
initialize_thetas <- list(
  laplacian_eigenmaps = function(sample_X,n){
    initialize_laplacian_eigenmaps_thetas(sample_X,n,max_K = max_K_graph,dist_to_r = dist_to_r)
  },
  laplacian_smoothing = function(sample_X,n){
    initialize_laplacian_smoothing_thetas(sample_X,n,max_K = max_K_graph,dist_to_r = dist_to_r)
  },
  least_squares = function(sample_X,n){
    initialize_spectral_projection_thetas(sample_X,n,max_K = max_K_least_squares)
  },
  kernel_smoothing = function(sample_X,n){
    initialize_kernel_smoothing_thetas(sample_X,n,max_degree = max_degree_kernel_smoothing(n))
  }
)
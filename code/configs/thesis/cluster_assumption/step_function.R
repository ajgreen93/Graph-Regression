#------------------------------------------------#
# Configs for Cluster Assumption plot in thesis.
#------------------------------------------------#

### General configs.
d <- 1
s <- 1
M <- 2^(s - d)
ns <- round(seq(400,2000,length.out = 10))
iters <- 5

### Configs for sampling. ###
sample_X <- function(n){
  sep <- (log(n))^2/(2 * n)
  make_sample_uniform_mixture(d, sep = sep)(n)
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
initialize_thetas <- list(
  initialize_laplacian_eigenmaps_thetas,
  initialize_laplacian_smoothing_thetas,
  initialize_spectral_projection_thetas,
  initialize_kernel_smoothing_thetas
  )
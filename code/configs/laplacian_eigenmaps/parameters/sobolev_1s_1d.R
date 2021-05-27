#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 1.
#------------------------------------------------#

### General configs.
d <- 1
s <- 1
M <- 4
ns <- 1000
iters <- 100

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- function(d,n){
  g0 <- make_spectral_sobolev(d,n,s,M)
  f0 <- function(x){g0(x)}
}


### Methods. ###
methods <- list(
  laplacian_eigenmaps = make_laplacian_eigenmaps,
  laplacian_eigenmaps_plus_kernel_smoothing = make_laplacian_eigenmaps_plus_kernel_smoothing,
  spectral_projection = make_spectral_projection,
  least_squares = make_least_squares
)
initialize_thetas <- list(
  initialize_laplacian_eigenmaps_thetas,
  initialize_laplacian_eigenmaps_plus_kernel_smoothing_thetas,
  initialize_spectral_projection_thetas,
  initialize_spectral_projection_thetas
)
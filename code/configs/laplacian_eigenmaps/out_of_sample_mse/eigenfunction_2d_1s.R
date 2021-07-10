#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 2.
#------------------------------------------------#

### General configs.
d <- 2
s <- 1
M <- 2^(s - 2*d + 1)
ns <- round(seq(1000,4000,length.out = 10))
iters <- 200

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- function(d,n){
  g0 <- make_eigenfunction(d,round((M^2*n)^(d/(2*s + d))),M,s)
  f0 <- function(x){2 * g0(x)}
}

### Methods. ###
methods <- list(
  laplacian_eigenmaps_plus_kernel_smoothing = make_laplacian_eigenmaps_plus_kernel_smoothing,
  spectral_projection = make_spectral_projection,
  least_squares = make_least_squares
)
initialize_thetas <- list(
  initialize_laplacian_eigenmaps_plus_kernel_smoothing_thetas,
  initialize_spectral_projection_thetas,
  initialize_spectral_projection_thetas
)
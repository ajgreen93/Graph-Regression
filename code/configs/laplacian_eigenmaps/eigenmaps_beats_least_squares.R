#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 2.
#------------------------------------------------#

### General configs.
d <- 1
ns <- round(seq(400,2000,length.out = 10))
iters <- 50

### Configs for sampling. ###
sample_X <- make_sample_gaussian_mixture(d,mu = c(-1,1),pi = c(3/4, 1/4), sigma = c(.4,.4))
# sample_X <- make_sample_uniform(d)
make_f0 <- make_doppler


### Methods. ###
methods <- list(laplacian_eigenmaps = make_laplacian_eigenmaps,
                least_squares = make_least_squares,
                kernel_smoothing = make_kernel_smoothing)
initialize_thetas <- list(initialize_laplacian_eigenmaps_thetas,
                          initialize_spectral_projection_thetas,
                          initialize_kernel_smoothing_thetas)
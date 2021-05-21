#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 1.
#------------------------------------------------#

### General configs.
d <- 2
ns <- 1000
# ns <- round(seq(400,2000,length.out = 10))
iters <- 50

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- make_spectral_sobolev


### Methods. ###
methods <- list(laplacian_eigenmaps = make_laplacian_eigenmaps,
                spectral_projection = make_spectral_projection)
initialize_thetas <- list(initialize_laplacian_eigenmaps_thetas,
                          initialize_spectral_projection_thetas)
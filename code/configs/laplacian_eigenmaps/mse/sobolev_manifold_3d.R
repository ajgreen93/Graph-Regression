#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 1.
#------------------------------------------------#

### General configs.
D <- 10
d <- 3
# ns <- 2000
ns <- round(seq(400,2000,length.out = 10))
iters <- 50

### Configs for sampling. ###
sample_X <- make_sample_manifold(d,D)
make_f0 <- make_spectral_sobolev_manifold


### Methods. ###
methods <- list(laplacian_eigenmaps = make_laplacian_eigenmaps)
initialize_thetas <- list(initialize_laplacian_eigenmaps_thetas)
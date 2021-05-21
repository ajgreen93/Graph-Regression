#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 1.
#------------------------------------------------#

### General configs.
d <- 2
s <- 2
M <- 16
# ns <- 2000
ns <- round(seq(400,2000,length.out = 10))
iters <- 50

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- function(d,n){make_spectral_sobolev(d,max(ns),s,M)}


### Methods. ###
methods <- list(laplacian_eigenmaps = make_laplacian_eigenmaps,
                spectral_projection = make_spectral_projection,
                least_squares = make_least_squares)
initialize_thetas <- list(initialize_laplacian_eigenmaps_thetas,
                          initialize_spectral_projection_thetas,
                          initialize_spectral_projection_thetas)
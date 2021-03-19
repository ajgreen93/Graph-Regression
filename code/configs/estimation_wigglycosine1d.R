#------------------------------------------------#
# Configs for first estimation MSE experiment.
#------------------------------------------------#

### General configs.
d <- 1
ns <- 1000
iters <- 1

### Configs for sampling. ###
sample_X <- make_sample_gaussian(d)
make_f0 <- make_wiggly_gaussian


### Methods. ###
methods <- list(laplacian_smoothing = make_laplacian_smoothing_knn)
initialize_thetas <- list(initialize_laplacian_smoothing_knn_thetas)
#------------------------------------------------#
# Configs for second estimation MSE experiment.
#------------------------------------------------#

### General configs. ###
d <- 2
ns <- round( 10^{seq(2,3,length.out = 4)} )
iters <- 5

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- make_cosine_f0


### Methods. ###
methods <- list(make_laplacian_smoothing_knn)
initialize_thetas <- list(initialize_laplacian_smoothing_knn_thetas)

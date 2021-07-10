#------------------------------------------------#
# Configs for second estimation MSE experiment.
#------------------------------------------------#

### General configs. ###
d <- 4
ns <- round( 10^seq(2.5,3.5,length.out = 10) )
iters <- 5

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- make_linear_f0


### Methods. ###
methods <- list(
  knn = make_knn,
  laplacian_smoothing = make_laplacian_smoothing_knn
)

initialize_thetas <- list(
  knn = initialize_knn_thetas,
  laplacian_smoothing = initialize_laplacian_smoothing_knn_thetas
)
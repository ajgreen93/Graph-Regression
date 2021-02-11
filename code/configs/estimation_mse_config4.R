#------------------------------------------------#
# Configs for fourth estimation MSE experiment.
#------------------------------------------------#

### General configs. ###
d <- 3
ns <- round( 10^seq(3,4,length.out = 10) )
iters <- 5

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- make_cosine_f0


### Methods. ###
methods <- list(
  knn = make_knn,
  laplacian_smoothing = make_laplacian_smoothing_knn
)

initialize_thetas <- list(
  knn = initialize_knn_thetas,
  laplacian_smoothing = initialize_laplacian_smoothing_knn_thetas
)
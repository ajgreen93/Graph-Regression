#------------------------------------------------#
# Configs for second estimation MSE experiment.
#------------------------------------------------#

### General configs. ###
d <- 1
ns <- seq(200,5000,length.out = 10)
iters <- 5

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
make_f0 <- make_cosine_f0


### Methods. ###
methods <- list(
  knn = make_knn
 #make_laplacian_smoothing_knn
                )
initialize_thetas <- list(
  knn = initialize_knn_thetas
# initialize_laplacian_smoothing_knn_thetas
  )

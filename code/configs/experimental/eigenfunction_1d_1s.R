#------------------------------------------------#
# Configs for Laplacian Eigenmaps Figure 1.
#------------------------------------------------#

### General configs.
d <- 1
s <- 1
M <- 2^(s - d)
ns <- round(seq(1000,2000,length.out = 10))
iters <- 100

### Configs for sampling. ###
sample_X <- make_sample_uniform(d)
eigenfunction_indx <- round((M^2*ns[round(length(ns)/2)])^(d/(2*s + d)))
make_f0 <- function(d,n){
  g0 <- make_eigenfunction(d,eigenfunction_indx,M,s)
  f0 <- function(x){2 * g0(x)}
}


### Methods. ###
methods <- list(
  laplacian_eigenmaps = make_laplacian_eigenmaps,
  spectral_projection = make_spectral_projection,
  # least_squares = make_least_squares
  laplacian_smoothing = make_laplacian_smoothing,
  laplacian_smoothing = make_laplacian_smoothing
)

### Initialize parameters. ###
initialize_le <- function(sample_X,n){
  initialize_laplacian_eigenmaps_thetas(sample_X,n,
                                        n_rs = 15,
                                        max_degree = 40,
                                        max_K = 5*eigenfunction_indx,
                                        dist_to_r = function(x){quantile(x,.9)})}
initialize_sp <- function(sample_X,n){
  initialize_spectral_projection_thetas(sample_X,n,
                                        max_K = 5*eigenfunction_indx)
}
initialize_ls <- function(sample_X,n){
  initialize_laplacian_smoothing_thetas(sample_X,n,
                                        n_rs = 10,
                                        n_rhos = 20,
                                        max_degree = 15,
                                        max_K = 5*eigenfunction_indx,
                                        dist_to_r = function(x){quantile(x,.9)},
                                        rho_multiplier = 3)}
initialize_ls_oos <- function(sample_X,n){
  initialize_laplacian_smoothing_thetas(sample_X,n,
                                        n_rs = 10,
                                        n_rhos = 20,
                                        max_degree = 15,
                                        max_K = 5*eigenfunction_indx,
                                        dist_to_r = function(x){quantile(x,.9)},
                                        rho_multiplier = 3,
                                        labeled = F)}
initialize_thetas <- list(
  initialize_le,
  initialize_sp,
  # initialize_spectral_projection_thetas
  initialize_ls,
  initialize_ls_oos
)
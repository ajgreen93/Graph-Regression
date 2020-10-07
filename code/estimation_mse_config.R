#------------------------------------------------#
# Configs for first estimation MSE experiment.
#------------------------------------------------#

### General configs.
ns <- seq(100,1000,by = 100)
iters <- 5

### Configs for sampling. ###
sample_X <- function(n){
  X <- runif(n,0,1)
  return(X)
}

f0 <- function(x){
  
}

### Methods. ###
laplacian_smoothing <- function(Y,X,theta){
  
  # Build G_{n,r} over X.
  
  # Get L.
  
  # Solve (rho*L + I)f = Y for f.
  
}
initialize_laplacian_smoothing_thetas <- function(sample,n){
  # Choose a range of radii centered around r_connect, where
  # r_connect is the connectivity radius: i.e. the minimum radius such that 
  # observed minimum degree is log(n).
  thetas <- data.frame()
  r_connects <- numeric(10)
  for(iter in 1:10)
  {
    X <- sample(n)
    rneighborhood <- nn2(data = X, query = X, k = round(log(n)) + 1)
    r_connects[iter] <- max(rneighborhood$nn.dists[,round(log(n)) + 1])
  }
  r_s <- seq(mean(r_connects) - 3*sd(r_connects),mean(r_connects) + 3*sd(r_connects), length.out = 9)
  
  
  rho_s <- 
  return(thetas)
}
methods <- list(laplacian_smoothing)
initialize_thetas <- list(initialize_laplacian_smoothing_thetas)
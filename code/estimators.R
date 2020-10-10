make_laplacian_smoothing <- function(theta)
{
  r <- theta[["r"]]
  rho <- theta[["rho"]]
  laplacian_smoothing <- function(Y,X){
    
    # Build G_{n,r} over X.
    G <- neighborhood_graph(X,r)
    
    # Get L.
    L <- Laplacian(G)
    
    # Solve (rho*L + I)f = Y for f.
    f_hat <- as.numeric(solve(rho * L + diag(n),Y))
    
    return(f_hat)
  }
}

make_laplacian_smoothing_knn <- function(theta)
{
  k <- theta[["k"]]
  rho <- theta[["rho"]]
  laplacian_smoothing_knn <- function(Y,X){
    
    # Build G_{n,k} over X.
    G <- knn_graph(X,k)
    
    # Get L.
    L <- Laplacian(G)
    
    # Solve (rho*L + I)f = Y for f.
    f_hat <- as.numeric(solve(rho * L + diag(n),Y))
  }
}

initialize_laplacian_smoothing_thetas <- function(sample_X,n){
  thetas <- data.frame(r = numeric(), rho = numeric())
  
  # Choose a range of radii centered around r_connect, where
  # r_connect is the connectivity radius: i.e. the minimum radius such that 
  # observed minimum degree is 2.
  r_connects <- numeric(10)
  for(iter in 1:10)
  {
    X <- sample_X(n)
    rneighborhood <- nn2(data = X, query = X, k = 3)
    r_connects[iter] <- max(rneighborhood$nn.dists[,3]) # connectivity radius
  }
  # Choose a range of radii centered around r_connect.
  rs <- mean(r_connects)
  
  # Choose a range of penalty parameters.
  rho_optimal_theory <- function(r){1 / (n^{2/(2+d)} * n * r^{d + 2})} # bias-variance balance
  rho_optimals <- sapply(rs,rho_optimal_theory)
  rhos <- sapply(rho_optimals,FUN = function(rho){seq(rho, 10 * rho,length.out = 50)})
  
  # All combinations.
  for (ii in 1:length(rs))
  {
    for(jj in 1:nrow(rhos))
    {
      thetas <- thetas %>% add_row(r = rs[ii],rho = rhos[jj,ii])
    }
  }
  return(thetas)
}

initialize_laplacian_smoothing_knn_thetas <- function(sample_X,n){
  thetas <- data.frame(k = numeric(), rho = numeric())
  
  ks <- 5
  
  # Choose a range of penalty parameters.
  rho_optimal_theory <- function(k){n^{2/d} / (n^{2/(2+d)} * k^{(d + 2)/d})} # bias-variance balance
  rho_optimals <- sapply(ks,rho_optimal_theory)
  rhos <- sapply(rho_optimals,FUN = function(rho){seq(2^{d - 2} * rho, 5^d * rho,length.out = 50)})
  
  # All combinations.
  for (ii in 1:length(ks))
  {
    for(jj in 1:nrow(rhos))
    {
      thetas <- thetas %>% add_row(k = ks[ii],rho = rhos[jj,ii])
    }
  }
  return(thetas)
}
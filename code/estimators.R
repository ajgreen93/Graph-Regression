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

make_knn <- function(theta)
{
  k <- theta[["k"]]
  knn <- function(Y,X)
  {
    # Build G_{n,k} over X.
    G <- knn_graph(X,k)
    
    # normalize by row sums
    H <- G / rowSums(G)
    
    # kNN estimator
    f_hat <- as.numeric(H %*% Y)
    
    return(f_hat)
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
    rneighborhood <- nn2(data = X, query = X, k = round(log(n) + 1))
    r_connects[iter] <- max(rneighborhood$nn.dists[,round(log(n) + 1)]) # connectivity radius
  }
  # Choose a range of radii centered around r_connect.
  rs <- mean(r_connects) + sd(r_connects) * seq(-2,2,length.out = 10)
  
  
  # Choose a range of penalty parameters.
  rho_optimal_theory <- function(r){1 / (n^{2/(2+d)} * n * r^{d + 2})} # bias-variance balance
  rho_optimals <- sapply(rs,rho_optimal_theory)
  rhos <- sapply(rho_optimals,FUN = function(rho){seq(rho, 5 * rho,length.out = 25)})
  
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
  
  # Choose a range of connectivity parameters.
  ks <- seq(2, 17, by = 3)
  
  # Choose a range of penalty parameters.
  m <- choose_f0_frequency(d)
  a <- choose_f0_amplitude(d,m)
  rho_optimal_theory <- function(k){
    n^{2/d} / (n^{2/(2+d)} * k^{(d + 2)/d}) *
    (d * (a^2/2^d) * m^2 * pi^2)^{-2/(2+d)}
  } # bias-variance balance
  rho_optimals <- sapply(ks,rho_optimal_theory)
  rhos <- sapply(rho_optimals,FUN = function(rho){exp( seq(log(1/25 * rho),log(25 * rho),length.out = 25) )})
  
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

initialize_knn_thetas <- function(sample_X,n){
  m <- choose_f0_frequency(d)
  a <- choose_f0_amplitude(d,m)
  k_min <- 2
  k_max <- min(max(round(5 * (d * (a^2/2^d) * m^2 * pi^2)^{-d/(2+d)} * n^{2/(d + 2)}), 50), n/4)
  ks <- c(2:12, seq(13,k_max,by = 2))
  
  thetas <- data.frame(k = ks)
}
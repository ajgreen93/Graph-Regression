make_laplacian_eigenmaps <- function(theta)
{
  r <- theta[["r"]]
  K <- theta[["K"]]
  laplacian_eigenmaps <- function(Y,X){
    # If we've already computed the eigenvectors don't do it again
    if(exists("precomputed_data")){
      stopifnot(r %in% names(precomputed_data))
      L_eigenvectors <- precomputed_data[[which(names(precomputed_data) == r)]]
    } else{
      # Build G_{n,r} over X.
      G <- neighborhood_graph(X,r)
      
      # Get L.
      L <- Laplacian(G)
      
      # Get spectra.
      spectra <- get_spectra(L,K)
      L_eigenvectors <- spectra$vectors[,ncol(spectra$vectors):1]
    }
    V_K <- L_eigenvectors[,1:K] # Eigenvectors we need
    a_k <- t(V_K) %*% Y         # Empirical Fourier coefficients
    f_hat <- V_K %*% a_k        # Spectral projection
    
    return(f_hat)
  }
}
attr(make_laplacian_eigenmaps,"precompute_data") <- function(Y,X, theta_df){
  rs <- unique(theta_df$r)
  precomputed_data <- vector(mode = "list",length = length(rs))
  names(precomputed_data) <- rs
  for(kk in 1:length(rs)){
    r <- rs[kk]
    K <- thetas[[ii]][[jj]] %>% filter(r == !!r) %>% pull(K)
    
    # Build G_{n,r} over X.
    G <- neighborhood_graph(X,r)
    
    # Get L.
    L <- Laplacian(G)
    
    # Compute as many eigenvectors as we will need.
    spectra <- get_spectra(L,max(K))
    L_eigenvectors <- spectra$vectors[,ncol(spectra$vectors):1]
    
    precomputed_data[[kk]] <- L_eigenvectors
  }
  return(precomputed_data)
}

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

make_spectral_projection <- function(theta)
{
  K <- theta[["K"]]
  spectral_projection <- function(Y,X)
  {
    d <- ncol(X)
    n <- length(Y)
  
    if(exists("precomputed_data")){
      stopifnot(length(precomputed_data) == 1)
      psi_K <- precomputed_data[[1]][,1:K,drop = FALSE]
    } else{
      trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
      psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
    }
    a_k <- 1/n * (t(psi_K) %*% Y)   # Empirical Fourier coefficients
    f_hat <- psi_K %*% a_k        # Spectral projection
  }
} 
attr(make_spectral_projection,"precompute_data") <- function(Y,X, theta_df){
  precomputed_data <- vector(mode = "list",length = 1)
  
  K <- max(theta_df$K)
  trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
  psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
  
  precomputed_data[[1]] <- psi_K
  return(precomputed_data)
}

make_least_squares <- function(theta)
{
  K <- theta[["K"]]
  least_squares <- function(Y,X)
  {
    d <- ncol(X)
    n <- length(Y)
    
    if(exists("precomputed_data")){
      stopifnot(length(precomputed_data) == 1)
      psi_K <- precomputed_data[[1]][,1:K,drop = FALSE]
    } else{
      trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
      psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
    }
    
    f_hat <- lm.fit(x = psi_K,y = Y)$fitted.values # Least squares
  }
}
attr(make_least_squares,"precompute_data") <- function(Y,X, theta_df){
  precomputed_data <- vector(mode = "list",length = 1)
  
  K <- max(theta_df$K)
  trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
  psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
  
  precomputed_data[[1]] <- psi_K
  return(precomputed_data)
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

make_kernel_smoothing <- function(theta)
{
  r <- theta[["r"]]
  kernel_smoothing <- function(Y,X)
  {
    # Build G_{n,r} over X.
    G <- neighborhood_graph(X,r)
    
    # Add self loops
    diag(G) <- 1
    
    # normalize by row sums
    H <- G / rowSums(G)
    
    # kernel smoothing estimator
    f_hat <- as.numeric(H %*% Y)
    
    return(f_hat)
  }
}

#------------------------------------------------------#
# Functions to choose initial grid of tuning parameters.
#------------------------------------------------------#
initialize_laplacian_eigenmaps_thetas <- function(sample_X,n){
  # Choose a range of radii based on different desired minimum degree.
  # Currently, a range so that the minimum degree is between log(n) and 2*n^{1/2}.
  iters <- 10
  n_rs <- 10
  degrees <- round(exp(seq(log(1/4*log(n)),log(1/2*n^(1/2)),length.out = n_rs)) + 1)
  
  # ALDEN CHANGE
  # n_rs <- 40
  # degrees <- round(seq(4,80,length.out = n_rs))
  
  degree_dependent_rs <- matrix(ncol = n_rs, nrow = iters)
  for(iter in 1:iters)
  {
    X <- sample_X(n)
    for(jj in 1:n_rs)
    {
      rneighborhood <- nn2(data = X, query = X, k = degrees[jj]) 
      degree_dependent_rs[iter,jj] <- max(rneighborhood$nn.dists[,degrees[jj]])
    }
    # r_connects[iter] <-  # connectivity radius
  }
  # Choose a range of radii intended to achieve the desired min degree.
  rs <- colMeans(degree_dependent_rs)
  
  # Choose a huge range of eigenvectors
  Ks <- 1:round(n/10) # Maximum is 1/2 the recommendation of original Belkin + Niyogi paper.
  
  thetas <- expand.grid(r = rs, K = Ks)
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
  ks <- seq(2, 8, by = 2)
  
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

initialize_kernel_smoothing_thetas <- function(sample_X,n){
  # Choose a range of radii based on different desired minimum degree.
  # Currently, a range so that the minimum degree is between log(n) and n^{1/2}.
  iters <- 10 
  n_rs <- 10
  degrees <- round(exp(seq(log(1/4*log(n)),log(n^(1/2)),length.out = n_rs)) + 1)
  degree_dependent_rs <- matrix(ncol = n_rs, nrow = iters)
  for(iter in 1:iters)
  {
    X <- sample_X(n)
    for(jj in 1:n_rs)
    {
      rneighborhood <- nn2(data = X, query = X, k = degrees[jj]) 
      degree_dependent_rs[iter,jj] <- max(rneighborhood$nn.dists[,degrees[jj]])
    }
    # r_connects[iter] <-  # connectivity radius
  }
  # Choose a range of radii intended to achieve the desired min degree.
  rs <- colMeans(degree_dependent_rs)
  
  thetas <- data.frame(r = rs)
}

initialize_knn_thetas <- function(sample_X,n){
  # m <- choose_f0_frequency(d)
  # a <- choose_f0_amplitude(d,m)
  # k_min <- 2
  # k_max <- min(max(round(5 * (d * (a^2/2^d) * m^2 * pi^2)^{-d/(2+d)} * n^{2/(d + 2)}), 50), n/4)
  ks <- c(2:13, seq(14,30,by = 2))
  
  thetas <- data.frame(k = ks)
}

initialize_spectral_projection_thetas <- function(sample_X,n){
  thetas <- data.frame(K = 1:(n/10))
}
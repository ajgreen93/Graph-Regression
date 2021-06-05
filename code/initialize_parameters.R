#------------------------------------------------------#
# Functions to choose initial grid of tuning parameters.
#------------------------------------------------------#
initialize_laplacian_eigenmaps_plus_kernel_smoothing_thetas <- function(sample_X,n){
  thetas_le <- initialize_laplacian_eigenmaps_thetas(sample_X,n)
  rs <- unique(thetas_le$r)
  Ks <- unique(thetas_le$K)
  thetas <- expand.grid(r = rs, K = Ks, h = rs)
}

initialize_laplacian_eigenmaps_test_thetas <- function(sample_X,n){
  max_K <- round(2*(M^2*n)^(2*d/(4*s + d)))
  initialize_laplacian_eigenmaps_thetas(sample_X,n,max_K)
}

initialize_laplacian_eigenmaps_thetas <- function(sample_X,n,max_K = NULL){
  # Choose a range of radii based on different desired minimum degree.
  # Currently, a range so that the minimum degree is between log(n) and 2*n^{1/2}.
  iters <- 10
  
  # ALDEN CHANGE
  # OPTION MSE
  # n_rs <- 10
  # degrees <- round(exp(seq(log(1/4*log(n)),log(1/2*n^(1/2)),length.out = n_rs)) + 1)
  
  # OPTION TUNING
  n_rs <- 10
  degrees <- round(seq(2,60,length.out = n_rs))
  
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
  
  # Choose same range of eigenvectors as used for spectral projection
  if(is.null(max_K)){
    Ks <- initialize_spectral_projection_thetas(sample_X,n) %>% pull(K)
  } else{
    Ks <- initialize_spectral_projection_thetas(sample_X,n,max_K) %>% pull(K)
  }
  thetas <- expand.grid(r = rs, K = Ks)
}

initialize_laplacian_smoothing_thetas <- function(sample_X,n, max_K = 2*round( (M^2*n)^(d/(2*s + d)))){
  # Choose a range of radii based on different desired minimum degree.
  iters <- 10
  
  # OPTION TUNING
  n_rs <- 10
  degrees <- round(seq(2,60,length.out = n_rs))
  
  degree_dependent_rs <- matrix(ncol = n_rs, nrow = iters)
  for(iter in 1:iters)
  {
    X <- sample_X(n)
    for(jj in 1:n_rs)
    {
      rneighborhood <- nn2(data = X, query = X, k = degrees[jj]) 
      degree_dependent_rs[iter,jj] <- max(rneighborhood$nn.dists[,degrees[jj]])
    }
  }
  # Choose a range of radii intended to achieve the desired min degree.
  rs <- colMeans(degree_dependent_rs)
  
  # Choose rho roughly proportional to 1/lambda_K.
  desired_K <- seq(1,max_K,length.out = 10) # indices of the eigenvalues you will consider.
  rhos <- vector(mode = "list",length = length(rs))
  for(jj in 1:length(rs))
  {
    r <- rs[jj]
    rhos[[jj]] <- matrix(NA,nrow = iters,ncol = length(desired_K))
    for(iter in 1:iters)
    {
      X <- sample_X(n)
      G <- neighborhood_graph(X,r)
      L <- Laplacian(G)
      spectra <- pmax(get_spectra(L,max(desired_K))$values[rev(desired_K)],1e-12)
      rhos[[jj]][iter,] <- 1/spectra
    }
    rhos[[jj]] <- apply(rhos[[jj]],2,FUN = median) %>% unique()
  }
  
  # Combine rhos and rs.
  thetas <- mapply(expand.grid,rs,rhos,SIMPLIFY = F) %>% bind_rows()
  names(thetas) <- c("r","rho")
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
  n_rs <- 30
  degrees <- round(seq(2,120,length.out = n_rs))
  
  degree_dependent_rs <- matrix(ncol = n_rs, nrow = iters)
  for(iter in 1:iters)
  {
    X <- sample_X(n)
    for(jj in 1:n_rs)
    {
      rneighborhood <- nn2(data = X, query = X, k = degrees[jj]) 
      degree_dependent_rs[iter,jj] <- max(rneighborhood$nn.dists[,degrees[jj]])
    }
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

initialize_spectral_projection_test_thetas <- function(sample_X,n){
  max_K <- round(2*(M^2*n)^(2*d/(4*s + d)))
  initialize_spectral_projection_thetas(sample_X,n,max_K)
}

initialize_spectral_projection_thetas <- function(sample_X,n,max_K = 2*round( (M^2*n)^(d/(2*s + d))) ){
  thetas <- data.frame(K = 1:max_K) 
}
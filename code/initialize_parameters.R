#------------------------------------------------------#
# Functions to choose initial grid of tuning parameters, for estimation.
#------------------------------------------------------#
initialize_laplacian_eigenmaps_plus_kernel_smoothing_thetas <- function(sample_X,n){
  thetas_le <- initialize_laplacian_eigenmaps_thetas(sample_X,n)
  rs <- unique(thetas_le$r)
  Ks <- unique(thetas_le$K)
  thetas <- expand.grid(r = rs, K = Ks, h = rs)
}

initialize_laplacian_eigenmaps_thetas <- function(sample_X,n,
                                                  n_rs = 10,
                                                  n_Ks = max_K,
                                                  max_degree = 60,
                                                  max_K = max(initialize_spectral_projection_thetas(sample_X,n)$K),
                                                  dist_to_r = max){
  
  degrees <- seq(2,max_degree + 1,length.out = n_rs) %>% round() %>% unique()
  
  # Choose parameters in a data-dependent way using Monte Carlo
  iters <- 10 
  
  # Choose a range of radii based on different desired degree.
  degree_dependent_rs <- matrix(ncol = n_rs, nrow = iters)
  for(iter in 1:iters)
  {
    X <- sample_X(n)
    for(jj in 1:n_rs)
    {
      rneighborhood <- nn2(data = X, query = X, k = degrees[jj]) 
      degree_dependent_rs[iter,jj] <- dist_to_r(rneighborhood$nn.dists[,degrees[jj]])
    }
  }
  # Choose a range of radii intended to achieve the desired min degree.
  rs <- colMeans(degree_dependent_rs) %>% sort() %>% unique()
  
  # Choose same range of eigenvectors as used for spectral projection
  if(is.null(max_K)){
    Ks <- initialize_spectral_projection_thetas(sample_X,n) %>% pull(K)
  } else{
    Ks <- initialize_spectral_projection_thetas(sample_X,n,max_K) %>% pull(K)
  }
  thetas <- expand.grid(r = rs, K = Ks)
}

initialize_laplacian_smoothing_thetas <- function(sample_X,n,
                                                  n_rs = 10,
                                                  n_rhos = 10,
                                                  max_degree = 60,
                                                  max_K = max(initialize_spectral_projection_thetas(sample_X,n)$K),
                                                  dist_to_r = max,
                                                  rho_multiplier = 1,
                                                  labeled = T){
  #------------------------------------------------------------#
  # Choose tuning parameters for Laplacian smoothing.
  # -- r is determined based on the minimum degree we want in the graph.
  # -- rho is chosen to be inversely proportional to the kth eigenvalue of the Laplacian,
  #    for different values of k, times a multiplier.
  #
  # Inputs:
  # -- dist_to_r: a function that goes from distances (in particular, the distance of each
  #               point to its degreeth neighbor) and chooses one value of r.
  # -- labeled: logical, indicating whether to fit
  # 
  #     min \|Y - f\|_n + f^{\top} L_n f (labeled)
  # 
  #    or 
  #    
  #     min \|Y - f\|_n + f^{\top} L_{2n} f (unlabeled)
  #
  #    where L_{2n} is the Laplacian over G_{2n}, the graph with vertices
  #    {1,...,n,1,...,n}.
  #------------------------------------------------------------#

  degrees <- seq(2,max_degree + 1,length.out = n_rs) %>% round() %>% unique()
  K <- seq(1,max_K,length.out = n_rhos) %>% round() %>% unique()
  
  # Choose parameters in a data-dependent way using Monte Carlo
  r_iters <- 10 
  rho_iters <- 5
  
  # Choose a range of radii based on different desired degree.
  degree_dependent_rs <- matrix(ncol = n_rs, nrow = r_iters)
  for(iter in 1:r_iters)
  {
    X <- sample_X(n)
    for(jj in 1:n_rs)
    {
      rneighborhood <- nn2(data = X, query = X, k = degrees[jj]) 
      degree_dependent_rs[iter,jj] <- dist_to_r(rneighborhood$nn.dists[,degrees[jj]])
    }
  }
  rs <- colMeans(degree_dependent_rs) %>% sort() %>% unique()
  
  # Choose rho equal to 1/lambda_k for different k.
  rhos <- vector(mode = "list",length = length(rs))
  for(jj in 1:length(rs))
  {
    r <- rs[jj]
    rhos[[jj]] <- matrix(NA,nrow = rho_iters,ncol = length(K))
    for(iter in 1:rho_iters)
    {
      X <- sample_X(n)
      G <- neighborhood_graph(X,r)
      L <- Laplacian(G)
      spectra <- pmax(get_spectra(L,max_K)$values[rev(K)],1e-12) # Machine zero.
      rhos[[jj]][iter,] <- 1/spectra * rho_multiplier
    }
    rhos[[jj]] <- apply(rhos[[jj]],2,FUN = median) %>% unique()
    if(!labeled){
      # Re-normalize, because f^{\top} L_{G_{2n}} f = 4 f^{\top} L f.
      rhos[[jj]] <- ifelse(rhos[[jj]] >= 1e12, rhos[[jj]], rhos[[jj]]/4)
    }
  }
  
  # Combine rhos, rs, and nfolds.
  thetas <- mapply(expand.grid,rs,rhos,SIMPLIFY = F) %>% 
    bind_rows() %>%
    cbind(labeled)
  names(thetas) <- c("r","rho","labeled")
  return(thetas)
}

initialize_kernel_smoothing_thetas <- function(sample_X,n,
                                               n_rs = 30,
                                               max_degree = 120){
  #------------------------------------------------------------#
  # Choose tuning parameters for kernel smoothing.
  # -- r is determined based on the minimum number of neighbors wanted for each X.
  #------------------------------------------------------------#
  degrees <- seq(2,max_degree + 1,length.out = n_rs) %>% round() %>% unique()
  
  # Choose parameters in a data-dependent way using Monte Carlo
  iters <- 10 
  
  # Choose a range of radii based on different desired **minimum** degree.
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

initialize_spectral_projection_thetas <- function(sample_X,n,max_K = n/50){
  thetas <- data.frame(K = 1:max_K) 
}

initialize_knn_thetas <- function(sample_X,n){
  # m <- choose_f0_frequency(d)
  # a <- choose_f0_amplitude(d,m)
  # k_min <- 2
  # k_max <- min(max(round(5 * (d * (a^2/2^d) * m^2 * pi^2)^{-d/(2+d)} * n^{2/(d + 2)}), 50), n/4)
  ks <- c(2:13, seq(14,30,by = 2))
  
  thetas <- data.frame(k = ks)
}

#------------------------------------------------------------------#
# Functions to choose initial grid of tuning parameters, for testing.
#------------------------------------------------------------------#
initialize_spectral_projection_test_thetas <- function(sample_X,n){
  max_K <- round(2*(M^2*n)^(2*d/(4*s + d)))
  initialize_spectral_projection_thetas(sample_X,n,max_K)
}

initialize_laplacian_eigenmaps_test_thetas <- function(sample_X,n){
  max_K <- round(2*(M^2*n)^(2*d/(4*s + d)))
  initialize_laplacian_eigenmaps_thetas(sample_X,n,max_K)
}

#------------------------------------------------------------------#
# Old.
#------------------------------------------------------------------#

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
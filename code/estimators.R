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
    
    # "Prediction" function
    predict <- function(X){f_hat}
    return(predict)
  }
}
attr(make_laplacian_eigenmaps,"precompute_data") <- function(Y,X,X_test,X_validate,theta_df){
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
  precomputed_data <- list(train = precomputed_data, test = NULL)
  return(precomputed_data)
}

make_laplacian_eigenmaps_plus_kernel_smoothing <- function(theta)
{
  r <- theta[["r"]]
  K <- theta[["K"]]
  h <- theta[["h"]]
  theta_kernel_smoothing <- data.frame(r = h)
  laplacian_eigenmaps_plus_kernel_smoothing <- function(Y,X){
    # If we've already computed the eigenvectors don't do it again
    if(exists("precomputed_data")){
      stopifnot("eigenvectors" %in% names(precomputed_data))
      L_eigenvectors <- precomputed_data$eigenvectors[[which(names(precomputed_data$eigenvectors) == r)]]
    } else{
      # Build G_{n,r} over X.
      G <- neighborhood_graph(X,r)
      
      # Get L.
      L <- Laplacian(G)
      
      # Get spectra.
      spectra <- get_spectra(L,K)
      L_eigenvectors <- spectra$vectors[,ncol(spectra$vectors):1,drop = FALSE]
    }
    V_K <- L_eigenvectors[,1:K] # Eigenvectors we need
    a_k <- t(V_K) %*% Y         # Empirical Fourier coefficients
    f_hat <- V_K %*% a_k        # Spectral projection
    
    # Prediction function
    predict <- function(X_new){
      # If we've already computed the smoohter matrix, don't recompute
      if(exists("precomputed_data"))
      {
        stopifnot("kernel_matrices" %in% names(precomputed_data))
        H <- precomputed_data$kernel_matrices[[which(names(precomputed_data$kernel_matrices) == h)]]
        as.numeric(H %*% f_hat)
      } else{
        kernel_estimator <- make_kernel_smoothing(theta_kernel_smoothing)
        kernel_estimate <- kernel_estimator(f_hat,X)
        kernel_estimate(X_new)
      }
    }
    return(predict)
  }
}
attr(make_laplacian_eigenmaps_plus_kernel_smoothing,"precompute_data") <- function(Y,X,X_test,X_validate,theta_df){
  # Parameters
  rs <- unique(theta_df$r)
  hs <- unique(theta_df$h)
  
  # Objects for storing pre-computed data
  precomputed_data_train_eigenvectors <- vector(mode = "list",length = length(rs))
  precomputed_data_train_kernel_matrices <- vector(mode = "list",length = length(hs))
  precomputed_data_validate_kernel_matrices <- vector(mode = "list",length = length(hs))
  precomputed_data_test_kernel_matrices <- vector(mode = "list",length = length(hs))
  names(precomputed_data_train_eigenvectors) <- rs
  names(precomputed_data_train_kernel_matrices) <- hs
  names(precomputed_data_validate_kernel_matrices) <- hs
  names(precomputed_data_test_kernel_matrices) <- hs
  
  # Pre-compute eigenvectors
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
    
    precomputed_data_train_eigenvectors[[kk]] <- L_eigenvectors
  }
  
  # Pre-compute kernel smoothing "hat" matrices.
  for(kk in 1:length(hs)){
    h <- hs[kk]
    G_train <- neighborhood_graph(X,h,loop = T)
    G_test <- neighborhood_graph(X,h,X_test)
    G_validate <-  neighborhood_graph(X,h,X_validate)
    precomputed_data_train_kernel_matrices[[kk]] <- G_train/pmax(rowSums(G_train),1) # Avoid divide by zero.
    precomputed_data_test_kernel_matrices[[kk]] <-  G_test/pmax(rowSums(G_test),1)
    precomputed_data_validate_kernel_matrices[[kk]] <-  G_validate/pmax(rowSums(G_test),1)
  }
  
  # Save it all.
  precomputed_data <- list(train = list(eigenvectors = precomputed_data_train_eigenvectors,
                                        kernel_matrices = precomputed_data_train_kernel_matrices),
                           validate = list(kernel_matrices = precomputed_data_validate_kernel_matrices),
                           test = list(kernel_matrices = precomputed_data_test_kernel_matrices))
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
    
    # "Prediction" function
    predict <- function(X){f_hat}
    return(predict)
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
    
    # "Prediction" function
    predict <- function(X){f_hat}
    return(predict)
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
      psi_K <- precomputed_data[,1:K,drop = FALSE]
    } else{
      trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
      psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
    }
    a_k <- 1/n * (t(psi_K) %*% Y)   # Empirical Fourier coefficients
    
    # Prediction function.
    predict <- function(X){
      
      if(exists("precomputed_data")){
        psi_K <- precomputed_data[,1:K,drop = FALSE]
      } else{
        trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
        psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
      }
      
      preds <- psi_K %*% a_k # Spectral projection
      return(preds)
    }
    
    return(predict)
  }
} 
attr(make_spectral_projection,"precompute_data") <- function(Y,X,X_test,X_validate,theta_df){
  precomputed_data <- vector(mode = "list",length = 3)
  names(precomputed_data) <- c("train","validate","test")
  
  K <- max(theta_df$K)
  trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
  precomputed_data[["train"]] <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
  precomputed_data[["validate"]] <- if(K == 1) apply(X_validate,1,trig_basis) else apply(X_validate,1,trig_basis) %>% t()
  precomputed_data[["test"]] <- if(K == 1) apply(X_test,1,trig_basis) else apply(X_test,1,trig_basis) %>% t()
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
      psi_K <- precomputed_data[,1:K,drop = FALSE]
    } else{
      trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
      psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
    }
    
    a_k <- lm.fit(x = psi_K,y = Y)$coefficients %>% as.matrix()
    
    # Prediction function.
    predict <- function(X){
      
      if(exists("precomputed_data")){
        psi_K <- precomputed_data[,1:K,drop = FALSE]
      } else{
        trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
        psi_K <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
      }
      
      preds <- psi_K %*% a_k # Spectral projection
      return(preds)
    }
  }
}
attr(make_least_squares,"precompute_data") <- function(Y,X,X_test,X_validate,theta_df){
  precomputed_data <- vector(mode = "list",length = 3)
  names(precomputed_data) <- c("train","validate","test")
  
  K <- max(theta_df$K)
  trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
  precomputed_data[["train"]] <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
  precomputed_data[["validate"]] <- if(K == 1) apply(X_validate,1,trig_basis) else apply(X_validate,1,trig_basis) %>% t()
  precomputed_data[["test"]] <- if(K == 1) apply(X_test,1,trig_basis) else apply(X_test,1,trig_basis) %>% t()
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
    predict <- function(X_new){
      # Build kernel matrix
      K <- neighborhood_graph(X,r,X_new,loop = T)
      
      # normalize by row sums
      H <- K / rowSums(K)
      
      # kernel smoothing estimator
      f_hat <- as.numeric(H %*% Y)
    }
    return(predict)
  }
}

#------------------------------------------------------#
# Functions to choose initial grid of tuning parameters.
#------------------------------------------------------#
initialize_laplacian_eigenmaps_plus_kernel_smoothing_thetas <- function(sample_X,n){
  thetas_le <- initialize_laplacian_eigenmaps_thetas(sample_X,n)
  rs <- unique(thetas_le$r)
  Ks <- unique(thetas_le$K)
  thetas <- expand.grid(r = rs, K = Ks, h = rs)
}

initialize_laplacian_eigenmaps_thetas <- function(sample_X,n){
  # Choose a range of radii based on different desired minimum degree.
  # Currently, a range so that the minimum degree is between log(n) and 2*n^{1/2}.
  iters <- 10
  
  # ALDEN CHANGE
  # OPTION MSE
  # n_rs <- 50
  # degrees <- round(exp(seq(log(1/4*log(n)),log(1/2*n^(1/2)),length.out = n_rs)) + 1)
  
  # OPTION TUNING
  n_rs <- 40
  degrees <- round(seq(4,80,length.out = n_rs))
  
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
  Ks <- initialize_spectral_projection_thetas(sample_X,n) %>% pull(K)
  
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
  # ALDEN CHANGE
  max_K <- 2*round( (M^2*n)^(d/(2*s + d))) 
  thetas <- data.frame(K = 1:max_K) 
}
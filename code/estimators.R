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
attr(make_laplacian_eigenmaps,"precompute_data") <- precompute_laplacian_eigenmaps_data

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
attr(make_spectral_projection,"precompute_data") <- precompute_least_squares_data

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
attr(make_least_squares,"precompute_data") <- precompute_least_squares_data

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


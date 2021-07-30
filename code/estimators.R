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
attr(make_laplacian_eigenmaps_plus_kernel_smoothing,"precompute_data") <- 
  precompute_laplacian_eigenmaps_plus_kernel_smoothing_data

make_laplacian_smoothing <- function(theta)
{
  r <- theta[["r"]]
  rho <- theta[["rho"]]
  labeled <- theta[["labeled"]]
  laplacian_smoothing <- function(Y,X){
    # Build G_{n,r} over X.
    G <- neighborhood_graph(X,r)
    
    # Get L.
    L <- Laplacian(G)
    
    # Solve for the estimator, 
    # unless rho is ``too large'', in which case we set equal to mean.
    if(rho >= 1e+12){
      f_hat <- rep(mean(Y),length(Y))
    } else
    {
      if(labeled)
      {
        f_hat <- as.numeric(solve(rho * L + Diagonal(n),Y))
      } else if(!labeled)
      {
        # Create a graph H where every vertex in G is duplicated, included 
        # once with its corresponding label and once without.
        # Solve the Laplacian smoothing problem over H.
        # Return the estimate which at the "unlabeled" points.
        H <- cbind(rbind(G,G + Diagonal(n)),rbind(G + Diagonal(n),G))
        L_H <- Laplacian(H)
        I_n <- Diagonal(2*n,x = c(rep(1,n),rep(0,n)))
        Y_n <- c(Y,rep(0,n))
        f_hat <- as.numeric(solve(rho * L_H + I_n,Y_n))[(n + 1):(2*n)]
      } 
    }
    
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

#------------------------------------------------------------#
# Not currently maintained code is below.
#------------------------------------------------------------#


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


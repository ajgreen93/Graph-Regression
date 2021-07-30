#----------------------------------------------------#
# Functions to precompute data which will be used repeatedly.
# E.g. graph eigenvectors.
#----------------------------------------------------#

precompute_least_squares_data <- function(Y,X,X_test,X_validate,theta_df){
  precomputed_data <- vector(mode = "list",length = 3)
  names(precomputed_data) <- c("train","validate","test")
  
  K <- max(theta_df$K)
  trig_basis <- get_trigonometric_basis(d,K) # Fourier basis
  precomputed_data[["train"]] <- if(K == 1) apply(X,1,trig_basis) else apply(X,1,trig_basis) %>% t()
  precomputed_data[["validate"]] <- if(K == 1) apply(X_validate,1,trig_basis) else apply(X_validate,1,trig_basis) %>% t()
  precomputed_data[["test"]] <- if(K == 1) apply(X_test,1,trig_basis) else apply(X_test,1,trig_basis) %>% t()
  return(precomputed_data)
}

precompute_laplacian_eigenmaps_data <- function(Y,X,X_test,X_validate,theta_df){
  if(!is.list(Y)) Y  <- list(Y)
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
    
    # Store eigenvectors
    precomputed_data[[kk]] <- L_eigenvectors 
  }
  precomputed_data <- list(train = precomputed_data, test = NULL)
  return(precomputed_data)
}

precompute_laplacian_eigenmaps_plus_kernel_smoothing_data <- 
function(Y,X,X_test,X_validate,theta_df){
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
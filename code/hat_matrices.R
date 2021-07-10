make_laplacian_smoothing_hat <- function(theta)
{
  r <- theta[["r"]]
  rho <- theta[["rho"]]
  kernel <- theta[["kernel"]]
  laplacian_smoothing <- function(X){
    
    # Build G_{n,r} over X.
    G <- neighborhood_graph(X,r)
    
    # Get L.
    L <- Laplacian(G)
    
    # Solve (rho*L + I)^{-1},
    # unless rho is ``too large'', in which case we set H to global mean.
    if(rho >= 1e+12){
      H <- matrix(rep(1/n,n^2),ncol = n,nrow = n)
    }else{
      H <- solve(rho * L + Diagonal(n))
    }
    
    # "Prediction" function
    return(H)
  }
}

make_laplacian_eigenmaps_hat <- function(theta)
{
  r <- theta[["r"]]
  K <- theta[["K"]]
  kernel <- theta[["kernel"]]
  laplacian_eigenmaps <- function(X){
    # Build G_{n,r} over X.
    G <- neighborhood_graph(X,r)
    
    # Get L.
    L <- Laplacian(G)
    
    # Get spectra.
    spectra <- get_spectra(L,K)
    L_eigenvectors <- spectra$vectors[,ncol(spectra$vectors):1,drop = F]
    
    V_K <- L_eigenvectors[,1:K,drop = F] # Eigenvectors we need
    H <- V_K %*% t(V_K)
    
    return(H)
  }
}
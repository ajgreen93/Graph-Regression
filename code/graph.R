neighborhood_graph <- function(X, r)
{
  # unpack necessary parameters
  N <- nrow(X)
  
  k_max <- round(min(4*N^(.5),N - 1))
  A <- knn_to_neighborhood_graph(X,r,k_max)
  return(A)
}

knn_to_neighborhood_graph <- function(X, r, k_max)
{
  #--------------------#
  # Input: X (n x d matrix)
  #        r connectivity parameter
  #        k_max maximum degree in graph
  # Output: A (m x n sparse adjacency matrix)
  # Function: compute adjacency matrix A
  # with A = (A_ij)_{i,j in 1:n}, A_ij = 1 iff ||X[i,] - X[j,]||_2 < r
  #--------------------#
  
  # unpack necessary parameters
  n <- nrow(X)
  d <- ncol(X)
  A <- Matrix(0, nrow = n, ncol = n)
  
  # compute distance matrix
  rneighbors <- nn2(data = X, query = X, searchtype = 'radius', k = k_max, radius = r)$nn.idx[,-1]
  rneighbors[rneighbors == 0] <- NA
  r_neighbors_list <- as.matrix(melt(rneighbors)[,-2])
  r_neighbors_list <- r_neighbors_list[rowSums(is.na(r_neighbors_list)) == 0,]
  A[r_neighbors_list] <- 1
  
  # if k_max was not high enough, at least one vertex in A will have exactly k_max - 1 neighbors.
  # Increase k_max and recalculate.
  while((max(rowSums(A)) == k_max - 1) & (k_max < n - 1))
  {
    A <- Matrix(0, nrow = n, ncol = n)
    k_max <- min(2*k_max,n - 1)
    rneighbors <- nn2(data = X, query = X, searchtype = 'radius', k = k_max, radius = r)$nn.idx[,-1]
    rneighbors[rneighbors == 0] <- NA
    r_neighbors_list <- as.matrix(melt(rneighbors)[,-2])
    r_neighbors_list <- r_neighbors_list[rowSums(is.na(r_neighbors_list)) == 0,]
    A[r_neighbors_list] <- 1
  }
  
  return(A)
}

knn_graph <- function(X, k)
{
  # unpack necessary parameters
  n <- nrow(X)
  d <- ncol(X)
  A <- Matrix(0, nrow = n, ncol = n)
  
  # find neighbors
  neighbors <- nn2(data = X, query = X, k = k + 1)$nn.idx[,-1]
  neighbors_list <- as.matrix(melt(neighbors)[,-2])
  A[neighbors_list] <- 1
  
  # symmetrize
  A <- (A + t(A)) / 2
  
  return(A)
}

Laplacian <- function(A)
{
  # unpack necessary parameters
  n <- nrow(A)
  
  # Laplacian matrix
  D <- Matrix(0, ncol = n, nrow = n, sparse = T)
  diag(D) <- rowSums(A)
  if(class(D - A) == "dsyMatrix")
  {
    L <- as(D - A,"matrix")
  } else{
    L <- as(D - A,"dgCMatrix")
  }
  return(L)
}

# A wrapper around eigs_sym, which takes care of some failure cases when the
# matrix is poorly conditioned. 
get_spectra <- function(A,K,sigma = 0){
  # Compute as many eigenvectors as we will need.
  spectra <- tryCatch(eigs_sym(A,K,sigma),
                      error = function(e){
                        test <- eigs_sym(A,K, which = "SM")
                        return(test)
                      })
  if(max(abs(spectra$values)) > 2e10)
  {
    spectra <- eigs_sym(A,K,which = "SM")
  }
  spectra
}

make_sample_uniform <- function(d)
{
  sample_X <- function(n){
    X <- matrix(runif(n*d,-1,1),ncol = d)
    return(X)
  }
}


make_cosine_f0 <- function(d,n)
{
  # Smallest choice of m such that first-order regression error dominates.
  m <- round(1/2 * n^{1/d * (2 - (4 + d)/(2 + d))})
  f0 <- function(x){
    return(1/m * prod(cos(pi * m * x)))
  }
  return(f0)
}
make_sample_uniform <- function(d)
{
  sample_X <- function(n){
    X <- matrix(runif(n*d,-1,1),ncol = d)
    return(X)
  }
}


make_cosine_f0 <- function(d,n)
{
  # Choose m such that the first-order regression error dominates.
  # m <- round(2 * n^{1/5} ) # TODO: you only know this works when d = 1.
  m <- choose_f0_frequency(d)
  a <- choose_f0_amplitude(d,m)
  # Smallest value of a for which LS estimator dominates 0 estimator at n = n_smallest.
  
  f0 <- function(x){
    return(a * prod(cos(pi * m * x)))
  }
  return(f0)
}

choose_f0_frequency <- function(d)
{
  if(d == 1)
  {
    m <- 8
  } else if(d == 2)
  {
    m <- 2
  } else{
    m <- 1
  }
  
  return(m)
}

choose_f0_amplitude <- function(d,m)
{
  # a <- 4^d * (d * m^2 * pi^2)^{d/4} * N_SMALLEST^{-.5}
  if (d <= 3) 
  {
    a <- 3
  } else {
    a <- 4
  }
    
  return(a)
}
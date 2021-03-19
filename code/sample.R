make_sample_gaussian <- function(d)
{
  sample_X <- function(n){
    X <- matrix(NA,nrow = n,ncol = d)
    while(any(is.na(X))){
      X[which(is.na(X))] <- rnorm(sum(is.na(X)))
      X[abs(X) > 1] <- NA
    }
    return(X)
  }
}

make_sample_uniform <- function(d)
{
  sample_X <- function(n){
    X <- matrix(runif(n*d,-1,1),ncol = d)
    return(X)
  }
}

make_wiggly_gaussian <- function(d,n)
{
  f0 <- function(x){
    return(prod((1/dnorm(x)) * cos(2*pi*x/dnorm(x)^{1.5})))
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

make_linear_f0 <- function(d,n)
{
  #f_0 = b*(X_1 + ... + X_d)
  
  # Choose slope
  b <- choose_linear_f0_slope(d)
  
  f0 <- function(x){
    return(b * sum(x))
  }
  return(f0)
}

choose_linear_f0_slope <- function(d)
{
  m <- choose_f0_frequency(d)
  a <- choose_f0_amplitude(d,m)
  b <- a*m*pi
  if(d == 3)
  {
    b <- b / 3
  } else if(d == 4)
  {
    b <- b / 4
  }
  return(b)
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
  } else if (d == 4){
    a <- 4
  } else if (d == 5){
    a <- 6
  }
    
  return(a)
}
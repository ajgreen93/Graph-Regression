make_sample_gaussian <- function(d, sigma = 1)
{
  sample_X <- function(n){
    X <- matrix(NA,nrow = n,ncol = d)
    while(any(is.na(X))){
      X[which(is.na(X))] <- rnorm(sum(is.na(X)),0,sigma)
      X[abs(X) > 1] <- NA
    }
    return(X)
  }
}

make_sample_gaussian_mixture <- function(d,mu = NULL,pi = NULL,sigma = NULL)
{
  if(is.null(mu)){mu <- 0}
  if(is.null(pi)){pi <- rep(1,length(mu))}
  if(is.null(sigma)){sigma <- rep(1,length(mu))}
  sample_X <- function(n){
    X <- matrix(NA,nrow = n,ncol = d)
    while(any(is.na(X))){
      class <- apply(rmultinom(sum(is.na(X)), 1, prob = pi),2,FUN = function(i){which(i != 0)})
      X[which(is.na(X))] <- sapply(class,FUN = function(l){rnorm(1,mu[l],sigma[l])})
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

#-----------------------------------------------#
# Regression functions
#-----------------------------------------------#

make_doppler <- function(d,n){
  stopifnot(d == 1)
  f0 <- function(x){
    if(x > 0){cos(4*pi/((abs(x)^{1/3})))}
    else{cos(pi/(abs(x)^{1/3}))}
  }
}

make_spectral_sobolev <- function(d,n){
  #------------------------------#
  # Fixes a vector beta that belongs to a Sobolev ellipsoid,
  # then builds the regression function
  # 
  # f0(x) = sum_{k} beta_k psi_k
  #
  # where psi_k(x) = cos(k*pi*x)
  stopifnot(d == 1) # Only do 1d for the moment
  c <- sqrt(sum((1:3)^2))
  beta <- c(rep(1/c,3),1/(sqrt(17)*(4:20)))
  f0 <- function(x){
    sum(sapply(1:length(beta), function(idx){beta[idx] * sqrt(2) * cos(idx/2*pi*(x + 1))}))
  }
}

make_wiggly_gaussian <- function(d,n)
{
  f0 <- function(x){
    return(1/2 * prod((1/dnorm(x)) * cos(2*pi*x/dnorm(x)^{1.5})))
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
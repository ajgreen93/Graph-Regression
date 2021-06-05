#------------------------------------#
# Design distributions.
#------------------------------------#
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

# p(x) = 1/2 * (1{0 < x < 1/2 - sep} + 1{1/2 + sep < x < 1}).
make_sample_uniform_mixture <- function(d,sep)
{
  stopifnot(d == 1) # For now, just 1d.
  stopifnot(sep < 1/2)
  L <- c(0,1/2 + sep); U = c(1/2 - sep,1) # Lower and upper extremes.
  sample_X <- function(n){
      class <- apply( rmultinom(n, 1, prob = c(1/2,1/2)),2,FUN = function(i){which(i != 0)})
      X <- sapply(class,FUN = function(l){runif(1,L[l],U[l])}) %>% as.matrix()
      X
  }
}

make_sample_manifold <- function(d,D)
{
  stopifnot(d <= D)
  sample_X <- function(n){
    Z <- matrix(runif(n*d,-1,1),ncol = d)
    X <- t(apply(Z,1,FUN = function(z){chart(z,D)}))
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

make_eigenfunction <- function(d,K,M,s){
  #----------------------------------#
  # Builds the regression function
  #
  # f0(x) = M/lambda_K^(s/2) * psi_K
  #
  # where psi_K(x) is the kth element in the trigonometric basis.
  #----------------------------------#
  lambda_K <- get_eigenvalue(d,K)
  eigenfunction_K <- get_eigenfunction(d,K)
  f0 <- function(x){
    M/lambda_K^(s/2) * eigenfunction_K(x)
  }
  attr(f0,"type") <- "eigenfunction"
  return(f0)
}

make_testing_eigenfunction_f0s <- function(d,n){
  # Configuration type parameter
  max_K <- round(2*(M^2*n)^(2*d/(4*s + d)))
  
  # Get eigenfunctions, and associated l2 norms
  g0s <- lapply(1:max_K,FUN = function(k){make_eigenfunction(d,k,M,s)})
  l2_norms <- sapply(g0s,FUN = get_l2_norm)
  
  # Remove anything with infinite l2 norm (constant functions), or
  # duplicates in l2 norm.
  g0s <- g0s[l2_norms < Inf]
  l2_norms <- l2_norms[l2_norms < Inf]
  g0s <- g0s[match(unique(l2_norms),l2_norms)]
  l2_norms <- l2_norms[match(unique(l2_norms),l2_norms)]
  
  # Construct regression functions.
  make_f0 <- function(g,l){
    f0 <- function(x){2 * g(x)}
    attr(f0,"l2_norm") <- 4 * l
    return(f0)
  }
  f0s <- mapply(FUN = make_f0, g0s,l2_norms)
  f0_null <- function(x){0}; 
  attr(f0_null,"l2_norm") <- 0
  f0s <- c(f0_null,f0s)
  return(f0s)
}

make_testing_eigenfunction_mixtures_f0s <- function(d,n){
  # Configuration type parameters
  max_K <- round(2*(M^2*n)^(2*d/(4*s + d)))
  n_f0s <- 50
  
  # Get eigenfunctions, and associated l2 norms
  g0s <- lapply(1:max_K,FUN = function(k){make_eigenfunction(d,k,M,s)})
  l2_norms <- sapply(g0s,FUN = get_l2_norm)
  
  # Remove anything with infinite l2 norm (constant functions), or
  # duplicates in l2 norm.
  g0s <- g0s[l2_norms < Inf]
  l2_norms <- l2_norms[l2_norms < Inf]
  g0s <- g0s[match(unique(l2_norms),l2_norms)]
  l2_norms <- l2_norms[match(unique(l2_norms),l2_norms)]
  
  # Construct the regression functions---which are mixtures of eigenfunctions---
  # sorted by decreasing l2 norm.
  f0s <- list()
  
  mixing_indxs <-  seq(1,length(l2_norms),length.out = n_f0s)
  
  f0s <- lapply(mixing_indxs,FUN = function(mixing_indx){
    if(mixing_indx == 1){
      f0 <- function(x){2 * g0s[[1]](x)}
      attr(f0,"l2_norm") <- 4 * l2_norms[1]
    } else{
      f0 <- function(x){
        2 * (1 - (mixing_indx - floor(mixing_indx))) * g0s[[floor(mixing_indx)]](x) + 
        2 * (mixing_indx - floor(mixing_indx)) * g0s[[ceiling(mixing_indx)]](x)
      }
      attr(f0,"l2_norm") <- 
        (2 * (1 - (mixing_indx - floor(mixing_indx))))^2 * l2_norms[floor(mixing_indx)] +
        (2 * (mixing_indx - floor(mixing_indx)))^2 * l2_norms[ceiling(mixing_indx)]
    }
    return(f0)
  })
  l2_norm_order <- order(sapply(f0s,FUN = function(f){attributes(f)$l2_norm}),decreasing = T)
  f0s <- f0s[l2_norm_order]
  f0_null <- function(x){0}; 
  attr(f0_null,"l2_norm") <- 0
  f0s <- c(f0_null,f0s)
}

make_spectral_sobolev <- function(d,n,s = 1,M = 2^s,beta = NULL){
  #------------------------------#
  # Fixes a vector beta that belongs to a Sobolev ellipsoid,
  # then builds the regression function
  # 
  # f0(x) = sum_{k} beta_k psi_k
  #
  # where psi_k(x) is the kth element in the trigonometric basis.
  #------------------------------#
  if(is.null(beta)){
    n_coefs <- n/20
    beta <- M/sqrt(n_coefs) * c(1/sqrt( (1:n_coefs)^(2*s/d) ) ) # Sobolev norm = 2^s.
  }
  trig_basis <- get_trigonometric_basis(d,length(beta))
  f0 <- function(x){
    t(beta) %*% trig_basis(x)
  }
}

make_spectral_sobolev_manifold <- function(d,n,beta = NULL){
  g0 <- make_spectral_sobolev(d,n,beta)
  f0 <- function(x){
    g0(inverse_chart(x,d))
  }
}

make_wiggly_gaussian <- function(d,n)
{
  f0 <- function(x){
    return(1/2 * prod((1/dnorm(x)) * cos(2*pi*x/dnorm(x)^{1.5})))
  }
}

make_stepfunction <- function(height){
  f0 <- function(x){
    return(height*((0 < x & x < 1/2) - (1/2 < x & x < 1)))
  }
}

#---------------------------------------------------#
# Old.
#---------------------------------------------------#

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
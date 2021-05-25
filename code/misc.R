#------------------------------------------#
# Helper functions
#------------------------------------------#

get_trigonometric_basis <- function(d,K){
  # Get the indices for all tensor product polynomials.
  max_indx <- ceiling( sqrt(d*K^(2/d)) )
  indx <- expand.grid(rep(list(0:(max_indx - 1)), d)) %>% as.matrix()
  eigenvalues <- apply(indx,1,FUN = function(k){sum(k^2)})
  indx <- indx[order(eigenvalues),,drop = FALSE][1:K,,drop = FALSE]
  
  # NOTE: following code is optimized for speed, rather than ease of understanding.
  trig_basis <- function(x){
    basis_by_dim <- sqrt(2) * t(cos(t(indx)*pi*(x + 1)/2)) * (1/sqrt(2))^(indx == 0)
    basis <- rep(1,K)
    for(ii in 1:d){
      basis <- basis * basis_by_dim[,ii]
    }
    return(basis)
  }
  return(trig_basis)
}

get_eigenfunction <- function(d,K){
  max_indx <- ceiling( sqrt(d*K^(2/d)) )
  indx <- expand.grid(rep(list(0:(max_indx - 1)), d))
  eigenvalues <- apply(indx,1,FUN = function(k){sum(k^2)})
  k <- indx[order(eigenvalues)[K],] %>% as.numeric()
  # eigenfunction <- function(x){cos(k*pi*(x + 1)/2)}
  eigenfunction <- function(x){1/sqrt(2)^(sum(k == 0)) * prod(sqrt(2) * cos(k*pi*(x + 1)/2))}
  return(eigenfunction)
}

get_eigenvalue <- function(d,K){
  max_indx <- ceiling(K^(1/d))
  indx <- expand.grid(rep(list(0:(max_indx - 1)), d))
  eigenvalues <- apply(indx,1,FUN = function(k){sum(k^2)})
  k <- indx[order(eigenvalues)[K],,drop = FALSE]
  eigenvalue <- sum(k^2)
  return(eigenvalue)
}


find_best_K <- function(mse,thetas){
  best_K <- matrix(nrow = length(ns),ncol = length(thetas[[1]]))
  for(ii in 1:length(ns))
  {
    for(jj in 1:length(mse[[ii]]))
    {
      best_K[ii,jj] <- thetas[[ii]][[jj]][which.min(rowSums(mse[[ii]][[jj]])),"K"]
    }
  }
  best_K
}

find_best_radius <- function(mse,thetas){
  best_r <- matrix(nrow = length(ns),ncol = length(thetas[[1]]))
  for(ii in 1:length(ns))
  {
    for(jj in 1:length(mse[[ii]]))
    {
      if(!("r" %in% names(thetas[[ii]][[jj]]))){next}
      best_r[ii,jj] <- thetas[[ii]][[jj]][which.min(rowSums(mse[[ii]][[jj]])),"r"]
    }
  }
  best_r
}

find_best_mse <- function(mse,methods){
  best_mse <- vector(mode = "list", length = length(methods))
  names(best_mse) <- names(methods)
  
  for(jj in 1:length(methods))
  {
    best_mse_jj <- numeric()
    sd_best_mse_jj <- numeric()
    for(ii in 1:length(ns))
    {
      mse_ii_jj <- mse[[ii]][[jj]]
      best_mse_jj[ii] <- min(rowMeans(mse_ii_jj))
      sd_best_mse_jj[ii] <- apply(mse_ii_jj,1,sd)[which.min(rowMeans(mse_ii_jj))]/sqrt(ncol(mse_ii_jj))
    }
    # Rescale minimax mse to match intercept with best_mse
    best_mse[[jj]] <- data.frame(x = ns, y = best_mse_jj,sd = sd_best_mse_jj) 
  }
  return(best_mse)
}

get_bias_spectral_sobolev <- function(Ks){
  beta <- environment(make_f0(d,n))$beta
  sapply(Ks,function(K){sum(beta[-(1:K)]^2)})
}

theory_K_spectral_sobolev <- function(ns){
  max_K <- length(environment(make_f0(d,n))$beta)
  temp <- matrix(nrow = max_K,ncol = length(ns))
  for(ii in 1:length(ns))
  {
    temp[,ii] <- (1:max_K)/ns[ii] + get_bias_spectral_sobolev(1:max_K)
  }
  apply(temp,2,which.min)
}

mse_slope <- function(mse,ns){
  log_mse <- log(mse)
  log_ns <- log(ns)
  lm(log_mse ~ log_ns)$coefficients[2]
}
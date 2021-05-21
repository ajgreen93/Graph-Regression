#------------------------------------------#
# Helper functions
#------------------------------------------#

get_trigonometric_basis <- function(d,K){
  max_indx <- ceiling(K^(1/d))
  indx <- expand.grid(rep(list(0:(max_indx - 1)), d))
  eigenvalues <- apply(indx,1,FUN = function(k){sum(k^2)})
  indx <- indx[order(eigenvalues),,drop = FALSE][1:K,,drop = FALSE]
  trig_basis <- function(x){apply(indx,1,
                                  function(k){1/sqrt(2)^(sum(k == 0)) * prod(sqrt(2) * cos(k*pi*(x + 1)/2))}
  )}
  return(trig_basis)
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
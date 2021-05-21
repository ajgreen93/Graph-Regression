chart <- function(z,D){
  # Embed m-dimensional data into a higher dimensional space
  m <- length(z)
  rep(z^(1:d),ceiling(D/d))[1:D]
}

inverse_chart <- function(x,d){
  # Inverse of chart()
  sign(x[1:d]) * abs(x[1:d])^(1/(1:d))
}
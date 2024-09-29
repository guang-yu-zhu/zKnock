# Covariance matrices scaled to be approximately 1/n on the diagonal

cov_diag <- function(n, p, rho=NULL) {
  # Diagonal covariance
  s <- diag(p)
  return(s)
}

cov_equi <- function(n, p, rho = 0.5) {
  # Equicorrelated covariance
  s <- (diag(1 - rho, p, p) + rho)
  return(s)
}

cov_ar1 <- function(n, p, rho = 0.5) {
  # AR(1) covariance
  s <- toeplitz(rho^(0:(p - 1)))
  return(s)
}


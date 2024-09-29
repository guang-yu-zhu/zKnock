#' Simulate Gaussian and binary covariate predictors
#'
#' This function simulates Gaussian predictors with mean zero
#' and a covariance structure determined by the "cov_type" argument.
#' Then, `p_b` randomly selected columns are dichotomized.
#'
#' @param n Number of observations (rows of X).
#' @param p Total number of covariates (columns of X), both continuous and binary.
#' @param p_b Number of binary covariates (0 <= p_b <= p).
#' @param cov_type Character string specifying the covariance function:
#'  - "cov_diag" (independent columns).
#'  - "cov_equi" (equi-correlated columns).
#'  - "cov_ar1" (AR(1)-correlated columns).
#' @param rho Correlation parameter; input to the cov_type function.
#'
#' @details This function simulates a data frame whose rows are multivariate Gaussian with mean zero
#' and covariance structure determined by the "cov_type" argument. Then, `p_b` randomly selected columns are
#' dichotomized with the function \(1(x > 0)\). The continuous columns are of class "numeric",
#' and the binary columns are set to class "factor".
#'
#' @return A simulated data frame with `n` rows and `p` columns, of which `p_b` are binary and `p - p_b` are Gaussian.
#' Each column is either of class "numeric" or "factor".
#' @family generate
#' @export
#'
#' @examples
#' # All columns are continuous:
#' X <- generate_X(n=80, p=100, p_b=0, cov_type="cov_equi", rho=0.5)
#'
#' # Two of the columns are dichotomized (and set to class factor):
#' X <- generate_X(n=80, p=100, p_b=2, cov_type="cov_equi", rho=0.5)
#'
generate_X <- function(n, p, p_b=0, cov_type="cov_ar1", rho=0.5) {

  # Generate covariance matrix based on the specified covariance structure
  sigma_z <- do.call(get(cov_type), args=list(n=n, p=p, rho=rho))

  # Simulate multivariate Gaussian data
  X <- data.frame(matrix(rnorm(n * p), nrow=n) %*% chol(sigma_z))

  # Random indices for the binary columns
  inds_b <- sample(1:p, size=p_b, replace=FALSE)

  # Dichotomize selected columns
  X[, inds_b] <- lapply(X[, inds_b], function(col) as.factor(ifelse(col > 0, 1, 0)))

  # Scale continuous columns
  X <- dplyr::mutate_if(X, is.numeric, ~ scale(.) %>% as.numeric)
  return(X)
}

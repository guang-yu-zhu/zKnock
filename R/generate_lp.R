#' Generate linear predictor with first p_nn beta coefficients = a, all other = 0
#'
#' @param X data.frame with numeric and factor columns only
#' @param p_nn number of non-null covariate predictors.
#' The regression coefficients (beta) corresponding to
#' columns 1:p_nn of x will be non-zero, all other are set to zero.
#' @param a amplitude of non-null regression coefficients
#'
#' @return linear predictor X%*%beta
#' @family generate
#' @export
generate_lp <- function(X, p_nn, a) {

  x <- model.matrix(~., data=X)[,-1]

  n <- nrow(x)
  p <- ncol(x)

  # Standardize design matrix to zero mean, "variance" one
  x_centered <- apply(x, 2, function(x) x - mean(x))
  x <- apply(x_centered, 2, function(x) x / sqrt(mean(x^2)))

  beta = matrix(0, p, 1)
  beta[1:p_nn] = sample(c(-1, 1) * a, p_nn, replace = TRUE)

  lp <- as.numeric(x %*% beta)

  return(lp)

}

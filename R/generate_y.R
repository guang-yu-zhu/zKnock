#' Simulate Gaussian response from a sparse regression model
#'
#' @param X data.frame with numeric and factor columns only.
#' @param p_nn number of non-null covariate predictors.
#' The regression coefficients (beta) corresponding to
#' columns 1:p_nn of x will be non-zero, all other are set to zero.
#' @param a amplitude of non-null regression coefficients
#'
#' @details This function takes as input data.frame X (created with the function \code{generate_X})
#' that may consist of both numeric and binary factor columns. This data frame is then expanded
#' to a model matrix x (with the model.matrix function) and subsequently scaled in the same way as
#' LASSO scaling. Next we simulate y ~ N(x%*%beta,I) where the first p_nn beta coefficients are equal to a, while
#' the remaining coefficients (p_nn+1):ncol(x) are set to zero.
#'
#' @return simulated Gaussian response from regression model y = x%*%beta+epsilon, where epsilon ~ N(0,I) and
#' x is the (scaled) model.matrix of X.
#' @family generate
#' @export
#'
#' @examples
#' set.seed(1)
#' # Simulate 4 Gaussian and 2 binary covariate predictors:
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#' # Simulate response from model y = 2*X[,1] + 2*X[,2] + epsilon, where epsilon ~ N(0,1)
#' y <- generate_y(X, p_nn=2, a=2)
generate_y <- function(X, p_nn, a) {

  x <- model.matrix(~., data=data.frame(X))[,-1]

  n <- nrow(x)
  p <- ncol(x)

  # Standardize design matrix to zero mean, "variance" one
  x_centered <- apply(x, 2, function(x) x - mean(x))
  x <- apply(x_centered, 2, function(x) x / sqrt(mean(x^2)))

  #beta <- rep(c(a, 0), c(p_nn, p - p_nn))
  beta = matrix(0, p, 1)
  beta[1:p_nn] = sample(c(-1, 1) * a, p_nn, replace = TRUE)

  mu <- as.numeric(x %*% beta)

  y <- mu + rnorm(n)

  return(y)

}

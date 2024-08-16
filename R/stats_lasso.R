#' Importance statistics based on the lasso
#' 
#' Fit the lasso path and computes the difference statistic
#'   \deqn{W_j = Z_j - \tilde{Z}_j}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of the 
#' regularization parameter \eqn{\lambda} at which the jth variable 
#' and its knockoff enter the penalized linear regression model, respectively.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. It should be numeric.
#' @param ... additional arguments specific to `glmnet` (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#' 
#' @details This function uses `glmnet` to compute the lasso path
#' on a fine grid of \eqn{\lambda}'s and is a wrapper around the more general
#' [stat.glmnet_lambdadiff].
#' 
#' The `nlambda` parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of `nlambda` is `500`.
#' 
#' Unless a lambda sequence is provided by the user, this function generates it on a 
#' log-linear scale before calling `glmnet` (default 'nlambda': 500).
#' 
#' For a complete list of the available additional arguments, see [glmnet::glmnet()]
#' or [lars::lars()].
#' 
#' @family statistics
#' 
#' @examples
#' # Synthetic Data
#' set.seed(2024)
#' p=200; n=100; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = 1:k
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#'
#' # Knockoff Procedure
#' Xk = create.knockoff(X = X, type = 'shrink', num = 2)
#' res = knockoff.filter(X,y,Xk,statistic = stat.lasso_lambdadiff)
#' res$s
#' 
#' @rdname stat.lasso_lambdadiff
#' @export
stat.lasso_lambdadiff <- function(X, X_k, y, ...) {
  if( is.numeric(y) ){
    y = as.vector(y)
  } else {
    stop('Knockoff statistic stat.lasso_lambdadiff requires the input y to be a numeric vector')
  }
  
  stat.glmnet_lambdadiff(X, X_k, y, family='gaussian', ...)
}

#' Penalized linear regression statistics for knockoff
#' 
#' Computes the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the penalized linear regression model.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. It should be numeric.
#' @param ... additional arguments specific to `glmnet` or `lars` (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#'   
#' @details This function uses `glmnet` to compute the regularization path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The additional `nlambda` 
#' parameter can be used to control the granularity of the grid of \eqn{\lambda} values. 
#' The default value of `nlambda` is `500`.
#' 
#' Unless a lambda sequence is provided by the user, this function generates it on a 
#' log-linear scale before calling `glmnet` (default 'nlambda': 500).
#' 
#' This function is a wrapper around the more general 
#' [stat.glmnet_lambdadiff()].
#' 
#' For a complete list of the available additional arguments, see [glmnet::glmnet()].
#' 
#' @examples
#' # Synthetic Data
#' set.seed(2024)
#' p=200; n=100; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = 1:k
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#'
#' # Knockoff Procedure
#' Xk = create.knockoff(X = X, type = 'shrink', num = 2)
#' res = knockoff.filter(X,y,Xk,statistic = stat.lasso_lambdasmax)
#' res$s
#' 
#' @rdname stat.lasso_lambdasmax
#' @export
stat.lasso_lambdasmax <- function(X, X_k, y, ...) {
  if( is.numeric(y) ){
    y = as.vector(y)
  } else {
    stop('Knockoff statistic stat.lasso_lambdasmax requires the input y to be a numeric vector')
  }

  stat.glmnet_lambdasmax(X, X_k, y, family='gaussian', ...)
}

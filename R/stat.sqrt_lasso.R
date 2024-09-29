#' Importance statistics based on the square-root lasso
#' 
#' Computes the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the SQRT lasso model.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables of numeric type.
#' @param ... additional arguments specific to `slim`.
#' @return A vector of statistics \eqn{W} of length p.
#' 
#' @details With default parameters, this function uses the package `RPtests`
#' to run the SQRT lasso. By specifying the appropriate optional parameters, 
#' one can use different Lasso variants including Dantzig Selector, LAD Lasso,
#' SQRT Lasso and Lq Lasso for estimating high dimensional sparse linear models.
#' 
#' For a complete list of the available additional arguments, see [RPtests::sqrt_lasso()].
#' 
#' @family statistics
#' 
#' @examples
#' set.seed(2024)
#' n=80; p=100; k=10; Ac = 1:k; Ic = (k+1):p
#' X = generate_X(n=n,p=p)
#' y <- generate_y(X, p_nn=k, a=3)
#' Xk = create.shrink_Gaussian(X = X, n_ko = 10)
#' res1 = knockoff.filter(X, y, Xk, statistic = stat.sqrt_lasso,
#'                        offset = 1, fdr = 0.1)
#' res1
#' perf_eval(res1$shat,Ac,Ic)
#' 
#' @rdname stat.sqrt_lasso
#' @importFrom RPtests sqrt_lasso
#' @export
#' @md
stat.sqrt_lasso <- function(X, X_k, y, ...) {
  if (!requireNamespace('RPtests', quietly=T))
    stop('RPtests is not installed', call.=F)
  if (!(is.vector(y) && is.numeric(y)))  {
    stop('Knockoff statistic stat.sqrt_lasso requires the input y to be a numeric vector')
  }
  p = ncol(X)
  
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)
  
  # Compute statistics
  Z = RPtests::sqrt_lasso(cbind(X.swap, Xk.swap), as.numeric(y), ...)
  p = ncol(X)
  orig = 1:p
  W = pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
  W = as.vector(W)
  return(W)
}

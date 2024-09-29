#' Importance statistics based on a GLM
#' 
#' Fits a generalized linear model via penalized maximum likelihood and
#' computes the difference statistic
#'   \deqn{W_j = Z_j - \tilde{Z}_j}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of the 
#' regularization parameter \eqn{\lambda} at which the jth variable 
#' and its knockoff enter the model, respectively.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y Response variable vector of length n. 
#'   Quantitative for family = "gaussian" or "poisson".
#'   - For family = "binomial", y should be either:
#'   	- a two-level factor,
#'   	- a two-column matrix of counts, or
#'   	- proportions.
#'   - For family = "multinomial", y can be a factor with at least two levels or a matrix.
#'   - For family = "cox", y should be a two-column matrix with 'time' and 'status'.
#'   - For family = "mgaussian", y is a matrix of quantitative responses.
#' @param family response type (see above).
#' @param ... additional arguments specific to `glmnet` (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#' 
#' @details This function uses `glmnet` to compute the regularization path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The `nlambda` parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of `nlambda` is `500`.
#' 
#' If the family is 'binomial' and a lambda sequence is not provided by the user, 
#' this function generates it on a log-linear scale before calling 'glmnet'.
#' 
#' The default response family is 'gaussian', for a linear regression model.
#' Different response families (e.g. 'binomial') can be specified by passing an
#' optional parameter 'family'.
#' 
#' For a complete list of the available additional arguments, see [glmnet::glmnet()].
#' 
#' @family statistics
#' 
#' @examples
#' set.seed(2024)
#' n=80; p=100; k=10; Ac = 1:k; Ic = (k+1):p
#' X = generate_X(n=n,p=p)
#' y <- generate_y(X, p_nn=k, a=3)
#' Xk = create.shrink_Gaussian(X = X, n_ko = 10)
#' res1 = knockoff.filter(X, y, Xk, statistic = stat.glmnet_lambdadiff,
#'                        offset = 1, fdr = 0.1)
#' res1
#' perf_eval(res1$shat,Ac,Ic)
#' 
#' @rdname stat.glmnet_lambdadiff
#' @import glmnet
#' @export
#' @md
stat.glmnet_lambdadiff <- function(X, X_k, y, family='gaussian', ...) {
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)
  
  # Compute statistics
  Z = lasso_max_lambda(cbind(X.swap, Xk.swap), y, method='glmnet', family=family, ...)
  p = ncol(X)
  orig = 1:p
  W = Z[orig] - Z[orig+p]
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
  W = as.vector(W)
  return(W)
}

#' GLM statistics for knockoff
#' 
#' Computes the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the generalized linear model.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y Response variable vector of length n. 
#'   Quantitative for family = "gaussian" or "poisson".
#'   - For family = "binomial", y should be either:
#'   	- a two-level factor,
#'   	- a two-column matrix of counts, or
#'   	- proportions.
#'   - For family = "multinomial", y can be a factor with at least two levels or a matrix.
#'   - For family = "cox", y should be a two-column matrix with 'time' and 'status'.
#'   - For family = "mgaussian", y is a matrix of quantitative responses.
#' @param family response type (see above).
#' @param ... additional arguments specific to `glmnet` (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#'   
#' @details This function uses `glmnet` to compute the regularization path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The additional `nlambda` 
#' parameter can be used to control the granularity of the grid of \eqn{\lambda} values. 
#' The default value of `nlambda` is `500`.
#' 
#' If the family is 'binomial' and a lambda sequence is not provided by the user, 
#' this function generates it on a log-linear scale before calling 'glmnet'.
#' 
#' For a complete list of the available additional arguments, see [glmnet::glmnet()].
#' 
#' @examples
#' set.seed(2024)
#' n=80; p=100; k=10; Ac = 1:k; Ic = (k+1):p
#' X = generate_X(n=n,p=p)
#' y <- generate_y(X, p_nn=k, a=3)
#' Xk = create.shrink_Gaussian(X = X, n_ko = 10)
#' res1 = knockoff.filter(X, y, Xk, statistic = stat.glmnet_lambdasmax,
#'                        offset = 1, fdr = 0.1)
#' res1
#' perf_eval(res1$shat,Ac,Ic)
#' 
#' @rdname stat.glmnet_lambdasmax
#' @export
#' @md
stat.glmnet_lambdasmax <- function(X, X_k, y, family='gaussian', ...) {
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)
  
  # Compute statistics
  Z = lasso_max_lambda(cbind(X.swap, Xk.swap), y, method='glmnet', family=family, ...)
  p = ncol(X)
  orig = 1:p
  W = pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
  W = as.vector(W)
  return(W)
}


#' @keywords internal
lasso_max_lambda_lars <- function(X, y, ...) {
  if (!requireNamespace('lars', quietly=T))
    stop('lars is not installed', call.=F)
  
  fit <- lars::lars(X, y, normalize=T, intercept=F, ...)
  lambda <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    entry <- fit$entry[j]
    if (entry > 0) lambda[j] <- fit$lambda[entry]
  }
  return(lambda)
}

#' @keywords internal
lasso_max_lambda_glmnet <- function(X, y, nlambda=500, intercept=T, standardize=T, ...) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  
  # Standardize the variables
  if( standardize ){
    X = scale(X)
  }
    
  n = nrow(X); p = ncol(X)
  if (!methods::hasArg(family) ) family = "gaussian"
  else family = list(...)$family
  
  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      if(!is.numeric(y)) {
        stop('Input y must be numeric.')
      }
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }

  fit <- glmnet::glmnet(X, y, lambda=lambda, intercept=intercept, standardize=F, standardize.response=F, ...)
  
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  if(family=="multinomial") {
      indices <- sapply(fit$beta, function(beta) apply(beta, 1, first_nonzero))
      indices <- apply(indices, 1, min)
  } else {
      indices <- apply(fit$beta, 1, first_nonzero)
  }
  names(indices) <- NULL
  ifelse(is.na(indices), 0, fit$lambda[indices] * n)
}

#' Maximum lambda in lasso model
#' 
#' Computes the earliest (largest) lambda's for which predictors enter the
#' lasso model.
#' 
#' @param X matrix of predictors
#' @param y response vector
#' @param method either 'glmnet' or 'lars'
#' @return vector of maximum lambda's
#' 
#' @keywords internal
lasso_max_lambda <- function(X, y, method=c('glmnet','lars'), ...) {
  switch(match.arg(method), 
         glmnet = lasso_max_lambda_glmnet(X,y,...),
         lars = lasso_max_lambda_lars(X,y,...)
         )
}

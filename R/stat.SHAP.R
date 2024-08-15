#' Importance statistics based on xgboost
#'
#' Computes the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the shap (Hapley Additive exPlanations)
#' of the jth variable and its knockoff, respectively.
#' saabas is a vector of an individualized heuristic feature attribution method, which can be considered as an approximation for shap.
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. If a factor, classification is assumed,
#' otherwise regression is assumed.
#' @param ... additional arguments specific to \code{ranger} (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#'
#' @details .
#'
#' @family statistics
#'
#' @examples
#' # Synthetic Data
#' set.seed(2022)
#' p=200; n=100; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = 1:k
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#'
#' # Knockoff Procedure
#' Xk = create.knockoff(X = X, type = 'shrink', num = 2)
#' res = knockoff.filter(X,y,Xk,statistic = stat.SHAP)
#'
#' @rdname stat.SHAP
#' @export
stat.SHAP<-function(X, X_k, y,nrounds=2){
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)

  x = cbind(X.swap, Xk.swap)
  dtrain <- xgboost::xgb.DMatrix(x, label = y)
  fit <- xgboost::xgb.train(data = dtrain, nrounds = nrounds)

  # predcontrib = TRUE: Requests the SHAP (SHapley Additive exPlanations) values or contributions for the predictions.
  # SHAP values help understand the contribution of each feature to the prediction for each sample.
  # approxcontrib = FALSE: Indicates that exact SHAP values should be computed rather than approximated.
  shap_values <- predict(object = fit, newdata = dtrain, predcontrib = TRUE, approxcontrib = FALSE)
  shap_values <- shap_values[, -c(ncol(shap_values))] # Remove the last column for bias
  mean_abs_shap_values <- colMeans(abs(shap_values))

  #saabas_values <- predict(object = fit, newdata = dtrain, predcontrib = TRUE, approxcontrib = TRUE)
  #saabas_values <- saabas_values[, -c(ncol(saabas_values))] # Remove the last column for bias
  #mean_abs_saabas_values <- colMeans(abs(saabas_values))

  Z <- mean_abs_shap_values
  #Z <- importance.score(fit = fit.model, Y = as.factor(y),X = x)$shap

  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])

  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
  W = as.vector(W)
  return(W)
}

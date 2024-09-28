#' @title Importance statistics based on XGBoost
#'
#' @description
#' Computes the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the SHAP (SHapley Additive exPlanations)
#' of the jth variable and its knockoff, respectively.
#'
#' @details
#'  - In XGBoost, SHAP (SHapley Additive exPlanations) values provide a way to interpret the model's predictions by breaking down the contribution of each feature to the final prediction. SHAP values show how much each feature increases or decreases the prediction compared to the average.
#'  - XGBoost uses the Tree SHAP algorithm, which efficiently computes these values for tree-based models. This helps in understanding both global feature importance (how features influence the model overall) and local explanations (how features impact individual predictions).
#'  - Key benefits include transparency, detailed feature importance, and the ability to explain complex models in a clear, interpretable way.
#'  - saabas is a vector of an individualized heuristic feature attribution method, which can be considered as an approximation for SHAP.
#'
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. If a factor, classification is assumed,
#' otherwise regression is assumed.
#' @param nrounds Number of boosting rounds for training the XGBoost model. Default is 2.
#' @param ... additional arguments specific to `xgboost` (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#'
#' @family statistics
#'
#' @examples
#' # Synthetic Data
#' set.seed(2024)
#' p=30; n=200; k=15
#' mu = rep(0, p); Sigma = diag(p)
#' X = matrix(rnorm(n * p), n)
#' nonzero = 1:k
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#'
#' # Knockoff Procedure
#' Xk = create.knockoff(X = X, type = 'shrink', n_ko = 2)
#' res = knockoff.filter(X, y, Xk, statistic = stat.SHAP)
#' res$shat
#'
#' @rdname stat.SHAP
#' @import xgboost
#' @export
#' @md
stat.SHAP <- function(X, X_k, y, nrounds = 2, ...) {
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X), 1, 0.5)
  swap.M = matrix(swap, nrow = nrow(X), ncol = length(swap), byrow = TRUE)
  X.swap = X * (1 - swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1 - swap.M)

  x = cbind(X.swap, Xk.swap)
  dtrain <- xgboost::xgb.DMatrix(x, label = y)
  fit <- xgboost::xgb.train(data = dtrain, nrounds = nrounds, ...)

  # SHAP values
  shap_values <- predict(object = fit, newdata = dtrain, predcontrib = TRUE, approxcontrib = FALSE)
  shap_values <- shap_values[, -ncol(shap_values)]  # Remove the last column for bias
  mean_abs_shap_values <- colMeans(abs(shap_values))

  Z <- mean_abs_shap_values
  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig + p])

  # Correct for swapping of columns of X and Xk
  W = W * (1 - 2 * swap)
  W = as.vector(W)
  return(W)
}

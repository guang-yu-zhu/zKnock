#' Importance statistics based on xgboost
#'
#' Computes the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the xgboost feature importances
#' of the jth variable and its knockoff, respectively.
#'
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. If a factor, classification is assumed,
#' otherwise regression is assumed.
#' @param family specifies the type of model to be fit: 'gaussian' for regression or 'binomial' for classification.
#' @param nrounds number of boosting rounds for xgboost.
#' @param ... additional arguments specific to xgboost (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr right_join mutate
#' @importFrom tibble tibble
#' @importFrom xgboost xgb.DMatrix xgb.train xgb.importance
#' @rdname stat.xgboost
#' @export
stat.xgboost <- function(X, X_k, y, family = 'gaussian', nrounds = 10, ...) {
  # Randomly swap columns of X and X_k
  swap = rbinom(ncol(X), 1, 0.5)
  swap.M = matrix(swap, nrow = nrow(X), ncol = length(swap), byrow = TRUE)
  X.swap = X * (1 - swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1 - swap.M)

  # Set parameters for xgboost
  if (family == 'gaussian') {
    dtrain <- xgboost::xgb.DMatrix(data = cbind(X.swap, Xk.swap), label = y)
    params <- list(
      booster = "gbtree",
      objective = "reg:squarederror",
      eta = 0.2,
      max_depth = 3,
      eval_metric = "rmse"
    )
  } else if (family == 'binomial') {
    dtrain <- xgboost::xgb.DMatrix(data = cbind(X.swap, Xk.swap), label = as.numeric(y) - 1)
    params <- list(
      booster = "gbtree",
      objective = "binary:logistic",
      eta = 0.2,
      max_depth = 3,
      eval_metric = "logloss"
    )
  }

  # Train the xgboost model
  xgb_model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )

  # Extract feature names
  X.names <- colnames(X)
  p <- ncol(X)
  if (is.null(X.names)) X.names <- 1:p

  # Get feature importance
  all_features <- c(X.names, paste0('knock-', X.names))
  importance_matrix <- xgboost::xgb.importance(feature_names = all_features, model = xgb_model)

  # Ensure all features are included with zero importance for unused ones
  importance_full <- importance_matrix %>%
    dplyr::right_join(tibble::tibble(Feature = all_features), by = "Feature") %>%
    dplyr::mutate(
      Gain = ifelse(is.na(Gain), 0, Gain),
      Cover = ifelse(is.na(Cover), 0, Cover),
      Frequency = ifelse(is.na(Frequency), 0, Frequency)
    )

  # Compute statistics
  Z <- importance_full$Gain
  W <- abs(Z[1:p]) - abs(Z[(p + 1):(2 * p)])

  # Correct for swapping of columns of X and X_k
  W <- W * (1 - 2 * swap)
  return(as.vector(W))
}

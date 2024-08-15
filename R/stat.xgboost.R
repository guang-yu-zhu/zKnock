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
#' res = knockoff.filter(X,y,Xk,statistic = stat.xgboost,family='gaussian')
#'
#' @rdname stat.xgboost
#' @export
stat.xgboost<-function(X, X_k, y,family='gaussian',nrounds=10){
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)


  if(family=='gaussian'){
    dtrain <- xgboost::xgb.DMatrix(data = cbind(X.swap, Xk.swap), label = y)
    params <- list(
      booster = "gbtree",
      objective = "reg:squarederror",  # Regression with squared error (default for regression)
      eta = 0.2,                       # Learning rate
      max_depth = 3,                   # Maximum tree depth
      eval_metric = "rmse"             # Root Mean Squared Error (alternative: "mae" for Mean Absolute Error)
    )
  }else if(familly=='binomial'){
    dtrain <- xgboost::xgb.DMatrix(data = cbind(X.swap, Xk.swap), label = as.numeric(y)-1)
    params <- list(
      booster = "gbtree",
      objective = "binary:logistic",  # Logistic regression
      eta = 0.2,                      # Learning rate
      max_depth = 3,                  # Maximum tree depth
      eval_metric = "logloss"         # Evaluation metric for binary classification
    )
  }

  xgb_model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,                  # Number of boosting rounds
    verbose = 0                     # Print the training log
  )

  X.names = colnames(X)
  ncol(X)
  if (is.null(X.names)) X.names=1:p

  all_features <- c(X.names,paste0('knock-',X.names))
  importance_matrix <- xgboost::xgb.importance(feature_names =all_features,model = xgb_model)
  # Ensure all features are included, with zero importance for those not used
  importance_full <- importance_matrix %>%
    right_join(tibble(Feature = all_features), by = "Feature") %>%
    mutate(Gain = ifelse(is.na(Gain), 0, Gain),
           Cover = ifelse(is.na(Cover), 0, Cover),
           Frequency = ifelse(is.na(Frequency), 0, Frequency))
  # Compute statistics
  Z=importance_full$Gain
  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])

  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
  W = as.vector(W)
  return(W)
}

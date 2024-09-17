#' Importance statistics based on a GLM with cross-validation
#'
#' Fits a generalized linear model via penalized maximum likelihood with cross-validation
#' and computes the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the
#' jth variable and its knockoff, respectively. The regularization parameter
#' \eqn{\lambda} is selected by cross-validation and computed with `glmnet`.
#'
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y Response variable vector of length n. Quantitative for family="gaussian" or "poisson".
#' For family="binomial", y should be either a two-level factor, a two-column matrix of counts,
#' or proportions. For family="multinomial", y can be a factor with at least two levels or a matrix.
#' For family="cox", y should be a two-column matrix with 'time' and 'status'. For family="mgaussian",
#' y is a matrix of quantitative responses.
#' @param family Response type, one of 'gaussian', 'binomial', 'multinomial', 'cox', or 'mgaussian'.
#' @param cores Number of cores to use for parallel computation. Defaults to 2 if available.
#' @param ... Additional arguments specific to `glmnet` (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#'
#' @details This function uses the `glmnet` package to fit a GLM model via penalized maximum likelihood.
#' The value of \eqn{\lambda} is chosen by 10-fold cross-validation unless specified otherwise.
#'
#' The optional `nlambda` parameter can control the number of \eqn{\lambda} values, and the default
#' is 500. For family="binomial", if the lambda sequence isn't provided, a log-linear sequence is generated.
#'
#' For more details, see [glmnet::cv.glmnet()] and [glmnet::glmnet()].
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
#' res = knockoff.filter(X, y, Xk, statistic = stat.glmnet_coefdiff, family='gaussian')
#' res$shat
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @rdname stat.glmnet_coefdiff
#' @export
stat.glmnet_coefdiff <- function(X, X_k, y, family = 'gaussian', cores = 2, ...) {
  # Ensure glmnet is installed
  if (!requireNamespace('glmnet', quietly = TRUE)) {
    stop('glmnet is not installed', call. = FALSE)
  }

  # Set up parallelization
  parallel <- TRUE
  if (!requireNamespace('doParallel', quietly = TRUE)) {
    warning('doParallel is not installed. Statistics will be slower without parallelization.', call. = FALSE, immediate. = TRUE)
    parallel <- FALSE
  }
  if (!requireNamespace('parallel', quietly = TRUE)) {
    warning('parallel is not installed. Statistics will be slower without parallelization.', call. = FALSE, immediate. = TRUE)
    parallel <- FALSE
  }

  # Detect available cores and adjust if necessary
  if (parallel) {
    ncores <- parallel::detectCores(all.tests = TRUE, logical = TRUE)
    cores <- min(cores, ncores)
    if (cores > 1) {
      doParallel::registerDoParallel(cores = cores)
    } else {
      parallel <- FALSE
    }
  }

  # Randomly swap columns of X and X_k
  swap <- rbinom(ncol(X), 1, 0.5)
  swap.M <- matrix(swap, nrow = nrow(X), ncol = length(swap), byrow = TRUE)
  X.swap <- X * (1 - swap.M) + X_k * swap.M
  Xk.swap <- X * swap.M + X_k * (1 - swap.M)

  p <- ncol(X)

  # Compute GLM coefficients using cross-validation
  glmnet.coefs <- cv_coeffs_glmnet(cbind(X.swap, Xk.swap), y, family = family, parallel = parallel, ...)

  # Extract coefficients and compute statistics
  if (family == "multinomial") {
    Z <- abs(glmnet.coefs[[1]][2:(2 * p + 1)])
    for (b in 2:length(glmnet.coefs)) {
      Z <- Z + abs(glmnet.coefs[[b]][2:(2 * p + 1)])
    }
  } else if (family == "cox") {
    Z <- glmnet.coefs[1:(2 * p)]
  } else {
    Z <- glmnet.coefs[2:(2 * p + 1)]
  }

  # Calculate the difference statistics W
  W <- abs(Z[1:p]) - abs(Z[(p + 1):(2 * p)])
  W <- W * (1 - 2 * swap)

  # Stop parallel cluster if applicable
  if (parallel && cores > 1) {
    doParallel::stopImplicitCluster()
  }

  return(as.vector(W))
}


#' @keywords internal
cv_coeffs_glmnet <- function(X, y, nlambda=500, intercept=T, parallel=T, ...) {
  # Standardize variables
  X = scale(X)

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

  cv.glmnet.fit <- glmnet::cv.glmnet(X, y, lambda=lambda, intercept=intercept,
                                     standardize=F,standardize.response=F, parallel=parallel, ...)

  coef(cv.glmnet.fit, s = "lambda.min")
}

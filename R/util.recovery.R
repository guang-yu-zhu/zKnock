
#' Simulate from glmnet penalized regression model
#'
#' @param y response vector (either "numeric" or "factor") that gets passed to cv.glmnet
#' @param X data.frame of covariates that are passed to cv.glmnet
#' @param ... other parameters passed to the function cv.glmnet
#'
#' @return simulated response
#'
#' @examples
#' set.seed(1)
#'
#' X = data.frame(matrix(rnorm(100 * 20), 100, 20))
#' y = X[,1] + rnorm(100)
#'
#' # simulate from elastic-net regression:
#' ysim = glmnet.recovery(y=y, X=X)
#' @family recovery
#' @export
#' @md
glmnet.recovery <- function(y, X, ...) {

  x <- model.matrix(~., data = X)[,-1]

  if (is.factor(y)) {

    classes <- levels(y)

    K <- length(classes)

    gm.cv <- glmnet::cv.glmnet(y=y, x=x, family="multinomial", intercept=TRUE, alpha=1, ...)

    # Beta coefficients (excluding intercept)
    beta.coefs <- as.numeric(coef(gm.cv, s = "lambda.1se")[[2]])[-1]

    mu <- predict(gm.cv, newx=x, type="response", s="lambda.1se")

    mat.multinom <- apply(mu, 1, function(prob) rmultinom(n=1, size=1, prob=prob))

    y.sim <- classes[apply((1:K)*mat.multinom, 2, max)]

    y.sim <- factor(y.sim, levels=classes)

    rmse <- NULL

  } else {

    if(!is.numeric(y)) stop("class(y) needs to be either 'numeric' or 'factor'")

    gm.cv <- glmnet::cv.glmnet(y=y, x=x, family="gaussian", intercept=TRUE, alpha=1, ...)

    # Beta coefficients (excluding intercept)
    beta.coefs <- as.numeric(coef(gm.cv, s = "lambda.1se"))[-1]

    # columns of predictor matrix corresponding to non-zero beta.coefs:
    non.zero.cols <- which(beta.coefs != 0)

    # Total number of non-zero parameters (including intercept, hence + 1)
    s.lambda = length(non.zero.cols) + 1

    mu <- predict(gm.cv, newx=x, type="response", s="lambda.1se")

    rmse = sqrt(sum((y-mu)^2)/(length(y) - s.lambda))

    y.sim <- rnorm(n=length(y), mean=mu, sd=rmse)

  }

  return(y.sim)

}



#' Simple knockoff generator
#'
#' This function generates knockoff variables using the following approaches:
#' - Continuous variables: Ordinary least squares (OLS) regression.
#' - Factor variables: Multinomial logistic regression.
#' If \code{X} is empty, knockoffs are sampled from the marginal distribution of \code{y}.
#'
#' @param y A response vector that can be either "numeric" or "factor". The vector is passed to \code{cv.glmnet} for fitting the model.
#' @param X A data frame of covariates that are passed to \code{cv.glmnet} for model fitting.
#'
#' @return A simulated response vector \code{y.sim} with the same type and length as \code{y}.
#' @importFrom nnet multinom
#' @family recovery
#' @export
#' @md
simple.recovery <- function(y, X) {

  dataset <- cbind(y, X)

  if (is.factor(y)) {

    # Handle factor response with multinomial logistic regression
    classes <- levels(y)

    if (length(classes) < 2) stop("A categorical feature with less than two different levels was provided. This feature is uninformative and should thus be removed.")

    # Fit multinomial model
    fit <- nnet::multinom(y ~ ., data = dataset)
    mu <- predict(fit, type = "probs")

    # Handle binary factor case separately
    if (length(classes) == 2) {
      # Complement probabilities if y is binary
      mu <- matrix(c(mu, 1 - mu), ncol = 2)
    }

    # Simulate from the fitted multinomial distribution
    mat.multinom <- apply(mu, 1, function(prob) rmultinom(n = 1, size = 1, prob = prob))
    y.sim <- factor(classes[apply((1:length(classes)) * mat.multinom, 2, max)], levels = classes)

  } else {
    # Handle numeric response with ordinary least squares (OLS)
    fit <- lm(y ~ ., data = dataset)
    mu <- predict(fit, type = "response")
    sigma <- summary(fit)$sigma
    # Simulate knockoffs from normal distribution with estimated parameters
    y.sim <- rnorm(n = length(y), mean = mu, sd = sigma)
  }

  return(y.sim)
}


#' Calculate \eqn{\hat{X}} by fitting PLS regression on its neighbours
#'
#' This function fits a partial least squares (PLS) or sparse partial least squares (sPLS) regression model to predict \eqn{\hat{Y}} using the predictor matrix \eqn{X}. The fitted values \eqn{\hat{Y}} are returned.
#'
#' @param Y A numeric vector representing the response variable.
#' @param X A numeric matrix or data frame representing the predictor variables.
#' @param ncomp An integer specifying the number of components to use in the PLS regression model.
#' @param keepX A numeric vector specifying the number of variables to keep in each dimension. By default, all variables are used.
#'
#' @return A numeric vector \eqn{\hat{Y}} containing the predicted values from the PLS regression model.
#'
#' @details This function performs PLS regression using the \code{mixOmics} package. It allows for both standard PLS and sparse PLS, where \code{keepX} controls the number of variables retained in each dimension. The function returns the predicted values \eqn{\hat{Y}} for the response variable \eqn{Y} using all the components specified by \code{ncomp}.
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' Y <- rnorm(10)
#' Y.hat <- pls.recovery(Y, X, ncomp = 2, keepX = rep(5, 2))
#'
#' @rdname pls.recovery
#' @family recovery
#' @export
#' @md
pls.recovery <- function(Y, X, ncomp, keepX = rep(ncol(X))) {
  X <- as.data.frame(X)
  n <- nrow(X)
  p <- ncol(X)

  # Fit the sparse PLS model using mixOmics package
  pls <- mixOmics::spls(X, Y, ncomp = ncomp, scale = F, keepX = keepX, mode = "regression")

  # Extract predictions using all ncomp components
  Y.hat <- predict(pls, X)$predict[1:n, 1, ncomp]

  return(Y.hat)
}


#' Calculate \eqn{\hat{X}} by fitting sparse PLS
#'
#' This function fits a sparse partial least squares (sPLS) model to the response variable \eqn{Y} using the predictor matrix \eqn{X}. The fitted values \eqn{\hat{Y}} are returned.
#'
#' @param Y A numeric vector representing the response variable.
#' @param X A numeric matrix or data frame representing the predictor variables.
#' @param ncomp An integer specifying the number of components to include in the sPLS model.
#' @param eta A numeric value between 0 and 1 that controls the sparsity of the PLS loadings.
#'
#' @return A numeric vector \eqn{\hat{Y}} containing the fitted values from the sparse PLS model.
#'
#' @details This function fits a sparse PLS model using the \code{spls} package. The sparsity level is controlled by the \code{eta} parameter, where higher values lead to sparser loadings. The number of components is specified by \code{ncomp}. The function returns the fitted values of \eqn{Y} from the model.
#'
#' @importFrom spls spls
#' @rdname spls.recovery
#' @family recovery
#' @export
#' @md
spls.recovery <- function(Y, X, ncomp, eta) {
  X <- as.data.frame(X)

  # Fit the sparse PLS model
  fit <- spls::spls(X, Y, K = ncomp, eta = eta)

  # Predict the values of Y using the fitted model
  Y.hat <- predict(fit)

  return(Y.hat)
}



#' Calculate \eqn{\hat{X}} by fitting OLS regression on its neighbors
#'
#' This function performs ordinary least squares (OLS) regression, estimating the coefficients that best fit the response variable \eqn{Y} given the predictor matrix \eqn{X}. The fitted values \eqn{\hat{Y}} are then returned.
#'
#' @param Y A numeric vector representing the response variable.
#' @param X A numeric matrix or data frame representing the predictor variables. The function removes any duplicated columns from \code{X}.
#'
#' @return A numeric vector \eqn{\hat{Y}} containing the fitted values from the OLS regression.
#'
#' @details The function first converts the input matrix \code{X} to a data frame to remove any duplicate columns, ensuring that the matrix used for regression is unique in terms of predictors. It then applies OLS regression by solving the normal equations to estimate the coefficients \eqn{\beta}. The fitted values \eqn{\hat{Y}} are computed as \eqn{X \beta}.
#'
#' @family recovery
#' @export
#' @md
ols.recovery <- function(Y, X) {
  X <- as.data.frame(X)
  X <- X[!duplicated(as.list(X))]
  X <- as.matrix(X)
  betas <- solve(crossprod(X)) %*% t(X) %*% Y
  Y.hat <- X %*% betas
  return(Y.hat)
}






#' Calculate the Empirical Number of Components for PLS Regression
#'
#' This function calculates the empirical number of components used for PLS regression
#' by applying the \eqn{PC_p1} criterion from Bai and Ng (2002). It determines the optimal number of factors based on the minimum value of the criterion. The method and function is adapted from Y Fan et al. (2019).
#'
#' @param X A numeric matrix or data frame where rows represent observations and columns represent variables.
#' @param rmax An integer specifying the maximum number of factors to consider. Default is 10.
#'
#' @return An integer representing the optimal number of components based on the \eqn{PC_p1} criterion.
#'
#' @references
#' - Fan Y, Lv J, Sharifvaghefi M et al. IPAD: Stable Interpretable Forecasting with Knockoffs Inference. Journal of the American Statistical Association 2020;115:1822–34.
#' - Bai J, Ng S. Determining the Number of Factors in Approximate Factor Models. Econometrica 2002;70:191–221.
#' @keywords internal
#' @md
r_criterion <- function(X, rmax = 10){
  # Number of observations (n) and variables (p)
  n = nrow(X)
  p = ncol(X)

  # Initialize a matrix to store the PC_p1 criterion values for each possible number of components
  PC = matrix(NA,rmax+1,1)

  SX = scale(X)
  mXX = SX %*%t(SX)

  # Calculate the criteria through the number of components from rmax to 0
  for(k in rmax:0){
    if (k == 0){
      PC[k+1,] = sum(SX^2/(n*p))
    } else
      eig <- RSpectra::eigs_sym(mXX, k, which = "LM", sigma = NULL, lower = TRUE)
    meigvec <- eig$vectors
    mF <- sqrt(n) * meigvec       # estimated factors
    Lam <- (t(mF) %*% SX)/n
    if(k==rmax){
      sigma2 = sum((SX-mF %*% Lam)^2)/(n*p)
    }
    PC[k+1,] = sum((SX-mF %*% Lam)^2)/(n*p) + k*sigma2*((n+p)/(n*p))*log((n*p)/(n+p))
  }
  r = which(PC == min(PC))[1]-1
}

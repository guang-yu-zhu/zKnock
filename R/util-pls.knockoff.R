#' Calculate \eqn{\hat{X}} by fitting PLS regression on its neighbours
#'
#' @rdname pls.recovery.generator
#' @keywords internal
#'
pls.recovery.generator <- function(Y, X, ncomp, keepX = rep(ncol(X))){
  X <- as.data.frame(X)
  #X <- X[!duplicated(as.list(X))]
  n <- nrow(X)
  p <- ncol(X)
  pls <- mixOmics::spls(X, Y, ncomp = ncomp, scale = F, keepX = keepX, mode = "regression")
  Y.hat <- predict(pls, X)$predict[1:n, 1, ncomp] # the prediction use all ncomp
  return(Y.hat)
}

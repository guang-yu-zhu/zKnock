#' @title Generate Knockoff Matrix
#' @description Generate different types of knockoff matrices given an original one.
#'
#' @param X A numeric matrix representing the original design matrix.
#' @param type The knockoff type to be generated. Available options for `type` are:
#' 1. "shrink" for the shrink Gaussian knockoff;
#' 2. "sparse" for the sparse Gaussian knockoff;
#' 3. "pc" for the principal component (PC) knockoff;
#' 4. "pls" for the partial least squares (PLS) knockoff;
#' 5. "zpls" for the GZ's sparse partial least squares knockoff.
#' @param num A positive integer indicating the number of knockoff matrices to be created.
#' @param ncomp The number of principal components to be used for generating knockoff matrices. Default is 10.
#' @param verbose Logical; if TRUE, messages about the progress are printed. Default is FALSE.
#' @param ... Additional arguments to be passed to specific knockoff creation functions.
#'
#' @return A list of created knockoff matrices. Each element of the list corresponds to a generated knockoff matrix.
#' @family create
#'
#' @importFrom spcov spcov
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(100), nrow = 10)
#' Xk <- create.knockoff(X = X, type = "shrink", num = 5)
#'
#' @export
#' @md
create.knockoff <- function(X, type, num, ncomp = 10, verbose = FALSE, ...) {
  # Input validation
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be a numeric matrix.")
  if (!is.numeric(num) || length(num) != 1 || num <= 0) stop("num must be a positive integer.")

  # Initialize result list to store knockoff matrices
  result <- vector(mode = "list", length = num)

  # Type "shrink" for shrink Gaussian knockoffs
  if (type == "shrink") {
    for (i in 1:num) {
      if (verbose) cat('--Generate', i, 'knockoffs\n')
      result[[i]] <- create.second_order(X = X, method = "sdp")
    }

    # Type "sparse" for sparse Gaussian knockoffs
  } else if (type == "sparse") {
    mu <- colMeans(X)
    S <- stats::cov(X)
    p <- ncol(X)
    step.size <- 100
    P <- matrix(1, p, p)
    diag(P) <- 0
    lam <- 0.06
    mm <- spcov::spcov(Sigma = diag(diag(S)), S = (S + 0.1 * diag(1, p)), lambda = lam * P, step.size = step.size)

    for (i in 1:num) {
      if (verbose) cat('--Generate', i, 'knockoffs\n')
      result[[i]] <- create.gaussian(X = X, mu = mu, Sigma = mm$Sigma, method = "sdp")
    }

    # Type "pc" for principal component knockoffs
  } else if (type == "pc") {
    for (i in 1:num) {
      if (verbose) cat('--Generate', i, 'knockoffs\n')
      result[[i]] <- create.pc.knockoff(X = X, pc.num = ncomp)
    }

    # Type "pls" for partial least squares knockoffs
  } else if (type == "pls") {
    for (i in 1:num) {
      if (verbose) cat('--Generate', i, 'knockoffs\n')
      result[[i]] <- create.pls.knockoff(X = X, ncomp = ncomp)
    }
  } else if (type == "zpls") {
    for (i in 1:num) {
      if (verbose) cat('--Generate', i, 'knockoffs\n')
      result[[i]] <- create.zpls.knockoff(X = X, ncomp = ncomp, ...)
    }
  } else {
    stop("Invalid knockoff type!")
  }

  return(result)
}

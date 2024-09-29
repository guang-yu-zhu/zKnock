#' Generate knockoff variables using PLS regression (PLSKO)
#'
#' This function generates knockoff variables using Partial Least Squares (PLS) regression, following the PLSKO algorithm. It is useful for generating knockoff variables for high-dimensional data.
#'
#' @param X A numeric matrix or data frame. The original design data matrix with \eqn{n} observations as rows and \eqn{p} variables as columns.
#' @param n_ko An integer specifying the number of knockoff variables to generate. Default is 1.
#' @param nb.list Optional. A list of length \eqn{p} or an adjacency matrix of \eqn{p \times p} that defines the neighbor relationships among variables.
#'   - A list of length \eqn{p} should include the neighbors' indices of each variable from \eqn{X_1} to \eqn{X_p} in order. The \eqn{i^{th}} element in the list includes the indices of the neighbor variables of \eqn{X_i}, or \code{NULL} when no neighbors.
#'   - An adjacency matrix should be symmetric with binary elements. \eqn{M_{ij} = 1} indicates that \eqn{X_i} and \eqn{X_j} are neighbors; \eqn{M_{ij} = 0} indicates no neighbor relationship or diagonal entries.
#'   - If not provided or \code{NULL}, neighborhoods are determined based on correlations.
#' @param threshold.abs Optional. A value between \eqn{0} and \eqn{1} to specify an absolute correlation threshold for defining neighborhoods.
#' @param threshold.q Optional. A numeric value between 0 and 1 indicating the quantile of the correlation values to use as a threshold for defining neighborhoods. Default is 0.9.
#' @param ncomp Optional. An integer specifying the number of components to use in the PLS regression. Default is \code{NULL}, in which case the number of components is chosen empirically.
#' @param sparsity Optional. A numeric value between 0 and 1 specifying the sparsity level in the PLS regression. Default is 1 (no sparsity).
#' @param verbose Logical. Whether to display progress information during the knockoff generation. Default is TRUE.
#'
#' @return A list of generated knockoff matrices, where each matrix has \eqn{n} rows (observations) and \eqn{p} columns (variables).
#'
#' @family create
#' @references
#' Yang, Guannan, et al. "PLSKO: a robust knockoff generator to control false discovery rate in omics variable selection." bioRxiv (2024): 2024-08.
#'
#' @importFrom mixOmics spls
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(100), nrow = 10)
#' Xk <- create.pls(X = X, ncomp = 3)
#'
#' @export
#' @md
create.pls <- function(X, n_ko = 1, ncomp = NULL, sparsity = 1, nb.list = NULL, threshold.abs = NULL, threshold.q = 0.9, verbose = FALSE) {
  #n_ko = 1; ncomp = NULL; sparsity = 1; nb.list = NULL; threshold.abs = NULL; threshold.q = 0.9; verbose = FALSE;
  # Initialize dimensions
  n <- nrow(X)
  p <- ncol(X)

  # Convert data frame to matrix if necessary
  if (is.data.frame(X)) {
    X.names <- names(X)
    X <- as.matrix(X)
  } else if (is.matrix(X)) {
    X.names <- colnames(X)
  } else {
    stop('Input X must be a numeric matrix or data frame')
  }

  # Validate input
  if (!is.numeric(n_ko) || length(n_ko) != 1 || n_ko <= 0 || floor(n_ko) != n_ko) stop("n_ko must be a positive integer.")
  if (!is.numeric(ncomp) || length(ncomp) != 1 || ncomp <= 0 || floor(ncomp) != ncomp) stop("ncomp must be a positive integer.")
  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')
  # Assign default column names if necessary
  if (is.null(X.names)) {
    X.names <- paste0('X', 1:p)
    colnames(X) <- X.names
  }
  # initial X_k
  X_k <- matrix(NA, nrow = n, ncol = p)
  rownames(X_k) <- rownames(X)
  colnames(X_k) <- paste0(X.names, ".tilde")


  # Neighborhood determination based on correlations or user input
  if (is.null(nb.list)) {
    sample.order.mat <- diag(1, p, p)[, sample(p)] # Randomly reorder columns
    X <- X %*% sample.order.mat
    mu <- colMeans(X)
    X <- scale(X, center = TRUE, scale = FALSE)
    cor_mat <- cor(X)
    diag(cor_mat) <- 0
    threshold <- if (!is.null(threshold.abs)) threshold.abs else quantile(unlist(abs(cor_mat)), prob = threshold.q, names = FALSE)

    # Define neighborhoods based on correlation threshold
    neighborhoods <- vector("list", length = p)
    for (i in 1:p) {
      cols <- which(abs(cor_mat[i, ]) > threshold)
      if (length(cols) > 0) neighborhoods[[i]] <- cols
    }
  } else {
    # User-provided neighborhoods
    if (is.list(nb.list)) {
      neighborhoods <- nb.list
    } else if (is.matrix(nb.list) & isSymmetric(nb.list)) {
      neighborhoods <- vector("list", length = p)
      for (i in 1:p) {
        cols <- which(nb.list[i, ] == 1)
        if (length(cols) > 0) neighborhoods[[i]] <- cols
      }
    }
  }

  # Empirical determination of ncomp if not provided
  if (is.null(ncomp)) r_emp <- r_criterion(X)

  # Initialize output
  result <- vector("list", length = n_ko)
  for (ko_index in 1:n_ko) {
    if (verbose) cat('--Generate', ko_index, 'knockoffs\n')

    for (i in 1:p) {
      nb <- neighborhoods[[i]]
      k.nb <- nb[nb < i]
      X.nb <- X[, nb]
      X_k.nb <- X_k[, k.nb]

      if (length(nb) == 0) {
        # No neighbors, zero regression fit
        Y <- X[, i]
        Y.hat <- 0
      } else if (length(nb) == 1) {
        # Ordinary least squares with single neighbor
        X.run <- cbind(X.nb, X_k.nb)
        Y <- X[, i]
        Y.hat <- ols.recovery(Y, X.run)
      } else {
        # PLS regression with multiple neighbors
        X.run <- cbind(X.nb, X_k.nb)
        Y <- X[, i]
        if (is.null(ncomp)) ncomp <- ceiling(min(ncol(X.nb) / 2, r_emp))
        this.ncomp <- min(ncomp, ncol(X.run) - 1)
        keepX <- rep(round(sparsity * ncol(X.run)), this.ncomp)
        Y.hat <- pls.recovery(Y, X.run, ncomp = this.ncomp, keepX = keepX)
      }

      # Generate knockoffs
      Y.res <- Y - Y.hat
      res <- sample(Y.res)
      X_k[, i] <- Y.hat + res
    }

    # Reintroduce the mean
    X_k <- apply(X_k, 1, function(x) { x + mu })
    X_k <- t(X_k)

    # Restore the original column order if neighborhoods were generated by correlation
    if (is.null(nb.list)) X_k <- X_k %*% t(sample.order.mat)

    result[[ko_index]] <- X_k
  }

  return(result)
}

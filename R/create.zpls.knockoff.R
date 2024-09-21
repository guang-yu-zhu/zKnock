#' Generate a knockoff variable set with spls regression
#' @details
#' Neighborhood Generation:
#'  - If `threshold.abs` is given: That absolute value is used directly.
#'  - If `threshold.q` is given: The threshold is set based on the quantile of the absolute correlation values.
#'  - If neither is provided: The function defaults to the 90th percentile of the absolute correlation values, which corresponds to using the strongest 10% of correlations to define neighborhoods.
#' Calculate the fitted value of \eqn{X_j} with  `spls::spls`, see [spls::spls()] for more details.
#'
#' @param X A numeric matrix or data frame. The original design data matrix with \eqn{n} observations as rows and \eqn{p} variables as columns.
#' @param nb.list Optional. A list of length \eqn{p} or adjacency matrix of \eqn{p \times p} that defines the neighbourship of variables.
#' A list of length \eqn{p} should include the neighbours' index of each variable from \eqn{X_1} to \eqn{X_p} in order;
#' The \eqn{i^{th}} element in the list includes the indices of the neighbour variables of \eqn{X_i}, or \code{NULL} when no neighbours.
#' A adjacency matrix should be symmetric with only binary element and where \eqn{M_{ij} = 1} when \eqn{X_i} and \eqn{X_j} are neighbours;
#' otherwise \eqn{M_{ij} = 0} when not neighbour or on diagonal (i.e. \eqn{i = j}).
#' If not provided or NULL, the neighborhoods are determined based on correlations.
#' @param threshold.abs Optional. A value between \eqn{0} and \eqn{1}. A numeric value specifying an absolute correlation threshold to define neighborhoods.
#' @param threshold.q Optional. A numeric value between 0 and 1 indicating the quantile of the correlation values to use as a threshold. Default is 0.9.
#' @param ncomp Optional. An integer specifying the number of components to use in the PLS regression. Default is 2.
#' @param sparsity Optional. A numeric value between 0 and 1 specifying the sparsity level in the PLS regression. Default is 1 (no sparsity).
#'
#' @return A matrix of generated knockoff variables of \eqn{n \times p}.
#'
#' @family create
#' @references
#' Yang, Guannan, et al. "PLSKO: a robust knockoff generator to control false discovery rate in omics variable selection." bioRxiv (2024): 2024-08.
#'
#' @import spls
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(100), nrow = 10)
#' Xk <- create.zpls.knockoff(X = X, ncomp = 3,eta=0.3)
#' @export
#' @md
create.zpls.knockoff <- function(X,ncomp = NULL, eta=0.3, nb.list = NULL, threshold.abs = NULL, threshold.q = 0.9) {
  #nb.list = NULL; threshold.abs = NULL; threshold.q = 0.9; ncomp = 4; sparsity = 1
  n <- nrow(X)
  p <- ncol(X)

  # Input type validation
  if (is.data.frame(X)) {
    X.name <- names(X)
    X <- as.matrix(X)
  } else if (is.matrix(X)) {
    X.name <- colnames(X)
  } else {
    stop('Input X must be a numeric matrix or data frame')
  }

  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')

  # If nb.list is not provided, generate a neighborhood list based on correlations
  if (is.null(nb.list)) {
    # Randomly swap columns of X
    sample.order.mat <- diag(1, p, p)[, sample(p)]
    X <- X %*% sample.order.mat

    mu <- colMeans(X)
    X <- scale(X, center = TRUE, scale = FALSE)

    # Calculate the correlation matrix and replace diagonal with 0
    cor_mat <- cor(X)
    diag(cor_mat) <- 0

    # Determine the correlation threshold
    if (!is.null(threshold.abs)) {
      threshold <- threshold.abs
    } else{
      threshold <- quantile(unlist(abs(cor_mat)), prob = threshold.q, names = FALSE)
    }

    # Create a list of neighborhoods based on the correlation threshold
    neighborhoods <- vector(mode = "list", length = p)
    for (i in 1:p) {
      cols <- which(abs(cor_mat[i, ]) > threshold)
      if (length(cols) > 0) {
        neighborhoods[[i]] <- cols
      }
    }
  } else {
    # If a neighborhood list is provided by the user
    if (is.list(nb.list)) {
      neighborhoods <- nb.list
    } else if (is.matrix(nb.list) & isSymmetric(nb.list)) {
      neighborhoods <- vector(mode = "list", length = p)
      for (i in 1:p) {
        cols <- which(nb.list[i, ] == 1)
        if (length(cols) > 0) {
          neighborhoods[[i]] <- cols
        }
      }
    }
  }

  # Create knockoff matrix
  X_k <- matrix(NA, nrow = n, ncol = p)
  X_k <- as.data.frame(X_k)
  rownames(X_k) <- rownames(X)
  colnames(X_k) <- paste0(colnames(X), "k")

  # Determine empirical number of components if ncomp is not provided
  if (is.null(ncomp)) {
    r_emp <- r_criterion(X)
  }

  # Generate knockoff variables
  for (i in 1:p) {
    nb <- neighborhoods[[i]]
    k.nb <- neighborhoods[[i]][neighborhoods[[i]] < i]
    X.nb <- X[, nb]
    X_k.nb <- X_k[, k.nb]

    if (length(nb) == 0) {
      Y <- X[, i]
      Y.hat <- 0
    } else if (length(nb) == 1) {
      X.run <- cbind(X.nb, X_k.nb)
      Y <- X[, i]
      Y.hat <- linear.regression.generator(Y, X.run)
    } else {
      X.run <- cbind(X.nb, X_k.nb)
      Y <- X[, i]
      if (is.null(ncomp)) {
        ncomp <- ceiling(min(ncol(X.nb) / 2, r_emp))
      }
      this.ncomp <- min(ncomp, ncol(X.run) - 1)
      #keepX <- rep(round(sparsity * ncol(X.run)), this.ncomp)
      #cat('this.ncomp:',this.ncomp,',  eta:',eta,'\n')
      Y.hat <- spls.recovery.generator(Y, X.run, ncomp = this.ncomp, eta=eta)
    }

    # Calculate residuals and permute
    Y.res <- Y - Y.hat
    res <- sample(Y.res)
    X_k[, i] <- Y.hat + res
  }

  # Add the mean back to knockoff variables generated from the centered data
  X_k <- apply(X_k, 1, function(x) {x + mu})
  X_k <- t(X_k)

  # Swap columns back to original order if necessary
  if (is.null(nb.list)) {
    X_k <- X_k %*% t(sample.order.mat)
  }

  return(X_k)
}

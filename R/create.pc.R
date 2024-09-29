#' Create PC Knockoffs
#'
#' A sequential algorithm to create non-parametric knockoffs based on principal component regression and residual permutation.
#'
#' @details
#' For each original variable \eqn{\mathbf{x}_j}, where \eqn{j = 1, \ldots, p}, the following steps are performed to generate knockoff variables:
#'
#' 1. Conduct PCA on the matrix \eqn{\left(\mathbf{X}_{-j}, \mathbf{Z}_{1:j-1}\right)}.
#' 2. For a fixed \eqn{K}, fit \eqn{\mathbf{x}_j} on \eqn{K} principal components (PCs). The choice of \eqn{K} involves a tradeoff: a larger \eqn{K} makes the knockoff more similar to the original variable, leading to a lower type I error but weaker power of the test.
#' 3. Compute the residual vector \eqn{\varepsilon_n = \mathbf{x}_j - \hat{\mathbf{x}}_j}.
#' 4. Permute \eqn{\varepsilon_n} randomly and denote the permuted residuals as \eqn{\varepsilon_n^*}.
#' 5. Set \eqn{\mathbf{z}_j = \hat{\mathbf{x}}_j + \varepsilon_n^*} and combine it with the current knockoff matrix \eqn{\mathbf{Z}_{1:j-1}}.
#'
#' @param X A numeric matrix representing the original design matrix.
#' @param n_ko An integer specifying the number of knockoff matrices to generate. Default is 1.
#' @param ncomp The number of principal components to use in the knockoff generation process.
#' @param verbose Logical. Whether to display progress information during knockoff generation. Default is TRUE.
#'
#' @return A list of principal component knockoff matrices. Each matrix corresponds to a generated knockoff.
#'
#' @references
#' Jiang, Tao, Yuanyuan Li, and Alison A. Motsinger-Reif. "Knockoff boosted tree for model-free variable selection." Bioinformatics 37.7 (2021): 976-983.
#'
#' Shen, A. et al. (2019). "False discovery rate control in cancer biomarker selection using knockoffs." Cancers, 11, 744.
#'
#' @family create
#' @import stats
#' @examples
#' X <- matrix(rnorm(100), nrow = 10)
#' Xk <- create.pc(X = X, ncomp = 5)
#'
#' @export
#' @md
create.pc <- function(X, n_ko = 1, ncomp, verbose = FALSE) {
  # Initialize dimensions
  n <- nrow(X)
  p <- ncol(X)
  bound <- min(n, p)
  if(is.null(ncomp)) ncomp=floor(p/2)

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
  # Knockoff matrix initialization
  Xk <- matrix(NA, nrow = n, ncol = p)
  rownames(Xk) <- rownames(X)
  colnames(Xk) <- paste0(X.names, ".tilde")

  if (ncomp >= bound) {
    ncomp <- bound - 1
    warning(paste0("Reset ncomp to ", ncomp, " as it exceeds the matrix dimensions."))
  }

  # Initialize result list to store knockoff matrices
  result <- vector(mode = "list", length = n_ko)

  for (ko_index in 1:n_ko) {
    if (verbose) cat('-- Generating knockoff matrix', ko_index, '\n')

    for (j in 1:p) {
      # Perform PCA on X without the j-th column
      pca <- stats::prcomp(X[, -j], center = TRUE, scale. = TRUE)
      PCs <- pca$x[, 1:ncomp]

      # Fit the j-th variable on the PCs
      fit <- stats::.lm.fit(x = PCs, y = X[, j])
      res <- fit$residuals

      # Generate knockoff by permuting residuals and adding fitted values
      res_permuted <- sample(res)
      Xk[, j] <- PCs %*% fit$coefficients + res_permuted

      # Update the design matrix by appending the new knockoff variable
      X <- cbind(X, Xk[, j])
    }

    result[[ko_index]] <- Xk
  }

  return(result)
}

#' Sparse sequential knockoff generation algorithm
#'
#' This function generates sparse sequential knockoff copies of the input data frame \code{X}. Sparse sequential knockoffs first estimate the adjacency matrix of \code{X} (which identifies the zeros/non-zeros in the precision matrix of \code{X}). Then, a modified sequential knockoff algorithm is applied where each regression includes only covariates that correspond to non-zero elements in the precision matrix of \code{X}. The use of a sparse model reduces the number of covariates per regression, improving efficiency.
#'
#' To enhance speed, least squares regression is used by default, unless the number of covariates exceeds half the number of observations (i.e., when \code{p > n/2}), in which case elastic net regularized regression is applied.
#'
#' @param X A data frame or tibble with numeric and factor columns only. The number of columns, \code{ncol(X)}, must be greater than 2.
#' @param n_ko Integer. The number of knockoff matrices to generate. Default is 1.
#' @param adjacency.matrix Optional. A user-specified adjacency matrix (binary indicator matrix corresponding to the non-zero elements of the precision matrix of \code{X}). If not provided, it is estimated within the function.
#' @param verbose Logical. Whether to display progress information during the knockoff generation. Default is TRUE.
#'
#' @return A list of data frames or tibbles, each being a sparse sequential knockoff copy of \code{X}, with the same type and dimensions as \code{X}.
#'
#' @examples
#' set.seed(1)
#' X <- generate_X(n = 100, p = 6, p_b = 2, cov_type = "cov_equi", rho = 0.5)
#' Xk <- create.sparse_seq(X)
#'
#' @family create
#' @export
create.sparse_seq <- function(X, n_ko = 1, adjacency.matrix = NULL, verbose = TRUE) {

  # Ensure input design matrix is valid
  check_design(X)

  # Calculate the adjacency matrix if not provided
  if (is.null(adjacency.matrix)) {
    adjacency.matrix <- glasso_adjacency_matrix(X)
  }

  diag(adjacency.matrix) <- 0 # Ensure no self-adjacency

  # Create a copy of X with tilde names for knockoffs
  Xk <- X
  names(Xk) <- paste0(names(X), ".tilde")

  # Initialize result list to store knockoff matrices
  result <- vector("list", length = n_ko)

  for (ko_index in 1:n_ko) {
    if (verbose) cat('-- Generating knockoff matrix', ko_index, '\n')

    # Randomly shuffle column indices of X
    shf <- sample(ncol(X))

    # Loop through columns of X in random order
    loop.count <- 1
    for (i in shf) {
      y <- X[[i]]

      # Identify the covariates (non-zero adjacency)
      j <- which(adjacency.matrix[, i] != 0)
      Xp <- X[, j, drop = FALSE]

      # Include knockoffs of previously processed columns
      if (loop.count > 1) {
        Xp <- cbind(Xk[, intersect(shf[1:(loop.count - 1)], j)], Xp)
      }

      # Generate knockoffs based on the number of predictors
      if (ncol(Xp) < length(y) / 2) {
        Xk[[i]] <- simple.recovery(y = y, X = Xp)
      } else {
        Xk[[i]] <- glmnet.recovery(y = y, X = Xp)
      }

      loop.count <- loop.count + 1
    }

    # Store the generated knockoff matrix
    result[[ko_index]] <- Xk
  }

  return(result)
}

#' Create Sparse Gaussian Knockoffs
#'
#' Generates sparse Gaussian knockoff variables using the second-order method with the semi-definite programming (SDP) approach.
#'
#' @param X A numeric matrix of predictors for which knockoffs will be generated. Each column represents a variable, and each row corresponds to an observation.
#' @param n_ko An integer specifying the number of knockoff matrices to generate. Default is 1.
#' @param verbose Logical. Whether to display progress information during the knockoff generation. Default is TRUE.
#'
#' @return A list of matrices containing sparse Gaussian knockoff variables corresponding to the original matrix \code{X}. Each matrix in the list has the same dimensions as \code{X}.
#'
#' @details This function generates knockoff variables using the semi-definite programming (SDP) method for second-order knockoffs with a sparse covariance matrix. Sparsity is controlled through the \code{spcov} package, which allows for shrinkage methods in variable selection tasks.
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' Xk <- create.sparse_Gaussian(X)
#'
#' @family create
#' @importFrom spcov spcov
#' @export
create.sparse_Gaussian <- function(X, n_ko=1,verbose = FALSE){
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

  # Assign default column names if necessary
  if (is.null(X.names)) {
    X.names <- paste0('X', 1:p)
    colnames(X) <- X.names
  }

  # Validate input types and parameters
  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')
  if (!is.numeric(n_ko) || length(n_ko) != 1 || n_ko <= 0 || floor(n_ko) != n_ko) stop("n_ko must be a positive integer.")

  # Calculate mean vector and covariance matrix
  mu <- colMeans(X)
  S <- stats::cov(X)
  p <- ncol(X)

  # Parameters for sparse covariance estimation
  step.size <- 100
  lam <- 0.06

  # Create the sparsity pattern matrix P (off-diagonals set to 1)
  P <- matrix(1, p, p)
  diag(P) <- 0

  # Estimate sparse covariance matrix using spcov package
  mm <- spcov::spcov(Sigma = diag(diag(S)), S = (S + 0.1 * diag(1, p)), lambda = lam * P, step.size = step.size)

  # Initialize a list to store knockoff matrices
  result <- vector(mode = "list", length = n_ko)

  # Generate knockoff matrices
  for (ko_index in 1:n_ko) {
    if (verbose) cat('-- Generating knockoff matrix', ko_index, '\n')
    result[[ko_index]] <- create.gaussian(X = X, mu = mu, Sigma = mm$Sigma, method = "sdp")
  }
  return(result)
}

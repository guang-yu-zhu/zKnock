#' Create Shrink Gaussian Knockoffs
#'
#' This function generates Gaussian knockoff variables using the second-order method with the semi-definite programming (SDP) approach. It is particularly useful for variable selection tasks involving shrinkage methods.
#'
#' @param X A numeric matrix or data frame of predictors for which knockoffs will be generated. Each column represents a variable, and each row corresponds to an observation.
#' @param n_ko An integer specifying the number of knockoff matrices to generate. Default is 1.
#' @param verbose Logical. If TRUE, displays progress information during knockoff generation. Default is TRUE.
#'
#' @return A list of matrices containing Gaussian knockoff variables corresponding to the original matrix \code{X}. Each matrix in the list has the same dimensions as \code{X}.
#'
#' @details This function uses the semi-definite programming (SDP) method to generate second-order Gaussian knockoffs. It is useful in statistical settings, such as controlling the false discovery rate in variable selection tasks that involve shrinkage methods.
#'
#' @examples
#' # Generate a random matrix of predictors
#' X <- matrix(rnorm(100), nrow = 10, ncol = 10)
#'
#' # Create shrink Gaussian knockoffs
#' Xk <- create.shrink_Gaussian(X)
#'
#' @family create
#' @export
create.shrink_Gaussian <- function(X, n_ko = 1, verbose = FALSE) {

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

  # Initialize a list to store the knockoff matrices
  result <- vector(mode = "list", length = n_ko)

  # Generate knockoff matrices
  for (ko_index in 1:n_ko) {
    if (verbose) cat('-- Generating knockoff matrix', ko_index, '\n')
    result[[ko_index]] <- create.second_order(X = X, method = "sdp")
  }

  return(result)
}

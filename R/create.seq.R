#' Sequential knockoffs for continuous and categorical variables
#'
#' This function generates sequential knockoff copies of the input data frame \code{X}. Sequential knockoffs simulate new variables for each column of \code{X} using the \code{seq_simulator}, which by default is \code{glmnet.recovery}. This default method fits elastic-net models to sequentially generate knockoffs.
#'
#' @param X A data frame or tibble with numeric and factor columns only. The number of columns, \code{ncol(X)}, must be greater than 2.
#' @param n_ko Integer. The number of knockoff matrices to generate. Default is 1.
#' @param seq_simulator Function that simulates sequential knockoffs. Default is \code{glmnet.recovery}, which simulates responses from an estimated elastic-net model.
#' @param verbose Logical. Whether to display progress information during the knockoff generation. Default is FALSE.
#' @param ... Additional parameters passed to \code{seq_simulator}. For the default elastic-net method, these are passed to \code{cv.glmnet}.
#'
#' @details \code{create.seq} performs sequential knockoff simulations using elastic-net regression as the default method. It loops over the columns of \code{X} and generates knockoffs for each column, using the other columns as predictors.
#'
#' @return A list of data frames or tibbles, each containing the sequential knockoff copies of \code{X}, with the same type and dimensions as the original \code{X}.
#'
#' @examples
#' set.seed(1)
#' X <- generate_X(n = 100, p = 6, p_b = 2, cov_type = "cov_equi", rho = 0.5)
#' Xk <- create.seq(X)
#'
#' @family create
#' @export
create.seq <- function(X, n_ko = 1, seq_simulator = glmnet.recovery, verbose = FALSE, ...) {

  # Ensure X is a data frame
  if (is.matrix(X)) X <- as.data.frame(X)

  # Validate input design
  check_design(X)

  # Create knockoff names for the columns of X
  Xk <- X
  names(Xk) <- paste0(names(Xk), ".tilde")

  # Initialize result list to store knockoff matrices
  result <- vector("list", length = n_ko)

  for (ko_index in 1:n_ko) {
    if (verbose) cat('-- Generating knockoff matrix', ko_index, '\n')

    # Randomly shuffle column indices of X
    shf <- sample(ncol(X))

    # Loop through columns of X in random order
    loop.count <- 1
    for (i in shf) {
      y <- X[[i]]  # i-th column serves as response
      Xp <- X[, -i]  # All columns except i serve as predictors

      # Include previously generated knockoffs for earlier columns
      if (loop.count > 1) Xp <- cbind(Xk[, shf[1:(loop.count - 1)]], Xp)

      # Generate knockoffs using the sequential simulator
      Xk[[i]] <- seq_simulator(y = y, X = Xp, ...)

      loop.count <- loop.count + 1
    }

    # Store generated knockoff matrix in the result list
    result[[ko_index]] <- Xk
  }

  return(result)
}

#' The Knockoff Filter
#'
#' This function runs the Knockoff procedure, selecting variables relevant for predicting
#' the outcome of interest using knockoffs and test statistics.
#'
#' @param X A numeric n-by-p matrix or data frame of predictors.
#' @param y A response vector of length n.
#' @param Xk A list of Knockoff copys.
#' @param statistic A function to compute test statistics (default: `stat.glmnet_coefdiff`).
#' @param aggregate Function to aggregate results from multiple knockoffs (default: `agg_Freq`).
#' @param fdr Target false discovery rate (default: 0.1).
#' @param offset Offset for threshold computation (0 or 1; default: 1).
#' @param verbose Logical; if TRUE, prints progress messages during knockoff generation and statistic calculation (default: FALSE).
#' @param ... Additional arguments passed to the `statistic` function.
#'
#' @return An object of class `knockoff.filter`, containing:
#' - Single Knockoff Case:
#'    - **call**: The matched call of the function.
#'    - **W**: The test statistics for the original variables.
#'    - **threshold**: The computed selection threshold.
#'    - **shat**: The indices of variables selected based on the threshold.
#' - Multiple Knockoff Case:
#'    - **call**: The matched call of the function.
#'    - **shat**: The aggregated indices of selected variables.
#'    - **Ws**: The matrix of test statistics for multiple knockoff copies.
#'    - **thresholds**: A vector of thresholds for each knockoff.
#'    - **shat_list**: A list where each element contains the indices of selected variables for a corresponding knockoff copy.
#'    - **shat_mat**: A binary matrix where each row indicates the selected variables for a specific knockoff copy (1 for selected, 0 for not selected).
#'
#' @references
#' Candes, E., Fan, Y., Janson, L., & Lv, J. (2018). Panning for gold: ‘model-X’ knockoffs for high dimensional controlled variable selection. Journal of the Royal Statistical Society Series B: Statistical Methodology, 80(3), 551-577.
#'
#' @examples
#' # Linear Regression
#' set.seed(2024)
#' n=80; p=100; k=10; Ac = 1:k; Ic = (k+1):p
#' X = generate_X(n=n,p=p)
#' y <- generate_y(X, p_nn=k, a=3)
#' Xk = create.shrink_Gaussian(X = X, n_ko = 10)
#' res1 = knockoff.filter(X, y, Xk, statistic = stat.glmnet_coefdiff,
#'                        offset = 1, fdr = 0.1)
#' res1
#' perf_eval(res1$shat,Ac,Ic)
#'
#' # Logistic Regression
#' lp <- generate_lp(X, p_nn=k, a=3)
#' pis <- plogis(lp)
#' Y <- factor(rbinom(n, 1, pis))
#' res2 = knockoff.filter(X, Y, Xk, statistic = stat.glmnet_coefdiff,
#'                        family = 'binomial', offset = 0, fdr = 0.2)
#' res2
#' perf_eval(res2$shat,Ac,Ic)
#'
#' @export
#' @md
knockoff.filter <- function(X, y, Xk = NULL,
                            statistic = stat.glmnet_coefdiff,
                            aggregate = agg_Freq,
                            fdr = 0.10,
                            offset = 1,
                            verbose = FALSE,
                            ...) {

  # Validate input types
  if (is.data.frame(X)) {
    X.names <- names(X)
    X <- as.matrix(X)
  } else if (is.matrix(X)) {
    X.names <- colnames(X)
  } else {
    stop('Input X must be a numeric matrix or data frame')
  }

  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')
  if ((!is.numeric(y) && !is.factor(y))) stop('Input y must be either numeric or factor type')
  if (is.numeric(y)) y <- as.vector(y)
  if (!offset %in% c(0, 1)) stop('Input offset must be either 0 or 1')
  if (!is.function(statistic)) stop('Input statistic must be a function')

  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(y) == n)

  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')
  # Assign default column names if necessary
  if (is.null(X.names)) {
    X.names <- paste0('X', 1:p)
    colnames(X) <- X.names
  }

  # Create knockoff variables if Xk is not provided
  if (is.null(Xk)) {
    stop('Input Xk must be provided.')
  }

  # Compute statistics for each knockoff
  Ws <- vector(mode = "list", length = length(Xk))
  for (i in seq_along(Xk)) {
    if (verbose) cat('-- Calculating knockoff statistics for copy', i, '.\n')
    Ws[[i]] <- statistic(X, Xk[[i]], y, ...)
  }

  # Convert list of statistics to matrix
  Ws_mat <- do.call(rbind, Ws)

  # Single knockoff case
  if (length(Ws) == 1) {
    W <- Ws_mat
    threshold <- knockoff.threshold(W, fdr, offset)
    shat <- which(W >= threshold)

    result <- list(
      call = match.call(),
      W = W,
      threshold = threshold,
      shat = shat
    )
  } else {  # Multiple knockoff case
    # Aggregate results from multiple knockoffs
    shat <- aggregate(Ws_mat, fdr = fdr, offset = offset)

    # Run separate filtering for each knockoff
    thresholds <- apply(Ws_mat, 1, knockoff.threshold, fdr = fdr, offset = offset)
    shat_mat <- sweep(Ws_mat, 2, thresholds, FUN = `>=`) * 1
    shat_list <- apply(shat_mat, 1, function(vec) which(vec == 1))

    result <- list(
      call = match.call(),
      shat = shat,
      Ws = Ws_mat,
      thresholds = thresholds,
      shat_list = shat_list,
      shat_mat = shat_mat
    )
  }

  class(result) <- 'knockoff.filter'
  return(result)
}

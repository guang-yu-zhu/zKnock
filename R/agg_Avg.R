#' Aggregated Knockoff with Average Test Statistics
#'
#' Performs Aggregated Knockoff analysis by averaging the test statistics across multiple knockoff filters.
#'
#' @param Ws_mat A matrix where each row represents the test statistics from a different knockoff filter and each column corresponds to a variable.
#' @param fdr A numeric value specifying the target false discovery rate (FDR) level. Default is \eqn{0.05}.
#' @param offset An integer (0 or 1) specifying the offset in the threshold calculation. Default is \eqn{0} (liberal control).
#'
#' @return A vector `shat` containing the indices of selected variables after averaging the test statistics and applying the knockoff filter.
#' @family aggregate
#' @examples
#' # Linear Regression
#' set.seed(2024)
#' p = 100; n = 80; k = 10; scale = 3
#' Ac = 1:k
#' rho = 0.3; SigmaX <- toeplitz(rho^(0:(p-1)))
#' SigmaXhalf = chol(SigmaX)
#' beta = matrix(0, p, 1)
#' beta[Ac] = sample(c(-1, 1) * scale, k, replace = TRUE)
#' X = matrix(rnorm(n * p), n) %*% SigmaXhalf
#' y = X %*% beta + rnorm(n)
#' Xk = create.knockoff(X = X, type = 'shrink', n_ko = 10)
#' res1 = knockoff.filter(X, y, Xk, statistic = stat.glmnet_coefdiff,
#'                        aggregate = agg_Freq,
#'                        offset = 1, fdr = 0.1)
#' agg_Avg(res1$Ws)
#'
#' @export
#' @md
agg_Avg <- function(Ws_mat, fdr = 0.05, offset = 0) {

  # Check that input is a matrix
  if (!is.matrix(Ws_mat)) stop('Input Ws_mat must be a matrix')

  # Calculate the average test statistic across knockoff copies
  W = colMeans(Ws_mat)

  # Compute the knockoff threshold
  thre <- knockoff.threshold(W, fdr, offset)

  # Select variables whose test statistics exceed the threshold
  shat <- which(W >= thre)

  return(shat)
}

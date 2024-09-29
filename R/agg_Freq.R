#' Aggregated Knockoff Using Selection Frequency
#'
#' Performs Aggregated Knockoff analysis by selecting variables based on their selection frequency across multiple knockoff filters.
#'
#' @param Ws_mat A matrix where each row represents the test statistics from a different knockoff filter, and each column corresponds to a variable.
#' @param fdr A numeric value specifying the target false discovery rate (FDR) level. Default is \eqn{0.05}.
#' @param offset An integer (0 or 1) specifying the offset in the threshold calculation. Default is \eqn{0} (liberal control).
#' @param tau A numeric value indicating the selection frequency threshold for variable inclusion. Default is \eqn{0.5}.
#'
#' @return A vector `shat` containing the indices of selected variables after applying selection frequency and the knockoff filter.
#' @family aggregate
#' @references
#' Ren, Z., Wei, Y., & Candès, E. (2023). Derandomizing knockoffs. Journal of the American Statistical Association, 118(542), 948-958.
#'
#' @examples
#' set.seed(2024)
#' p = 100; n = 80
#' X = generate_X(n=80,p=100)
#' y <- generate_y(X, p_nn=10, a=3)
#' Xk = create.shrink_Gaussian(X = X, n_ko = 10)
#' res1 = knockoff.filter(X, y, Xk, statistic = stat.glmnet_coefdiff,
#'                        aggregate = agg_BH,
#'                        offset = 1, fdr = 0.1)
#' res1
#' agg_Freq(res1$Ws)
#'
#' @export
#' @md
agg_Freq <- function(Ws_mat, fdr = 0.05, offset = 0, tau = 0.5) {

  # Check that input is a matrix
  if (!is.matrix(Ws_mat)) stop('Input Ws_mat must be a matrix')

  # Compute knockoff thresholds for each row (knockoff filter)
  thresholds <- apply(Ws_mat, 1, knockoff.threshold, fdr = fdr, offset = offset)

  # Create a binary matrix indicating selections based on thresholds
  shat_mat <- sweep(Ws_mat, 2, thresholds, FUN = `>=`) * 1

  # Select variables whose selection frequency exceeds tau
  shat <- which(colMeans(shat_mat) > tau)

  return(shat)
}

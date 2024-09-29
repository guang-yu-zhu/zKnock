#' Aggregated Knockoffs with AKO (Aggregation of Multiple Knockoffs)
#'
#' Performs AKO (Aggregation of Multiple Knockoffs)
#' @details
#' 1. **Calculate intermediate p-value** \eqn{\pi_j^{(b)}}, for all \eqn{j \in [p]} and \eqn{b \in [B]}:
#'    \deqn{
#'    \pi_j = \begin{cases}
#'    \frac{1 + \#\left\{k: W_k \leq -W_j \right\}}{p}, & \text{if } W_j > 0 \\
#'    1, & \text{if } W_j \leq 0
#'    \end{cases}}
#' 2. **Aggregate using the quantile aggregation procedure** (Meinshausen et al. 2009):
#'    \deqn{
#'    \bar{\pi}_j = \min \left\{1, \frac{q_\gamma\left(\left\{\pi_j^{(b)}: b \in [B]\right\}\right)}{\gamma}\right\}
#'    }
#' 3. **Control FDR using Benjamini-Hochberg step-up procedure** (BH, Benjamini & Hochberg 1995):
#'    - Order p-values: \eqn{\bar{\pi}_{(1)} \leq \bar{\pi}_{(2)} \leq \ldots \leq \bar{\pi}_{(p)}}.
#'    - Find: \eqn{\widehat{k}_{BH} = \max \left\{k: \bar{\pi}_{(k)} \leq \frac{k \alpha}{p}\right\}}.
#'    - Select: \eqn{\widehat{\mathcal{S}} = \left\{j \in [p]: \bar{\pi}_{(j)} \leqslant \bar{\pi}_{\left(\widehat{k}_{BH}\right)}\right\}}.
#' @param Ws_mat A matrix of test statistics from multiple knockoff filters, where each row represents one set of test statistics and each column represents a variable.
#' @param fdr A numeric value of the target false discovery rate (FDR) level. Default is \eqn{0.1}.
#' @param offset An integer (0 or 1) specifying the offset in the empirical p-value calculation. Default is \eqn{0}.
#' @param gamma A numeric value for quantile aggregation in the multiple knockoff p-value aggregation. Default is \eqn{0.3}.
#'
#' @return A vector `shat` containing the indices of selected variables after aggregating knockoff results.
#'
#' @references
#'  - Nguyen TB, Chevalier JA, Thirion B, Arlot S. Aggregation of multiple knockoffs. In: International Conference on Machine Learning. PMLR; 2020. p. 7283–93.
#'  - Tian P, Hu Y, Liu Z et al. Grace-AKO: a novel and stable knockoff filter for variable selection incorporating gene network structures. BMC Bioinformatics#' @keywords internal
#' @family aggregate
#'
#' @examples
#' set.seed(2024)
#' p = 100; n = 80
#' X = generate_X(n=80,p=100)
#' y <- generate_y(X, p_nn=10, a=3)
#' Xk = create.shrink_Gaussian(X = X, n_ko = 10)
#' res1 = knockoff.filter(X, y, Xk, statistic = stat.glmnet_coefdiff,
#'                        offset = 1, fdr = 0.1)
#' res1
#' agg_BH(res1$Ws)
#'
#' @export
#' @md
agg_BH <- function(Ws_mat, fdr = 0.1, offset = 0, gamma = 0.3) {
  # Check input
  if (!is.matrix(Ws_mat)) stop('Input Ws_mat must be a matrix')

  p = ncol(Ws_mat)
  n_ko = nrow(Ws_mat)

  # Compute empirical p-values for each knockoff
  pvals = apply(Ws_mat, MARGIN = 1, empirical_pval, offset = offset)

  # Aggregate p-values using quantile aggregation
  aggregated_pval = apply(pvals, 1, quantile_aggregation, gamma = gamma)

  # Compute the FDR control threshold
  threshold = bhq_threshold(aggregated_pval, fdr = fdr)

  # Select variables whose aggregated p-values are below the threshold
  shat <- which(aggregated_pval <= threshold)

  return(shat)
}

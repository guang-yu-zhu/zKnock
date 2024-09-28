#' Compute Empirical P-values from Test Statistics
#'
#' Computes empirical p-values from test statistics. Function Adapted from P Tian et al. (2022).
#'
#' @param test_score A numeric vector of test statistics.
#' @param offset A numeric value to adjust p-values. Options are `0` or `1`. Default is `1`.
#'
#' @return A numeric vector of empirical p-values.
#'
#' @references
#' - Nguyen TB, Chevalier JA, Thirion B, Arlot S. Aggregation of multiple knockoffs. In: International Conference on Machine Learning. PMLR; 2020. p. 7283–93.
#' - Tian P, Hu Y, Liu Z et al. Grace-AKO: a novel and stable knockoff filter for variable selection incorporating gene network structures. BMC Bioinformatics 2022;23:478.
#'
#' @keywords internal
#' @md
empirical_pval = function(test_score, offset = 1){
  pvals = c()
  n_features = length(test_score)
  if (offset !=0 && offset!=1){
    return("'offset' must be either 0 or 1")
  }
  else{
    test_score_inv = -test_score
    for (i in 1:n_features){
      if (test_score[i] <= 0){
        pvals = c(pvals, 1)
      }
      else{
        pvals = c(pvals,(offset+sum(test_score_inv[i] >= test_score))/n_features)
      }
    }
  }
  return (pvals)
}

#' Quantile Aggregation
#'
#' Aggregates multiple p-values using quantile aggregation based on Meinshausen et al. (2009). Function Adapted from P Tian et al. (2022).
#'
#' @param pvals A numeric vector of p-values.
#' @param gamma A numeric value for the quantile to use in aggregation. Default is `0.3`.
#'
#' @return A numeric value representing the aggregated p-value.
#' @references
#' - Nguyen TB, Chevalier JA, Thirion B, Arlot S. Aggregation of multiple knockoffs. In: International Conference on Machine Learning. PMLR; 2020. p. 7283–93.
#' - Tian P, Hu Y, Liu Z et al. Grace-AKO: a novel and stable knockoff filter for variable selection incorporating gene network structures. BMC Bioinformatics#' @keywords internal
#' @keywords internal
#' @md
quantile_aggregation = function(pvals, gamma=0.3){
  converted_score = (1 / gamma) *  quantile(pvals, gamma)
  return (min(1, converted_score))
}

#' BHQ Threshold for FDR Control
#'
#' Computes the Benjamini-Hochberg (BH) threshold for FDR control. Function Adapted from P Tian et al. (2022).
#'
#' @param pvals A numeric vector of p-values.
#' @param fdr A numeric value specifying the false discovery rate level. Default is `0.1`.
#'
#' @return A numeric value representing the threshold for FDR control.
#'
#' @references
#' - Nguyen TB, Chevalier JA, Thirion B, Arlot S. Aggregation of multiple knockoffs. In: International Conference on Machine Learning. PMLR; 2020. p. 7283–93.
#' - Tian P, Hu Y, Liu Z et al. Grace-AKO: a novel and stable knockoff filter for variable selection incorporating gene network structures. BMC Bioinformatics#' @keywords internal
#' @keywords internal
#' @md
bhq_threshold = function(pvals, fdr=0.1){
  n_features = length(pvals)
  pvals_sorted = sort(pvals)
  selected_index = 2 * n_features
  for (i in seq(n_features, 1, -1)){
    if (pvals_sorted[i] <= (fdr * i / n_features)){
      selected_index = i
      break
    }
  }
  if (selected_index <= n_features){
    return (pvals_sorted[selected_index])
  }
  else{
    return ('-1.0')
  }
}

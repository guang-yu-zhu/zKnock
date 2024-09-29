#' Evaluate False Discovery Proportion (FDP) and True Positive Proportion (TPP)
#'
#' This function calculates the false discovery rate and true positive rate
#' based on selected variables, known null variables (negatives), and known
#' non-null variables (positives).
#'
#' @param selected A vector of indices of selected variables.
#' @param positives A vector of indices of known non-null variables.
#' @param negatives A vector of indices of known null variables.
#'
#' @return A numeric vector of length 2 containing:
#' - `tpr`: The true positive rate.
#' - `fdr`: The false discovery rate.
#' @export
#'
#' @examples
#' perf_eval(selected = c(1, 2, 3, 5), positives = 1:4, negatives = 5:10)
perf_eval <- function(selected, positives, negatives) {
  # Calculate False Discovery Rate (FDR)
  fdr <- if (length(selected) > 0) {
    sum(selected %in% negatives) / length(selected)
  } else {
    0
  }

  # Calculate True Positive Rate (TPR)
  tpr <- sum(positives %in% selected) / length(positives)

  return(c(tpr, fdr))
}

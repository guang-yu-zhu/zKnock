#' Create PC Knockoffs
#' A sequential algorithm to create non-parametric knockoffs based on principal component regression and residuals permutation.
#'
#' @details
#' For each original variable \eqn{\mathbf{x}_j}, where \eqn{j = 1, \ldots, p}, the following steps are performed to generate knockoff variables:
#'
#' 1. Conduct PCA on the matrix \eqn{\left(\mathbf{X}_{-j}, \mathbf{Z}_{1: j-1}\right)}.
#' 2. For a fixed \eqn{K}, fit \eqn{\mathbf{x}_j} on \eqn{K} PCs. There is a tradeoff in that the larger the \eqn{K}, the more akin the knockoff will be to the original variables. This results in a smaller type 1 error but weaker power of the test.
#' 3. Compute a residual vector \eqn{\varepsilon_n=\left(\mathbf{x}_j-\hat{\mathbf{x}}_j\right)}.
#' 4. Permute \eqn{\varepsilon_n} randomly. Denote the permuted vector as \eqn{\varepsilon_n^*}.
#' 5. Set \eqn{\mathbf{z}_j=\hat{\mathbf{x}}_j+\varepsilon_n^*} and combine it with the current knockoff matrix \eqn{\mathbf{Z}_{1: j-1}}.
#'
#' @param X An input original design matrix.
#' @param pc.num The number of pricial components to be used for generating knockoff matrices.
#'
#' @return A principal component knockoff matrix.
#'
#' @references
#' Jiang, Tao, Yuanyuan Li, and Alison A. Motsinger-Reif. "Knockoff boosted tree for model-free variable selection." Bioinformatics 37.7 (2021): 976-983.
#'
#' Shen,A. et al. (2019) False discovery rate control in cancer biomarker selection using knockoffs. Cancers, 11, 744.
#' @family create
#' @export
#' @import stats
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(100), nrow = 10)
#' Xk <- create.pc.knockoff(X = X, pc.num = 5)
create.pc.knockoff <- function(X, pc.num) {
  p <- ncol(X)
  n <- nrow(X)
  bound <- min(n,p)

  if (pc.num >= bound) {
    pc.num <- (bound - 1)
    warning(paste0("Reset num.comp as ", pc.num, "\n"))
  }

  Z <- matrix(NA, ncol = p, nrow = n)
  for (i in 1:p) {
    pca <- stats::prcomp(X[,-i], center = TRUE, scale. = TRUE)
    PCs <- pca$x[, 1:pc.num]

    fit <- stats::.lm.fit(x = PCs, y = X[,i])
    res <- fit$residuals

    tmp <- sample(res) + PCs%*%fit$coefficients
    Z[,i] <- tmp
    X <- cbind(X, tmp)
  }

  return (Z)
}

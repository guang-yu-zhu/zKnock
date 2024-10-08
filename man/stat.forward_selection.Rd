% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat.forward_selection.R
\name{stat.forward_selection}
\alias{stat.forward_selection}
\title{Importance statistics based on forward selection}
\usage{
stat.forward_selection(X, X_k, y, omp = F)
}
\arguments{
\item{X}{n-by-p matrix of original variables.}

\item{X_k}{n-by-p matrix of knockoff variables.}

\item{y}{numeric vector of length n, containing the response variables.}

\item{omp}{whether to use orthogonal matching pursuit (default: F).}
}
\value{
A vector of statistics \eqn{W} of length p.
}
\description{
Computes the statistic
\deqn{W_j = \max(Z_j, Z_{j+p}) \cdot \mathrm{sgn}(Z_j - Z_{j+p}),}
where \eqn{Z_1,\dots,Z_{2p}} give the reverse order in which the 2p
variables (the originals and the knockoffs) enter the forward selection
model.
See the Details for information about forward selection.
}
\details{
In \emph{forward selection}, the variables are chosen iteratively to maximize
the inner product with the residual from the previous step. The initial
residual is always \code{y}. In standard forward selection
(\code{stat.forward_selection}), the next residual is the remainder after
regressing on the selected variable; when orthogonal matching pursuit
is used, the next residual is the remainder
after regressing on \emph{all} the previously selected variables.
}
\examples{
set.seed(2024)
n=80; p=100; k=10; Ac = 1:k; Ic = (k+1):p
X = generate_X(n=n,p=p)
y <- generate_y(X, p_nn=k, a=3)
Xk = create.shrink_Gaussian(X = X, n_ko = 10)
res1 = knockoff.filter(X, y, Xk, statistic = stat.forward_selection,
                       offset = 1, fdr = 0.1)
res1
perf_eval(res1$shat,Ac,Ic)

}
\seealso{
Other statistics: 
\code{\link{stat.SHAP}()},
\code{\link{stat.glmnet_coefdiff}()},
\code{\link{stat.glmnet_lambdadiff}()},
\code{\link{stat.random_forest}()},
\code{\link{stat.sqrt_lasso}()},
\code{\link{stat.stability_selection}()},
\code{\link{stat.xgboost}()}
}
\concept{statistics}

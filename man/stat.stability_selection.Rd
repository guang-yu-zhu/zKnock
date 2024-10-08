% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat.stability_selection.R
\name{stat.stability_selection}
\alias{stat.stability_selection}
\title{Importance statistics based on stability selection}
\usage{
stat.stability_selection(X, X_k, y, fitfun = stabs::lars.lasso, ...)
}
\arguments{
\item{X}{n-by-p matrix of original variables.}

\item{X_k}{n-by-p matrix of knockoff variables.}

\item{y}{response vector (length n)}

\item{fitfun}{fitfun a function that takes the arguments x, y as above,
and additionally the number of variables to include in each model q.
The function then needs to fit the model and to return a logical vector
that indicates which variable was selected (among the q selected variables).
The name of the function should be prefixed by 'stabs::'.}

\item{...}{additional arguments specific to 'stabs' (see Details).}
}
\value{
A vector of statistics \eqn{W} of length p.
}
\description{
Computes the difference statistic
\deqn{W_j = |Z_j| - |\tilde{Z}_j|}
where \eqn{Z_j} and \eqn{\tilde{Z}_j} are measure the importance
of the jth variable and its knockoff, respectively, based on the
stability of their selection upon subsampling of the data.
}
\details{
This function uses the \code{stabs} package to compute
variable selection stability. The selection stability of the j-th
variable is defined as its probability of being selected upon random
subsampling of the data. The default method for selecting variables
in each subsampled dataset is \code{\link[stabs:fitfuns]{stabs::lars.lasso()}}.

For a complete list of the available additional arguments, see \code{\link[stabs:stabsel]{stabs::stabsel()}}.
}
\examples{
set.seed(2024)
n=80; p=100; k=10; Ac = 1:k; Ic = (k+1):p
X = generate_X(n=n,p=p)
y <- generate_y(X, p_nn=k, a=3)
Xk = create.shrink_Gaussian(X = X, n_ko = 10)
res1 = knockoff.filter(X, y, Xk, statistic = stat.stability_selection,
                       offset = 1, fdr = 0.1)
res1
perf_eval(res1$shat,Ac,Ic)


}
\seealso{
Other statistics: 
\code{\link{stat.SHAP}()},
\code{\link{stat.forward_selection}()},
\code{\link{stat.glmnet_coefdiff}()},
\code{\link{stat.glmnet_lambdadiff}()},
\code{\link{stat.random_forest}()},
\code{\link{stat.sqrt_lasso}()},
\code{\link{stat.xgboost}()}
}
\concept{statistics}

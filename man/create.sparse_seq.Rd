% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.sparse_seq.R
\name{create.sparse_seq}
\alias{create.sparse_seq}
\title{Sparse sequential knockoff generation algorithm}
\usage{
create.sparse_seq(X, n_ko = 1, adjacency.matrix = NULL, verbose = TRUE)
}
\arguments{
\item{X}{A data frame or tibble with numeric and factor columns only. The number of columns, \code{ncol(X)}, must be greater than 2.}

\item{n_ko}{Integer. The number of knockoff matrices to generate. Default is 1.}

\item{adjacency.matrix}{Optional. A user-specified adjacency matrix (binary indicator matrix corresponding to the non-zero elements of the precision matrix of \code{X}). If not provided, it is estimated within the function.}

\item{verbose}{Logical. Whether to display progress information during the knockoff generation. Default is TRUE.}
}
\value{
A list of data frames or tibbles, each being a sparse sequential knockoff copy of \code{X}, with the same type and dimensions as \code{X}.
}
\description{
This function generates sparse sequential knockoff copies of the input data frame \code{X}. Sparse sequential knockoffs first estimate the adjacency matrix of \code{X} (which identifies the zeros/non-zeros in the precision matrix of \code{X}). Then, a modified sequential knockoff algorithm is applied where each regression includes only covariates that correspond to non-zero elements in the precision matrix of \code{X}. The use of a sparse model reduces the number of covariates per regression, improving efficiency.
}
\details{
To enhance speed, least squares regression is used by default, unless the number of covariates exceeds half the number of observations (i.e., when \code{p > n/2}), in which case elastic net regularized regression is applied.
}
\examples{
set.seed(1)
X <- generate_X(n = 100, p = 6, p_b = 2, cov_type = "cov_equi", rho = 0.5)
Xk <- create.sparse_seq(X)

}
\seealso{
Other create: 
\code{\link{create.fixed}()},
\code{\link{create.gaussian}()},
\code{\link{create.pc}()},
\code{\link{create.pls}()},
\code{\link{create.second_order}()},
\code{\link{create.seq}()},
\code{\link{create.shrink_Gaussian}()},
\code{\link{create.sparse_Gaussian}()},
\code{\link{create.zpls}()}
}
\concept{create}

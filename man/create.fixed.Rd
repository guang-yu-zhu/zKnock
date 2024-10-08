% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.fixed.R
\name{create.fixed}
\alias{create.fixed}
\title{Fixed-X knockoffs}
\usage{
create.fixed(
  X,
  method = c("sdp", "equi"),
  sigma = NULL,
  y = NULL,
  randomize = F
)
}
\arguments{
\item{X}{normalized n-by-p matrix of original variables.(\eqn{n \geq p}).}

\item{method}{either "equi" or "sdp" (default: "sdp").
This determines the method that will be used to minimize the correlation between the original variables and the knockoffs.}

\item{sigma}{the noise level, used to augment the data with extra rows if necessary (default: NULL).}

\item{y}{vector of length n, containing the observed responses.
This is needed to estimate the noise level if the parameter \code{sigma} is not provided,
in case \eqn{p \leq n < 2p} (default: NULL).}

\item{randomize}{whether the knockoffs are constructed deterministically or randomized (default: F).}
}
\value{
An object of class "knockoff.variables". This is a list
containing at least the following components:
\item{X}{n-by-p matrix of original variables (possibly augmented or transformed).}
\item{Xk}{n-by-p matrix of knockoff variables.}
\item{y}{vector of observed responses (possibly augmented). }
}
\description{
Creates fixed-X knockoff variables.
}
\details{
Fixed-X knockoffs assume a homoscedastic linear regression model for \eqn{Y|X}. Moreover, they only guarantee
FDR control when used in combination with statistics satisfying the "sufficiency" property.
In particular, the default statistics based on the cross-validated lasso does not satisfy this
property and should not be used with fixed-X knockoffs.
}
\references{
Barber and Candes,
Controlling the false discovery rate via knockoffs.
Ann. Statist. 43 (2015), no. 5, 2055--2085.
}
\seealso{
Other create: 
\code{\link{create.gaussian}()},
\code{\link{create.pc}()},
\code{\link{create.pls}()},
\code{\link{create.second_order}()},
\code{\link{create.seq}()},
\code{\link{create.shrink_Gaussian}()},
\code{\link{create.sparse_Gaussian}()},
\code{\link{create.sparse_seq}()},
\code{\link{create.zpls}()}
}
\concept{create}
\keyword{internal}

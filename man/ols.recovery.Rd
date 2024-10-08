% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.recovery.R
\name{ols.recovery}
\alias{ols.recovery}
\title{Calculate \eqn{\hat{X}} by fitting OLS regression on its neighbors}
\usage{
ols.recovery(Y, X)
}
\arguments{
\item{Y}{A numeric vector representing the response variable.}

\item{X}{A numeric matrix or data frame representing the predictor variables. The function removes any duplicated columns from \code{X}.}
}
\value{
A numeric vector \eqn{\hat{Y}} containing the fitted values from the OLS regression.
}
\description{
This function performs ordinary least squares (OLS) regression, estimating the coefficients that best fit the response variable \eqn{Y} given the predictor matrix \eqn{X}. The fitted values \eqn{\hat{Y}} are then returned.
}
\details{
The function first converts the input matrix \code{X} to a data frame to remove any duplicate columns, ensuring that the matrix used for regression is unique in terms of predictors. It then applies OLS regression by solving the normal equations to estimate the coefficients \eqn{\beta}. The fitted values \eqn{\hat{Y}} are computed as \eqn{X \beta}.
}
\seealso{
Other recovery: 
\code{\link{glmnet.recovery}()},
\code{\link{pls.recovery}()},
\code{\link{simple.recovery}()},
\code{\link{spls.recovery}()}
}
\concept{recovery}

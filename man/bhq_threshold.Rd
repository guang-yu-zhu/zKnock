% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.agg.R
\name{bhq_threshold}
\alias{bhq_threshold}
\title{BHQ Threshold for FDR Control}
\usage{
bhq_threshold(pvals, fdr = 0.1)
}
\arguments{
\item{pvals}{A numeric vector of p-values.}

\item{fdr}{A numeric value specifying the false discovery rate level. Default is \code{0.1}.}
}
\value{
A numeric value representing the threshold for FDR control.
}
\description{
Computes the Benjamini-Hochberg (BH) threshold for FDR control. Function Adapted from P Tian et al. (2022).
}
\references{
\itemize{
\item Nguyen TB, Chevalier JA, Thirion B, Arlot S. Aggregation of multiple knockoffs. In: International Conference on Machine Learning. PMLR; 2020. p. 7283–93.
\item Tian P, Hu Y, Liu Z et al. Grace-AKO: a novel and stable knockoff filter for variable selection incorporating gene network structures. BMC Bioinformatics#' @keywords internal
}
}
\keyword{internal}

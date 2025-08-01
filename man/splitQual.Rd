% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NADA_All.R
\name{splitQual}
\alias{splitQual}
\title{Internal: Split Qualified Numeric Vector Into Observations and Censoring}
\usage{
splitQual(v, qual.symbol = "<")
}
\arguments{
\item{v}{Character vector with qualifying symbols like "<0.5".}

\item{qual.symbol}{Character symbol used to indicate censoring (default "<").}
}
\value{
A list with elements \code{obs} (numeric) and \code{cen} (logical).
}
\description{
Internal: Split Qualified Numeric Vector Into Observations and Censoring
}
\keyword{internal}

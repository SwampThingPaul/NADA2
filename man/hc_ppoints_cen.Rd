% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NADA_ros.R
\name{hc_ppoints_cen}
\alias{hc_ppoints_cen}
\title{Plotting Positions for Censored Observations (Cohn Method)}
\usage{
hc_ppoints_cen(obs, censored, cn = NULL, na.action = getOption("na.action"))
}
\arguments{
\item{obs}{A numeric vector of observed values.}

\item{censored}{A logical vector indicating which observations are censored (\code{TRUE} for censored, \code{FALSE} otherwise).}

\item{cn}{An optional list containing Cohn grouping information (usually from \code{\link[=cohn]{cohn()}}); if missing, it will be computed internally.}

\item{na.action}{A function to handle missing values (default is \code{getOption("na.action")}).}
}
\value{
A numeric vector of plotting positions corresponding to the censored observations.
}
\description{
Computes plotting positions for censored observations using the
Cohn method grouping. This is primarily for hydrologic data with left-censoring.
}
\keyword{internal}

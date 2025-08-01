% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NADA_summary.R
\name{censtats}
\alias{censtats}
\title{Summary Statistics for Censored Data using ROS, MLE, and K-M}
\usage{
censtats(...)
}
\arguments{
\item{...}{to pass argument}
}
\value{
A data frame with the summary statistics.
}
\description{
A convenience function that produces a comparative table of
summary statistics obtained using the \code{ros()}, \code{cenmle()},
and \code{cenfit()} routines. These methods are Regression on
Order Statistics (ROS), Maximum Likelihood Estimation (MLE), and
Kaplan-Meier (K-M).
}
\details{
If the data do not fulfill the criteria for the application of
any method, no summary statistics will be produced.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfit.R
\name{cfit}
\alias{cfit}
\title{Compute an ECDF and Distribution Parameters for Censored Data}
\usage{
cfit(
  y1,
  y2,
  conf = 0.95,
  qtls = c(0.1, 0.25, 0.5, 0.75, 0.9),
  plot.pos = 0.375,
  q.type = 7,
  Cdf = TRUE,
  printstat = TRUE,
  Ylab = NULL
)
}
\arguments{
\item{y1}{Concentrations plus detection limits for indicator formatted data.}

\item{y2}{Censoring indicators (logical. 1 or \code{TRUE} = censored, 0 or FALSE = detected) for indicator formatted data.}

\item{conf}{The confidence coefficient for confidence intervals around the Kaplan-Meier mean and median. Default = 0.95.}

\item{qtls}{Probabilities for the quantiles to be estimated.  Defaults are (0.10, 0.25, 0.50, 0.75, 0.90).  You may add and/or substitute probabilities -- all must be between and not including 0 to 1.}

\item{plot.pos}{numeric scalar between 0 and 1 containing the value of the plotting position constant. The default value is \code{plot.pos=0.375}, the Blom plotting position}

\item{q.type}{an integer between 1 and 9 selecting one of the nine quantile algorithms detailed below to be used. See \code{stats::quantile} for more detail, default is set to 7.}

\item{Cdf}{Logical \code{TRUE}/\code{FALSE} indicator of whether to plot the empirical cumulative distribution function (cdf).}

\item{printstat}{Logical \code{TRUE}/\code{FALSE} option of whether to print the resulting statistics in the console window, or not.  Default is \code{TRUE.}}

\item{Ylab}{Optional input text in quotes to be used as the variable name on the ecdf plot.  The default is the name of the \code{y1} input variable.}
}
\value{
If \code{printstat=TRUE}: Based on the provided \code{conf} value, Kaplan-Meier summary statistics (\code{mean},\code{sd},\code{median}), lower and upper confidence intervals around the mean and median value, sample size and percent of censored samples are returned. The specified quantile values are also printed and returned.

If \code{Cdf=TRUE}: The ecdf of censored data is plotted.
}
\description{
Computes the empirical cumulative distribution function (ECDF) for censored data. Estimates parameters of the distribution, including the mean and quantiles.
}
\details{
Quantiles and parameters are estimated using the enparCensored and ecdfPlotCensored functions of the EnvStats package. This avoids a small bias in the mean produced by the NADA package's cenfit function, which uses the reverse Kaplan-Meier procedure, converting left-censored to right-censored data prior to computing the ecdf and mean. See Gillespie et al. for more discussion on the bias.

All printed values will also be output to an object if saved.  Values are character because of the possibility of a \verb{<1}, but if no \code{<} symbol can be converted to numeric using the \code{as.numeric(...)} command.  For data without censoring cfit will also return values.  In that case the values labeled "KM" are not both Kaplan-Meier results and standard arithmetic mean, t-interval CIs on the mean, and quantiles.  See ?quantiles for choosing the q.type; default q.type = 7.
}
\examples{

data(Brumbaugh)

cfit(Brumbaugh$Hg,Brumbaugh$HgCen)

}
\references{
Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.

Gillespie, B.W., et al., 2010.  Estimating Population Distributions When Some Data Are Below a Limit of Detection by Using a Reverse Kaplan-Meier Estimator. Epidemiology 21, 564-570.

Millard, S.P, 2013. EnvStats: An R Package for Environmental Statistics, 2nd ed. Springer Science+Business Media, USA, N.Y.  DOI 10.1007/978-1-4614-8456-1© Springer Science+Business Media New York 2013”

Excerpt From: Steven P. Millard. “EnvStats.” Apple Books.
}
\seealso{
\link[survival:survfit]{survival::survfit}
}

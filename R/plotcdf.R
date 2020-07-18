#' Censored Empirical Cumulative Distribution Function
#'
#' @description Plots cdfs of one or more groups of censored data.  Illustrates the differences between groups for group tests such as those done using cen1way or cenanova.
#' @param y.var The column of data values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param group Optional - grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param logscale Indicator of whether to plot the data in original y units, or on a log scale.  The default is to use the original units, which usually makes more sense for a cdf plot.
#' @param Ylab Optional – input text in quotes to be used as the variable name on the cdf plot.  The default is the name of the `y.var` input variable.
#' @keywords CDF
#' @return Plots an empirical cumulative distribution function , if `group=NULL` then a ECDF with 95% confidence interval is produced. If `group` is identified then ECDFs are produced for each group.
#' @export
#' @importFrom NADA cenfit
#'
#' @examples
#' #' library(NADA) #For example data
#' data(Golden)
#'
#' # with groups
#' with(Golden,plotcdf(Liver,LiverCen,DosageGroup))
#'
#' # all data
#' with(Golden,plotcdf(Liver,LiverCen))

plotcdf <- function(y.var, cen.var, group=NULL, logscale=FALSE, Ylab=varname) {
  varname <- deparse(substitute(y.var))
  log = ""
  if (logscale == TRUE) {log="x"}
  if ( is.null(group) ) {

    cenind <- as.logical(cen.var)
    cenfit.out <- cenfit(y.var, cenind)

    stripaxes <- function(..., axes) plot(...)
    stripaxes(cenfit.out, log=log, xlab=Ylab, lwd = 2, axes=F, main = "Cumulative Distribution Function")
    abline(h=1, lwd=3, col = "white")
    box(bty = "o")
    axis(side = 1)
    axis(side = 2)
  }
  else {
    Factor <- as.factor(group)
    cenind <- as.logical(cen.var)
    cenfit.out <- cenfit(y.var, cenind, Factor)
    grpname <- deparse(substitute(group))
    i <- length(levels(Factor))
    clrs <- c (1:i)

    if (logscale == TRUE)  {
        stripaxes <- function(..., axes) plot(...)
        P1=stripaxes(cenfit.out, log=log, xlab=Ylab, lwd = 2, col=clrs, axes=F, main = "Cumulative Distribution Function")
        abline(h=1, lwd=3, col = "white")
        box(bty = "o")
        axis(side = 1)
        axis(side = 2)
    }
    else {    # logscale is false
        stripaxes <- function(..., axes) plot(...)
        P1=stripaxes(cenfit.out, log=log, xlab=Ylab, lwd = 2, col=clrs, axes=F, main = "Cumulative Distribution Function")
        abline(h=1, lwd=3, col = "white")
        box(bty = "o")
        axis(side = 1)
        axis(side = 2)
    }
    legend("bottomright",levels(Factor), lty=1:i, lwd=2, text.col=clrs, col = clrs, title = grpname)
  }
  }

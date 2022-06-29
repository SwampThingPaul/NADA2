#' Censored Empirical Cumulative Distribution Function
#'
#' @description Plots ecdfs of one or more groups of censored data.  Illustrates the differences between groups for group tests such as those done using cen1way or cenanova.
#' @param x.var The column of data values plus detection limits
#' @param cens.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param xgroup Optional - grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param xlim Limits for the x (data) axis of the ecdf plot.  Default is 0 to the maximum of the y.var variable.  To change, use option xlim = c(0, 50) if 50 is to be the maximum on the plot.
#' @param Ylab Optional – input text in quotes to be used as the variable name on the ecdf plot.  The default is the name of the `y.var` input variable.
#' @keywords CDF ECDF
#' @return Plots an empirical cumulative distribution function. If `group=NULL` the ECDF of the entire dataset is produced. If `group` is identified then ECDFs are plotted for each group seperately and a legend added.
#' @export
#' @importFrom EnvStats ecdfPlotCensored ecdfPlot

#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#' Millard, S.P., 2013. EnvStats: An R Package for Environmental Statistics. Springer-Verlag, New York.
#'
#'@examples
#' data(PbHeron)
#'
#' # with groups
#' with(PbHeron, cen_ecdf(Liver, LiverCen, DosageGroup))
#'
#' # all data
#' with(PbHeron, cen_ecdf(Liver, LiverCen))

cen_ecdf <- function(x.var, cens.var, xgroup=NULL, xlim = c(0, max(y.var)), Ylab=varname) {
  varname <- deparse(substitute(x.var))
  if (is.null(xgroup) == TRUE) { ydat <- na.omit(data.frame(x.var, cens.var))
  y.var <- ydat[,1];  cen.var <- ydat[,2]; group <- xgroup
  }
  else {ydat <- na.omit(data.frame(x.var, cens.var, xgroup))
  y.var <- ydat[,1];  cen.var <- ydat[,2]; group <- ydat[,3]
  }

  if ( is.null(group) ) {
    if (sum(cen.var) != 0 )
    { ecdfPlotCensored(y.var, cen.var, main = "ECDF for Censored Data", xlab = Ylab, ecdf.lwd = 1, type = "s")}
    # no nondetects
    else {ecdfPlot(y.var, main = paste("ECDF for", varname), xlab = Ylab, ecdf.lwd = 1, type = "s")}
  }

  # plot data for multiple groups on same graph
  else {
    Factor <- as.factor(group)
    factorname <- deparse(substitute(xgroup))
    ngp <- length(levels(Factor))
    clrs <- c (1:ngp)
    groupnames <- as.character(levels(Factor))
    if (sum(cen.var[Factor == groupnames[1]]) != 0 )
    { ecdfPlotCensored(y.var[Factor==groupnames[1]], cen.var[Factor==groupnames[1]], main = "ECDF for Censored Data", xlab = Ylab, ecdf.lwd = 2, type = "s", xlim = xlim) }
    # no nondetects
    else {ecdfPlot(y.var[Factor==groupnames[1]], main = "ECDF for Censored Data", xlab = Ylab, ecdf.lwd = 2, type = "s", xlim = xlim) }

    for (i in 2:ngp)    {
      if (sum(cen.var[Factor == groupnames[i]]) != 0)
      { ecdfPlotCensored(y.var[Factor==groupnames[i]], cen.var[Factor==groupnames[i]], add = TRUE, ecdf.col = clrs[i], ecdf.lwd = 2, ecdf.lty = clrs[i], type = "s") }
      # no nondetects
      else {ecdfPlot(y.var[Factor==groupnames[i]], add = TRUE, ecdf.col = clrs[i], ecdf.lwd = 2, ecdf.lty = clrs[i], type = "s") }
    }
    legend("bottomright",levels(Factor), lty=1:ngp, lwd=2, text.col=clrs, col = clrs, title = factorname)
  }
}

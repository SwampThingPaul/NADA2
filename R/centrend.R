#' Trend analysis of residuals of censored data
#'
#' @description Trend analysis of residuals of censored data
#' @param y.var The column of y (response variable) values plus detection limits
#' @param y.cens The column of indicators, where 1 (or `TRUE`) indicates a detectionlimit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.var Column of a covariate (not time).  `y.var` will be smoothed versus `x.var` and residuals taken to subtract out the relationship between `y` and `x`.
#' @param time.var Column of the time variable, either a sequence of days or decimal times, etc.  Will be the scale used for time in ATS trend analysis.
#' @param link Default = `“identity”` which means the original data.
#' @param Smooth Type of smoother used in the GAM.Default is `“cs”`, shrinkage cubic regression splines.See details
#' @keywords trend analysis GAM spline
#' @export
#'
#' @importFrom mgcv gam
#' @importFrom NADA cenxyplot
#' @importFrom cenGAM tobit1
#' @return
#'
#' Prints three plots (Data with GAM Smooth, Residuals from GAM Smooth and ATS trend.)
#'
#' Returns GAM residuals, and ATS results
#'
#' @details
#'
#' Default `link` is identity, other options are available see `cenGAM::tobit1` for more options.
#'
#' Default `Smooth` is `"cs"` for shrinkage cubic regression splines. See `mgcv::smooth.terms` for other smoothing terms.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#'@seealso [mgcv::gam]
#'
#' @examples
#'
#' data(Brumbaugh)
#'
#' Brumbaugh$time=1:nrow(Brumbaugh)
#'
#' with(Brumbaugh,centrend(Hg,HgCen,SedTotHg,time.var=time))


centrend <- function(y.var, y.cens, x.var, time.var, link = "identity", Smooth = "cs") {
  yname <- deparse(substitute(y.var))
  xname <- deparse(substitute(x.var))
  tname <- deparse(substitute(time.var))
  y.txt <- paste(yname, "residuals")
  dl = y.cens * y.var * 1.001
  #  y.cens <- as.factor(y.cens)
  dat.all <- data.frame(y.var, x.var, time.var, dl, y.cens)
  dat.nonas <- na.omit(dat.all)
  colnames(dat.nonas) <- c(yname, xname, tname, "DL", "Cens = 1")

  # cs, ts, ds  top 3, but ds gives bends that may be overfitting.
  gam.y <- gam(dat.nonas[,1] ~ s(dat.nonas[,2], bs = Smooth), family = tobit1(link = link, left.threshold = dat.nonas[,4]))
  dat.out <- data.frame(gam.y$residuals, dat.nonas[,2])
  colnames(dat.out) <- c("GAMresidual", xname)
  # sorting fitted values
  o <- order(dat.nonas[,2] , dat.nonas [,1])
  #plot(dat.nonas[,1] ~ dat.nonas[,2], main = "1. Data and GAM Smooth", ylab = yname, xlab = xname)

  oldpar<- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(3,1))

  x.cens = rep(0, times=length(dat.nonas[,3]) )
  cenxyplot(dat.nonas[,2], as.logical(x.cens), dat.nonas[,1], as.logical(dat.nonas[,5]), main = "1. Data and GAM Smooth", ylab = yname, xlab = xname, pch = 19, cex = 0.7)
  # draw the smooth
  lines (dat.nonas[o,2], gam.y$fitted.values[o], col = 'red')
  plot(gam.y$residuals ~ dat.nonas[,2], main = "2. Residuals from GAM Smooth", xlab = xname, ylab = y.txt)

  ATS(gam.y$residuals, dat.nonas[,5], dat.nonas[,3], x.cens, LOG = FALSE, xlabel = tname, ylabel = y.txt)
  return (invisible(dat.out))
}

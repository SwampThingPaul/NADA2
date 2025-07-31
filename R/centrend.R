#' Trend analysis of censored data with a covariate
#'
#' @description Computes the ATS (Mann-Kendall trend test for censored data) after adjustment of censored data for a covariate.
#' @param y.var The column of y (response variable) values plus detection limits
#' @param y.cens The column of indicators, where 1 (or `TRUE`) indicates a detectionlimit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.var Column of a covariate (not time).  `y.var` will be smoothed versus `x.var` and residuals taken to subtract out the relationship between `y` and `x`.
#' @param time.var Column of the time variable, either a sequence of days or decimal times, etc.  Will be the scale used for time in ATS trend analysis.
#' @param link Default = `“identity”` which means it uses data in the original units. See details.
#' @param Smooth Type of smoother used in the GAM. Default is `“cs”`, shrinkage cubic regression splines. See details for other options.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @param stackplots logical `TRUE`/`FALSE` option to stack three plots that are output onto the same page instead of each on separate page.  Default is 'FALSE', each separate.
#' @param drawplot Logical `TRUE`/`FALSE` option of whether to draw plots or not. Default is `TRUE`
#' @keywords trend analysis GAM spline
#' @export
#'
#' @importFrom mgcv gam
#' @importFrom cenGAM tobit1
#' @return
#'
#' Prints three plots: Y data vs time with GAM Smooth, Residuals from GAM Smooth vs time, and ATS trend line of residuals vs time.
#'
#' Returns GAM residuals and ATS results on trend test of residuals (intercept, slope, Kendall's tau, p-value for trend)
#'
#' @details
#'
#' Default `link` = identity. The y variables are then used in their original units. Other options are available see `cenGAM::tobit1` for more options.
#'
#' Default `Smooth` is `"cs"` for shrinkage cubic regression splines. See `mgcv::smooth.terms` for other types of smoothing algorithms.  '"ts"' is a thin-plate regression spline and is also commonly used.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#'
#'@seealso [mgcv::gam]
#'
#' @examples
#'
#' data(Brumbaugh)
#'
#' Brumbaugh$time=1:nrow(Brumbaugh)
#'
#' with(Brumbaugh,centrend(Hg,HgCen,SedTotHg,time.var=time,drawplot=TRUE))

centrend <- function(y.var, y.cens, x.var, time.var, link = "identity",
                     Smooth = "cs", printstat=TRUE, stackplots = FALSE,
                     drawplot = TRUE) {
  yname <- deparse(substitute(y.var))
  xname <- deparse(substitute(x.var))
  tname <- deparse(substitute(time.var))
  y.txt <- paste(yname, "residuals")

  y.cens <- as.logical(y.cens)
  y.var <- as.numeric(y.var)

  DL <- y.cens * y.var * 1.001
  #  y.cens <- as.factor(y.cens)
  dat.all <- data.frame(y.var, x.var, time.var, DL, y.cens)
  dat.nonas <- na.omit(dat.all)
  colnames(dat.nonas) <- c(yname, xname, tname, "DL", "Cens")

  # cs, ts, ds  top 3, but ds gives bends that may be overfitting.
  gam.y <- gam(dat.nonas[,1] ~ s(dat.nonas[,2], bs = Smooth), family = tobit1(link = link, left.threshold = dat.nonas[,4]))
  # sorting fitted values
  o <- order(dat.nonas[,2] , dat.nonas [,1])
  #plot(dat.nonas[,1] ~ dat.nonas[,2], main = "1. Data and GAM Smooth", ylab = yname, xlab = xname)



  if (stackplots) {
    oldpar<- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(3,1))
  }
  # plot 1.  Y data vs covariate, with smooth
  x.cens = rep(0, times=length(dat.nonas[,3]) )
  if(drawplot == TRUE){
  cenxyplot(dat.nonas[,2], as.logical(x.cens), dat.nonas[,1], as.logical(dat.nonas[,5]), main = "1. Data and GAM Smooth", ylab = yname, xlab = xname, pch = 19, cex = 0.7)
  # draw the smooth
  lines (dat.nonas[o,2], gam.y$fitted.values[o], col = 'red')

  # plot 2. Y residuals from smooth vs. covariate
  plot(gam.y$residuals ~ dat.nonas[,2], main = "2. Residuals from GAM Smooth", xlab = xname, ylab = y.txt)
  abline (h=0, col = "blue")
  }

  if(printstat==TRUE){cat("Trend analysis of", yname, "adjusted for", xname, "\n")}

  # plot 3.  Y residuals vs time with ATS line
  ats.out <- ATS(gam.y$residuals, dat.nonas[,5], dat.nonas[,3], x.cens, LOG = FALSE, xlabel = tname, ylabel = y.txt,printstat = printstat,drawplot=drawplot)
  dat.out <- data.frame(gam.y$residuals, dat.nonas[,2])
  colnames(dat.out) <- c("GAMresidual", xname)
  dat.out <- list(dat.out, ats.out)
  return (invisible(dat.out))
}


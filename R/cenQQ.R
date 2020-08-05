#' Q-Q Plot censored data
#'
#' @description Plots quantile-quantile (Q-Q) plot of censored data fitted a data distribution
#' @param y.var The column of `y` (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param dist One of three distributional shapes to fit to your data:  lognormal (`lnorm`), normal (`norm`) or gamma (`gamma`).
#' @param Yname Optional – input text in quotes to be used as the variable name on the Q-Q plot.  The default is the name of the `y.var` input variable.
#' @importFrom EnvStats gofTestCensored distChooseCensored qqPlotCensored
#' @export
#'
#' @return A single Q-Q plot of data fitted to normal, lognormal or gamma distributions with Shapiro-Francia W value printed on plot.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#'
#' data(Brumbaugh)
#' cenQQ(Brumbaugh$Hg,Brumbaugh$HgCen)
#'
#' # User defined distribution
#' cenQQ(Brumbaugh$Hg,Brumbaugh$HgCen,dist="gamma")

cenQQ <- function(y.var, cen.var, dist = "lnorm", Yname = yname)  {
  #added to stop if dist is not from the list
  if(!(dist%in%c("norm","lnorm","gamma"))){stop(paste0(dist," distribution is not supported with this function, try again."))}

  yname <- deparse(substitute(y.var))
  cen.logical <- as.logical(cen.var)
  var.choose <- distChooseCensored(y.var, cen.logical)
  norm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$norm$statistic, 3) )
  lnorm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$lnorm$statistic, 3) )
  gamma.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$gamma$statistic, 3) )

  oldpar<- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (dist == "norm") {
    qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
    mtext(norm.text)
  }

  if (dist == "lnorm")  {
    ylabel <- paste ("ln (", Yname, ")", sep = "")
    qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = ylabel, main = "Lognormal Q-Q Plot")
    mtext(lnorm.text)
  }

  if (dist == "gamma")  {
    qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
    mtext(gamma.text)
  }
}

#' Q-Q Plot censored data
#'
#' @description Plots a quantile-quantile (Q-Q) plot of censored data versus a fitted data distribution
#' @param y.var The column of `y` (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param dist One of three distributional shapes to fit to your data:  lognormal (`lnorm`), normal (`norm`) or gamma (`gamma`).
#' @param Yname Optional – input text in quotes to be used as the variable name on the Q-Q plot.  The default is the `Yname` name of the `y.var` input variable.
#' @export
#'
#' @return A single Q-Q plot of data fitted by normal, lognormal or gamma distributions with Shapiro-Francia W value printed on plot.
#'
#' @importFrom EnvStats distChoose gofTestCensored qqPlotCensored distChooseCensored qqPlot
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#'\donttest{
#' data(Brumbaugh)
#' cenQQ(Brumbaugh$Hg,Brumbaugh$HgCen)
#'
#' # User defined distribution
#' cenQQ(Brumbaugh$Hg,Brumbaugh$HgCen,dist="gamma")
#'}



cenQQ <- function(y.var, cen.var, dist = "lnorm", Yname = yname)  {
  #added to stop if dist is not from the list
  #if(!(dist%in%c("norm","lnorm","gamma"))){stop(paste0(dist," distribution is not supported with this function, try again."))}

  yname <- deparse(substitute(y.var))
  cen.logical <- as.logical(cen.var)

  if (sum(as.integer(cen.var)) > 0)    # not all data are detects
  {var.choose <- distChooseCensored(y.var, cen.logical)
  norm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$norm$statistic, 3) )
  lnorm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$lnorm$statistic, 3) )
  gamma.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$gamma$statistic, 3) )

  if (dist == "norm") {
    qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
    mtext(norm.text)
    #  legend("bottomright", legend = norm.text)
  }

  if (dist == "lnorm")  {
    ylabel <- paste ("ln (", Yname, ")", sep = "")
    qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = ylabel, main = "Lognormal Q-Q Plot")
    mtext(lnorm.text)
    #  legend("bottomright", legend = lnorm.text)
  }

  if (dist == "gamma")  {
    qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
    mtext(gamma.text)
    #   legend("bottomright", legend = gamma.text)
  }
  }

  else    # all data are detects
  {var.choose <- distChoose(y.var, method = "sf", alpha = 0.05, choices = c("norm", "gamma", "lnorm"))
  norm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$norm$statistic, 3) )
  lnorm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$lnorm$statistic, 3) )
  gamma.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$gamma$statistic, 3) )

  if (dist == "norm") {
    qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
    mtext(norm.text)
  }

  if (dist == "lnorm")  {
    ylabel <- paste ("ln (", Yname, ")", sep = "")
    qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = ylabel, main = "Lognormal Q-Q Plot")
    mtext(lnorm.text)
  }

  if (dist == "gamma")  {
    qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
    mtext(gamma.text)
  } }
}

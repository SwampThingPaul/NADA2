#' Censored Q-Q Plot comparison
#'
#' @description Produces three quantile-quantile (Q-Q) plot, also called a probability plot based on three distribution (normal, lognormal and gamma distributions).
#' @param y.var The column of y (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param Yname Optional – input text in quotes to be used as the variable name on all plots.  The default is the name of the `y.var` input variable.
#' @importFrom EnvStats distChooseCensored qqPlotCensored
#' @importFrom EnvStats gofTestCensored
#' @export
#' @return Plots three Q-Q plots based on normal, lognormal and gamma distributions and prints the best-fit distribution.
#' @details Produces three Q-Q plots and reports which has the highest Shapiro-Francia test statistic (W).  The distribution with the highest W is the best fit of the three.
#'
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#' data(Brumbaugh)
#' cenCompareQQ(Brumbaugh$Hg,Brumbaugh$HgCen)

cenCompareQQ <- function(y.var, cen.var, Yname = yname)  {
  yname <- deparse(substitute(y.var))
  cen.logical <- as.logical(cen.var)
  var.choose <- distChooseCensored(y.var, cen.logical)
  norm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$norm$statistic, 3) )
  lnorm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$lnorm$statistic, 3) )
  gamma.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$gamma$statistic, 3) )
  all.W <- c(var.choose$test.results$norm$statistic, var.choose$test.results$lnorm$statistic, var.choose$test.results$gamma$statistic)
  all.W <- all.W/ max(all.W)
  best.text <- c("normal", "lognormal", "gamma")
  max.distrib <- best.text[all.W==1.0]

  if (var.choose$decision != "Nonparametric") {
    best.dist <- paste (var.choose$decision, "is a good fit")}
  else { best.dist <- paste ("Best of the three distributions is the", max.distrib)
  }
  cat(best.dist, "\n")
  oldpar<- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow=c(2,2))
  qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
  mtext(norm.text)
  #  legend("bottomright", legend = norm.text)

  qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = Yname, main = "Lognormal Q-Q Plot")
  mtext(lnorm.text)
  #  legend("bottomright", legend = lnorm.text)

  qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
  mtext(gamma.text)
  #   legend("bottomright", legend = gamma.text)

  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.47, y = 0.6, best.dist, pos = 1, cex = 1.2, col = "black", family="sans", font=1, adj=1)

}

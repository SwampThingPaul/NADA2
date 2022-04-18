#' Censored Q-Q Plot comparison
#'
#' @description Produces three quantile-quantile (Q-Q) plots, also called probability plots, based on three distributions (normal, lognormal and gamma distributions).
#' @param y.var The column of y (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param Yname Optional – input text in quotes to be used as the variable name on all plots.  The default `Yname` is the name of the `y.var` input variable.
#' @param printrslt Logical `TRUE`/`FALSE` option of whether to print the best distribution in the console window, or not.  Default is `TRUE.`
#' @export
#' @details Produces three Q-Q plots and reports which one has the highest Shapiro-Francia test statistic (W).  The distribution with the highest W is the best fit of the three.
#' @return Plots three Q-Q plots based on normal, lognormal and gamma distributions and prints the best-fit distribution.
#'
#' @importFrom EnvStats distChooseCensored distChoose gofTestCensored qqPlotCensored qqPlot
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Millard, S.P., 2013. EnvStats: An R Package for Environmental Statistics. Springer-Verlag, New York.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#'
#' data(Brumbaugh)
#' \donttest{cenCompareQQ(Brumbaugh$Hg,Brumbaugh$HgCen)}


cenCompareQQ <- function(y.var, cen.var, Yname = yname,printrslt=TRUE)  {
  yname <- deparse(substitute(y.var))
  if (sum(as.integer(cen.var)) > 0)    # not all data are detects

  { cen.logical <- as.logical(cen.var)
  var.choose <- EnvStats::distChooseCensored(x=y.var, censored=cen.logical)
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
  if(printrslt==TRUE){cat(best.dist, "\n")}

  par(mfrow=c(2,2))
  EnvStats::qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
  mtext(norm.text)
  #  legend("bottomright", legend = norm.text)

  EnvStats::qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = Yname, main = "Lognormal Q-Q Plot")
  mtext(lnorm.text)
  #  legend("bottomright", legend = lnorm.text)

  EnvStats::qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
  mtext(gamma.text)
  #   legend("bottomright", legend = gamma.text)

  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.47, y = 0.6, best.dist, pos = 1, cex = 1.2, col = "black", family="sans", font=1, adj=1)
  par(mfrow=c(1,1))
  }

else  # all data are detects
{ cen.logical <- as.logical(cen.var)
var.choose <- EnvStats::distChoose(x=y.var, method = "sf", alpha = 0.05, choices = c("norm", "gamma", "lnorm"))
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

if(printrslt==TRUE){cat(best.dist, "\n")}

par(mfrow=c(2,2))
EnvStats::qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
mtext(norm.text)

EnvStats::qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = Yname, main = "Lognormal Q-Q Plot")
mtext(lnorm.text)

EnvStats::qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
mtext(gamma.text)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.47, y = 0.6, best.dist, pos = 1, cex = 1.2, col = "black", family="sans", font=1, adj=1)
par(mfrow=c(1,1))
}
}


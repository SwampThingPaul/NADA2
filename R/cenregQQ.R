#' Q-Q plot of censored regression residuals
#'
#' @description Plots a quantile-quantile (Q-Q) plot of censored regression residuals for simple or multiple regression.
#' @param y.var The column of `y` (response variable) values plus detection limits. Alternatively, with interval-censord data, the column of the lower end of the interval.
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.  Alternatively, with interval-censored data the column of the high end of the interval.
#' @param x.vars One or more uncensored explanatory variable(s). See Details
#' @param LOG Indicator of whether to compute the regression in the original y units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @param intcens a logical value indicating the input data is interval-censored instead of a column of values plus a column of indicators.
#' @param main overall title for the plot. A default titel will be used if none is specified.
#' @export
#' @return Q-Q Plot of model residuals and Shapiro-Francia test results.
#'
#' @importFrom survival survreg Surv
#' @importFrom EnvStats qqPlotCensored gofTestCensored
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Millard, S.P., 2013. EnvStats: An R Package for Environmental Statistics. Springer-Verlag, New York.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#' data(Brumbaugh)
#'
#' # One variable
#' cenregQQ(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh$PctWetland)
#'
#' # More than one variable for demostration purposes
#' cenregQQ(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh[,c("PctWetland","SedLOI","Weight")])

cenregQQ <- function(y.var, cen.var, x.vars, LOG = TRUE, intcens = FALSE, main = NULL) {
  yname <- deparse(substitute(y.var))
  if (LOG == TRUE)  {lnvar <- log(y.var)
  flip.log <- max(lnvar) +1 - lnvar
  surv.log <- Surv(flip.log, as.logical(1-cen.var), type="right" )
  if (is.data.frame(x.vars))  {
    reg.out <- survreg(surv.log ~ ., data = x.vars, dist = "gaussian")

  }
  else{reg.out <- survreg(surv.log~x.vars, dist = "gaussian") }    # 1 x variable
  newcoeffs <- reg.out$coefficients * (-1)
  newcoeffs[1] <- max(lnvar)+1 + newcoeffs[1]
  ylog.pred <- max(lnvar) +1 - reg.out$linear.predictors
  ylog.resi <- lnvar - ylog.pred
  vtext<- paste("Quantiles of", yname, "residuals (log units)")
  testnorm <- gofTestCensored(ylog.resi,cen.var)
  if (is.null(main)) main <- "Lognormal Q-Q Plot of Residuals"
  ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
  qqPlotCensored(ylog.resi, cen.var, add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = main)
  mtext(ptext)
  }
  else{
    if (intcens == TRUE) {y.low <- y.var*(1-cen.var)
    surv.norm <- Surv(y.low, y.var, type="interval2")}
    else {flip.norm <- max(y.var) +1 - y.var
    surv.norm <- Surv(flip.norm, as.logical(1-cen.var) )}

    if (is.data.frame(x.vars))  {
      reg.out <- survreg(surv.norm ~ ., data = x.vars, dist = "gaussian")
    }

    else{reg.out <- survreg(surv.norm~x.vars, dist = "gaussian") }
    if (intcens == TRUE) {ynorm.pred <- reg.out$linear.predictors
    ynorm.resi <- y.low - ynorm.pred}
    else{ynorm.pred <- max(y.var) +1 - reg.out$linear.predictors
    ynorm.resi <- y.var - ynorm.pred}
    vtext<- paste("Quantiles of", yname, "residuals")
    testnorm <- gofTestCensored(ynorm.resi,cen.var)
    if (is.null(main)) main <- "Normal Q-Q Plot of Residuals"
    ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic, 5), "  p =", round(testnorm$p.value, 5))
    qqPlotCensored(ynorm.resi, cen.var, add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = main)
    mtext(ptext)
  }
}

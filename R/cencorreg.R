#' Correlation and Regression with censored data
#'
#' @description Computes three parametric correlation coefficients for one X variable and the corresponding R squared for multiple X variables, and a regression equation for censored data.
#' @param y.var The column of y (response variable) values plus detection limits.
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value is in `y.var`.
#' @param x.vars One or more uncensored explanatory variable(s). For multiple variables it must be a data frame of numeric, character and factor variables. See Details
#' @param LOG Indicator of whether to compute the regression in the original y units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @param verbose default `verbose=2`, see details.
#' @param pred.plot default is FALSE. Produces a plot of data and regression model predictions.  To do this the first (or only) x variable in the X dataframe must be a continuous (not factor) variable, and it becomes the x variable in the plot.
#' @param pred.col  default is "purple".  Changes the color of the predicted lines in the prediction plot.
#' @export
#' @return
#' When `x.vars` is one variable, likelihood, rescaled likelihood and McFaddens correlation coefficient (`R`) are printed.
#' When `x.vars` is a `data.frame` of more than one variable, likelihood, rescaled likelihood and McFaddens coefficent of determination (`R2`) are printed.
#'
#' Model coefficients (intercept and slopes), Chi-Squared statistic and p-value for the test that all slope coefficients equal zero (overall test), and model AIC and BIC are provided.
#'
#' A Q-Q plot of model residuals with corresponding Shapiro-Francia W and p-value are plotted for evaluation of model distributional assumptions when `verbose=2` (the default).
#'
#' @importFrom survival survreg Surv
#' @importFrom EnvStats gofTestCensored qqPlotCensored
#'
#' @details
#'
#' `x.vars`: If one x variable only, enter its name.  If multiple x variables, enter the name of a data frame of columns of the x variables. Only columns used as `X` variables in the regression are allowed. Create this by `x.frame <- data.frame (Temp, Flow, Time)` for 3 variables (temperature, flow and time) used as the `X` variables in the regression.  To produce a pred.plot plot of predicted values the first variable in the array must be a continuous (not a factor) variable.
#'
#' AIC and BIC are printed to help evaluate the ‘best’ regression model.  Lower values are better when comparing models with the same Y units and same data.Cannot be used to compare models with differing Y units (such as Y~X versus logY~X). Can be used to compare models with differing X units such as Y~X vs Y~logX.
#'
#' The default Y units are that the Y variable will be log transformed.  Change this with the LOG = option, setting LOG = FALSE.
#'
#' `verbose` option. Default is 2 which provides full output in the console and qqplots in a graphics window. A value of 1 only provides partial results in the console and no plots. A value of 0 provides no output; the returning computations will be stored in the specified object.
#'
#' The Y parameter in the model output (modelname$y) is -1 times the data that were input. This is due to the shift from left censored data to the required right censored data of the survreg function.  The original data are not changed and are used to draw the pred.plot.
#'
#' @seealso [survival::survreg]
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' @examples
#'
#' data(Brumbaugh)
#'
#' # One variable
#' cencorreg(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh$SedMeHg)
#'
#' # One variable with pred.plot=T
#' cencorreg(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh$SedMeHg,pred.plot=TRUE)
#'
#' # More than one variable for demostration purposes
#'cencorreg(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh[,c("SedMeHg","PctWetland")])
#'



cencorreg <- function(y.var, cen.var, x.vars, LOG = TRUE, verbose = 2, pred.plot = FALSE, pred.col = "purple")
{
  yname <- deparse(substitute(y.var))
  nonas <- na.omit(data.frame(y.var, cen.var, x.vars, stringsAsFactors = TRUE))
  xnona <- nonas[,-(1:2)]      # data frame of x variables

  # default verbose = 2 prints detailed output and QQ plots.
  #         verbose = 1 prints only one line of output.  No plots.
  #         verbose = 0 prints no output or plots.

  if (LOG == TRUE)  {lnvar <- log(nonas[,1])    # take logs of Y (default)
  flip.log <- -1 * lnvar
  surv.log <- Surv(flip.log, as.logical(1-nonas[,2]) )
  yplot.name <- paste("log(", yname, ")", sep = "")

  if (is.data.frame(xnona))  {             # multiple x variables
    reg.out <- survreg(surv.log ~ ., data = xnona, dist = "gaussian")
    cn <- names(reg.out$coefficients[-1])
    xvars.txt <- cn[1]
    for (i in 1:length(cn))  {j <-(i+1)
    if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
    }
    reg.out$call[3] <- xvars.txt
  }

  else { xname <- deparse(substitute(x.vars))          # 1 x variable
  xnona <- as.data.frame(xnona)
  names(xnona) <- xname
  reg.out <- survreg(surv.log~., data = xnona, dist = "gaussian")
  cn <- names(reg.out$coefficients[-1])
  reg.out$call[3] <- xname
  }

  ylog.pred <- -1 * reg.out$linear.predictors
  reg.out$call[2] <- paste("log(", yname, ")", sep = "")
  reg.out$coefficients <- reg.out$coefficients * (-1)
  ylog.resi <- lnvar - ylog.pred
  reg.out$linear.predictors <- ylog.pred
  reg.out$resids <- ylog.resi
  # reg.out$y <- reg.out$y * (-1)
  vtext<- paste("Quantiles of", yname, "residuals (log units)")

  if (verbose == 2) {
    testnorm <- gofTestCensored(ylog.resi,nonas[,2])
    ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
    qqPlotCensored(ylog.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Lognormal Q-Q Plot of residuals")
    mtext(ptext)
  } }  # end of taking logs of Y

  else                                          #  Y in original units, normal Q-Q plot
  { yplot.name <- yname
  if(min(nonas[,1] >= 0)) { y.low <- nonas[,1]*(1-nonas[,2])   #  0 for low end of all NDs
  surv.norm <- Surv(y.low, nonas[,1], type="interval2")
  if (is.data.frame(xnona))  {             # multiple x variables
    reg.out <- survreg(surv.norm ~ ., data = xnona, dist = "gaussian")
    cn <- names(reg.out$coefficients[-1])
    xvars.txt <- cn[1]
    for (i in 1:length(cn))  {j <-(i+1)
    if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
    }
    reg.out$call[3] <- xvars.txt
  }

  else { xname <- deparse(substitute(x.vars))          # 1 x variable
  xnona <- as.data.frame(xnona)
  names(xnona) <- xname
  reg.out <- survreg(surv.norm~., data = xnona, dist = "gaussian")
  cn <- names(reg.out$coefficients[-1])
  reg.out$call[3] <- xname
  }

  reg.out$call[2] <- yname
  ynorm.pred <- reg.out$linear.predictors
  ynorm.resi <- y.low - ynorm.pred
  # reg.out$linear.predictors <- ynorm.pred
  reg.out$resids <- ynorm.resi
  # reg.out$y <- reg.out$y * (-1)
  vtext<- paste("Quantiles of", yname, "residuals")

  if (verbose == 2) {
    testnorm <- gofTestCensored(ynorm.resi,nonas[,2])
    ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
    qqPlotCensored(ynorm.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Normal Q-Q Plot of residuals")
    mtext(ptext)
  } }  # end of low censored = 0

  #  negative y values.  Use flip variable to set -Inf as low end for censored values
  else{ flip.norm <- -1 * nonas[,1]
  surv.norm <- Surv(flip.norm, as.logical(1-nonas[,2]) )

  if (is.data.frame(xnona))  {             # multiple x variables
    reg.out <- survreg(surv.norm ~ ., data = xnona, dist = "gaussian")
    cn <- names(reg.out$coefficients[-1])
    xvars.txt <- cn[1]
    for (i in 1:length(cn))  {j <-(i+1)
    if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
    }
    reg.out$call[3] <- xvars.txt
  }

  else { xname <- deparse(substitute(x.vars))          # 1 x variable
  xnona <- as.data.frame(xnona)
  names(xnona) <- xname
  reg.out <- survreg(surv.norm~., data = xnona, dist = "gaussian")
  cn <- names(reg.out$coefficients[-1])
  reg.out$call[3] <- xname
  }

  reg.out$call[2] <- yname
  ynorm.pred <- -1 * reg.out$linear.predictors
  reg.out$coefficients <- reg.out$coefficients * (-1)
  ynorm.resi <- nonas[,1] - ynorm.pred
  reg.out$linear.predictors <- ynorm.pred
  reg.out$resids <- ynorm.resi
  # reg.out$y <- reg.out$y * (-1)
  vtext<- paste("Quantiles of", yname, "residuals")

  if (verbose == 2) {
    testnorm <- gofTestCensored(ynorm.resi,nonas[,2])
    ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
    qqPlotCensored(ynorm.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Normal Q-Q Plot of residuals")
    mtext(ptext)
  }
  }   # end of flipping
  } # end of Y in original units

  if (is.data.frame(xnona))       {   # multiple x variables.  Print r-squared.
    LRr2 <- signif(1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y)),4)
    McFr2 <- 1-(reg.out$loglik[2]/reg.out$loglik[1])
    Nag.r2 <- signif((1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))) /(1-exp(2*reg.out$loglik[1]/length(reg.out$y))),4)
    AIC <- -2*reg.out$loglik[2] + (2*reg.out$df +1)
    BIC <- -2*reg.out$loglik[2] + log(reg.out$df+reg.out$df.residual+1)*reg.out$df
    if (verbose > 0) cat(" Likelihood R2 =", LRr2, "                   ", "AIC =", AIC,"\n")
    if (verbose == 2) cat(" Rescaled Likelihood R2 =", Nag.r2, "           ", "BIC =", BIC, "\n", "McFaddens R2 =", McFr2, "\n","\n")
  }
  else {                                    # 1 x variable.  Print correlation coefficients
    LRcorr <- signif(sign(reg.out$coefficients[2])*sqrt(1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))),4)
    McFcorr <- signif(sign(reg.out$coefficients[2])*sqrt(1-reg.out$loglik[2]/reg.out$loglik[1]),4)
    Nag.cor <- signif(sign(reg.out$coefficients[2])*sqrt((1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))) /(1-exp(2*reg.out$loglik[1]/length(reg.out$y)))),4)
    AIC <- -2*reg.out$loglik[2] + (2*reg.out$df +1)
    BIC <- -2*reg.out$loglik[2] + log(reg.out$df+reg.out$df.residual+1)*reg.out$df
    if (verbose > 0)  cat(" Likelihood R =", LRcorr, "                   ", "AIC =", AIC,"\n")
    if (verbose == 2) cat("Rescaled Likelihood R =", Nag.cor, "           ", "BIC =", BIC, "\n", "McFaddens R =", McFcorr, "\n","\n")
  }
  if (pred.plot == TRUE & verbose == 2) {
    if (LOG == TRUE)  {z.1 <- data.frame(xnona[,1], reg.out$linear.predictors, lnvar, nonas[,2]) }
    else {z.1 <- data.frame(xnona[,1], reg.out$linear.predictors, nonas[,1], nonas[,2])}   # predicted values and original data
    z <- z.1[order(z.1[,1]),]     # data in order from low to high xnona[,1], the first (or only) x variable
    kenplot(z[,3], z[,4],z[,1], xcen = rep(0, times=length(reg.out$linear.predictors)), Title = "Data and Predictions from Censored Regression Model", xnam=cn[1], ynam=yplot.name, atsline = FALSE)
    lines(z, col = pred.col)
    mtext(paste("Predicted values in", pred.col))
  }
  return(reg.out)                   # returns reg.out object if function assigned to object
}

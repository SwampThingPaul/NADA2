#' Correlation and Regression with censored data
#'
#' @description Computes three parametric correlation coefficients for one X variable and the corresponding R2 for multiple X variables, and a regression equation for censored data.
#' @param y.var The column of y (response variable) values plus detection limits.
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value is in `y.var`.
#' @param x.vars One or more uncensored explanatory variable(s). See Details
#' @param LOG Indicator of whether to compute the regression in the original y units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @export
#' @return
#' When `x.vars` is one variable, likelihood, rescaled likelihood and McFaddens correlation coefficient (`R`) are printed.
#' When `x.vars` is more than one variable, likelihood, rescaled likelihood and McFaddens coefficent of determination (`R2`) are printed.
#'
#' Model coefficients (intercept and slopes), Chi-Squared statistic and p-value for test that all slope coefficients equal zero (overall test), and model AIC and BIC are provided.
#'
#' Q-Q plot of model residuals with corresponding Shapiro-Francia W and p-value are plotted for evaluation of model distributional assumptions.
#'
#' @importFrom survival survreg Surv
#' @importFrom EnvStats gofTestCensored qqPlotCensored
#'
#' @details
#'
#' `x.vars`: If 1 x variable only, enter its name.  If multiple x variables, enter the name of a data frame of columns of the x variables. No extra columns unused in the regression allowed. Create this by `x.frame <- data.frame (Temp, Flow, Time)` for 3 variables (temperature, flow and time).
#'
#' AIC and BIC are printed to help evaluate the ‘best’ regression model.
#'
#' The default is that the Y variable will be log transformed.
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
#' # More than one variable for demostration purposes
#'cencorreg(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh[,c("SedMeHg","PctWetland")])


cencorreg <- function(y.var, cen.var, x.vars, LOG = TRUE, verbose = TRUE) {
  yname <- deparse(substitute(y.var))
  nonas <- na.omit(cbind(y.var, cen.var, x.vars))
  xnona <- nonas[,-(1:2)]

  if (LOG == TRUE)  {lnvar <- log(nonas[,1])   # take logs of Y (default)
  flip.log <- max(lnvar) +1 - lnvar
  #  print(max(lnvar)+1)
  surv.log <- Surv(flip.log, as.logical(1-nonas[,2]) )

  if (is.data.frame(x.vars))  {             # multiple x variables
    reg.out <- survreg(surv.log ~ ., data = xnona, dist = "gaussian")
    cn <- names(reg.out$coefficients[-1])
    xvars.txt <- cn[1]
    for (i in 1:length(cn))  {j <-(i+1)
    if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
    }
    reg.out$call[3] <- xvars.txt
  }

  else { xname <- deparse(substitute(x.vars))          # 1 x variable
  x.df <- as.data.frame(xnona)
  names(x.df) <- xname
  reg.out <- survreg(surv.log~., data = x.df, dist = "gaussian")
  cn <- names(reg.out$coefficients[-1])
  reg.out$call[3] <- cn
  }

  ylog.pred <- max(lnvar) +1 - reg.out$linear.predictors
  reg.out$call[2] <- paste("log(", yname, ")", sep = "")
  reg.out$coefficients <- reg.out$coefficients * (-1)
  reg.out$coefficients[1] <- max(lnvar)+1 + reg.out$coefficients[1]   # coeffs in cenreg lognormal
  ylog.resi <- lnvar - ylog.pred
  reg.out$linear.predictors <- ylog.pred
  reg.out$resids <- ylog.resi
  vtext<- paste("Quantiles of", yname, "residuals (log units)")

  if (verbose == TRUE) {
    testnorm <- gofTestCensored(ylog.resi,nonas[,2])
    ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
    qqPlotCensored(ylog.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Lognormal Q-Q Plot of residuals")
    mtext(ptext)
  } }  # end of taking logs of Y

  else                                          #  Y in original units, normal Q-Q plot
  { if(min(nonas[,1] >= 0)) { y.low <- nonas[,1]*(1-nonas[,2])}   #  0 for low end of all NDs
    surv.norm <- Surv(y.low, nonas[,1], type="interval2")
    if (is.data.frame(x.vars))  {       # multiple x variables
      reg.out <- survreg(surv.norm ~ ., data = xnona, dist = "gaussian")
      cn <- names(reg.out$coefficients[-1])
      xvars.txt <- cn[1]
      for (i in 1:length(cn))  {j <-(i+1)
      if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
      }
      reg.out$call[3] <- xvars.txt
    }

    else { xname <- deparse(substitute(x.vars))          # 1 x variable
    x.df <- as.data.frame(xnona)
    names(x.df) <- xname
    reg.out <- survreg(surv.norm~., data = x.df, dist = "gaussian")
    cn <- names(reg.out$coefficients[-1])
    reg.out$call[3] <- cn
    }

    reg.out$call[2] <- yname
    ynorm.pred <- reg.out$linear.predictors
    ynorm.resi <- y.low - ynorm.pred
    reg.out$linear.predictors <- ynorm.pred
    reg.out$resids <- ynorm.resi
    vtext<- paste("Quantiles of", yname, "residuals")

    if (verbose == TRUE) {
      testnorm <- gofTestCensored(ynorm.resi,nonas[,2])
      ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
      qqPlotCensored(ynorm.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Normal Q-Q Plot of residuals")
      mtext(ptext)
    } }  # end of low censored = 0

  #  negative y values.  Use flip variable to set -Inf as low end for censored values
  else{ flip.norm <- max(nonas[,1]) +1 - nonas[,1]
  surv.norm <- Surv(flip.norm, as.logical(1-nonas[,2]) )

  if (is.data.frame(x.vars))  {       # multiple x variables
    reg.out <- survreg(surv.norm ~ ., data = xnona, dist = "gaussian")
    cn <- names(reg.out$coefficients[-1])
    xvars.txt <- cn[1]
    for (i in 1:length(cn))  {j <-(i+1)
    if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
    }
    reg.out$call[3] <- xvars.txt
  }

  else { xname <- deparse(substitute(x.vars))          # 1 x variable
  x.df <- as.data.frame(xnona)
  names(x.df) <- xname
  reg.out <- survreg(surv.norm~., data = x.df, dist = "gaussian")
  cn <- names(reg.out$coefficients[-1])
  reg.out$call[3] <- cn
  }

  reg.out$call[2] <- yname
  ynorm.pred <- max(nonas[,1]) +1 - reg.out$linear.predictors
  reg.out$coefficients <- reg.out$coefficients * (-1)
  reg.out$coefficients[1] <- max(nonas[,1]) +1 + reg.out$coefficients[1]   # coeffs in orig scale
  ynorm.resi <- nonas[,1] - ynorm.pred
  reg.out$linear.predictors <- ynorm.pred
  reg.out$resids <- ynorm.resi

  vtext<- paste("Quantiles of", yname, "residuals")

  if (verbose == TRUE) {
    testnorm <- gofTestCensored(ynorm.resi,nonas[,2])
    ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
    qqPlotCensored(ynorm.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Normal Q-Q Plot of residuals")
    mtext(ptext)
  }
  }   # end of flipping
} # end of Y in original units

if (is.data.frame(x.vars))  {             # multiple x variables.  Print r-squared.
  LRr2 <- signif(1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y)),4)
  McFr2 <- signif((1-reg.out$loglik[2]/reg.out$loglik[1]),4)
  Nag.r2 <- signif((1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))) /(1-exp(2*reg.out$loglik[1]/length(reg.out$y))),4)
  AIC <- -2*reg.out$loglik[2] + (2*reg.out$df +1)
  BIC <- -2*reg.out$loglik[2] + log(reg.out$df+reg.out$df.residual+1)*reg.out$df
  cat(" Likelihood R2 =", LRr2, "                   ", "AIC =", AIC,"\n")
  if (verbose == TRUE) cat(" Rescaled Likelihood R2 =", Nag.r2, "          ", "BIC =", BIC, "\n", "McFaddens R2 =", McFr2, "\n","\n")
}
else {                                    # 1 x variable.  Print correlation coefficients
  LRcorr <- signif(sign(reg.out$coefficients[2])*sqrt(1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))),4)
  McFcorr <- signif(sign(reg.out$coefficients[2])*sqrt(1-reg.out$loglik[2]/reg.out$loglik[1]),4)
  Nag.cor <- signif(sign(reg.out$coefficients[2])*sqrt((1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))) /(1-exp(2*reg.out$loglik[1]/length(reg.out$y)))),4)
  AIC <- -2*reg.out$loglik[2] + (2*reg.out$df +1)
  BIC <- -2*reg.out$loglik[2] + log(reg.out$df+reg.out$df.residual+1)*reg.out$df
  cat(" Likelihood R =", LRcorr, "                   ", "AIC =", AIC,"\n", "Rescaled Likelihood R =", Nag.cor, "          ", "BIC =", BIC, "\n", "McFaddens R =", McFcorr, "\n","\n")
}
return(reg.out)                   # returns reg.out object if function assigned to object
}

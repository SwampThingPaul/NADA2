#' Partial plots for censored MLE regression
#'
#' @description Draws a partial plot for each X variable in regression of a censored Y variable against multiple X variables.
#' @param y.var The column of y (response variable) values plus detection limits.
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.vars Multiple uncensored explanatory variable(s). See Details
#' @param LOG Indicator of whether to compute the censored regression in the original y units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original Y units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @param smooth.method Method for drawing a smooth on the partial plot.  Options are c("gam", "none"). "gam" is a censored generalized additive model using the cenGAM and mgcv packages.
#' @param gam.method	Method for computing the gam smooth.  See the mgcv package for options.  Default is a thinplate ("tp") spline.  "cs" is another good option.
#' @param multiplot  If TRUE, plots are drawn 6 per page.  If FALSE, all plots are drawn on a separate page.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @export
#' @return
#' When `x.vars` is one variable, a message is printed that partial plots are not possible with only one X variable and execution stops.
#' When `x.vars` is a data frame of more than one variable, partial plots are drawn for each X variable and text is printed comparing AICs for regression using the untransformed X variable with log and cube-root transforms of the X variable, as a supplement to evaluating linearity on the partial plots alone.
#'
#'
#' @importFrom survival survreg Surv
#' @importFrom mgcv gam
#' @importFrom stats lm
#' @importFrom graphics title
#'
#' @details
#' Partial plots for uncensored data often are drawn with superimposed smooths. At times looking only at the data values without a smooth can better enable the human eye to determine whether the overall pattern is linear or not.  If this is the best method for you, use the smooth.method = "none" option to not draw a smooth.  The most common smooth used for uncensored data is loess, which does not recognize censored data and so uses the detection limit (DL) value itself.  This results in biased-high smooths that incorrectly treat values at the DLs equal to uncensored (detected) data. The partplots function in NADA2 was written to provide a better alternative, smoothing the partial residual pattern with a censored generalized additive model (gam).  The censored gam recognizes the nondetects as left-censored data with a maximum at the DL when computing the smooth. DLs may vary with each observation -- multiple DLs in a dataset are not a problem in routines of the NADA2 package.
#'
#' 'y.var': The default is that the Y variable will be log transformed.
#'
#' `x.vars`: Enter the name of a data frame of columns of the x variables. No extra columns unused in the regression allowed. Create this by `x.frame <- data.frame (Temp, Flow, Time)` for 3 variables (temperature, flow and time).
#'
#' Gray open circles represent censored data and are the residual between the detection limit and the predicted value from the censored regression.  The GAM recognizes that the detection limit is an upper limit, predicted values on the regression line are most often below the detection limit, leading to positive residuals.  Note that the true residual for censored data could be anywhere below the plotted value.  That fact is recognized by the censored GAM but is difficult to represent on a plot.
#'
#' AIC for regression models with un-transformed X, log and cube-root transforms of X are printed to evaluate which of the three transformations results in the ‘best’ model.
#'

#' @seealso [survival::survreg]
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Cook, R.D., 1993. Exploring Partial Residual Plots, Technometrics 35, 351-362.
#'
#' @examples
#'
#' data(Brumbaugh)
#'
#' # For demostration purposes
#' partplots (Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh[,c("SedMeHg","PctWetland")])

partplots <- function (y.var, cen.var, x.vars, LOG = TRUE, smooth.method = "gam", gam.method = "tp", multiplot = TRUE,printstat=TRUE)
{  yname <- deparse(substitute(y.var))
nonas <- na.omit(cbind(y.var, cen.var, x.vars))
xnona <- data.frame(nonas[,-(1:2)])
pch.symb <- (nonas[, 2] + 16) - nonas[, 2] *16
col.symb <- (nonas[, 2] + 1) + nonas[, 2] *6
oldpar<- par(no.readonly=TRUE)
on.exit(par(oldpar))
if (multiplot == TRUE) {par(mfrow = c(3,2))}

# function to put data into format used by the gam function tobit1 model
dat4cengam <- function(y.orig, cen.orig) {
  y.gam <- y.orig*(1-cen.orig)     # low end of possible values. <DL data = 0.
  y.gam [y.gam == 0] <- -Inf       # <DL data = -Inf
  dl.gam <- y.orig*cen.orig         # 0 for detects.  DL for nondetects.
  dat.gam <- invisible(data.frame(y.gam, dl.gam))
  return(dat.gam)
}

if (LOG == TRUE)  {
  if (min(nonas[,1]) <= 0) stop('Cannot take logs of zero or negative values.  Use LOG=FALSE')
  lnvar <- log(nonas[,1])    # Y in log units (default)
  yname.log = paste("ln(", yname, ")", sep = "")
  flip.log <- max(lnvar) - lnvar
  surv.log <- Surv(flip.log, as.logical(1-nonas[,2]) )

  if (is.data.frame(x.vars))  {             # multiple x variables
    yname.part <- paste(yname.log, ".partial", sep="")
    for (i in 1:ncol(xnona))  {temp.x <- data.frame(xnona[, -i])
    xname.part <- paste(colnames(xnona[i]), ".partial", sep="")
    y.part <- survreg(surv.log ~ ., data = temp.x, dist = "gaussian")
    x.part <- lm(xnona[, i] ~ ., data = temp.x)
    flip.resi <- flip.log - y.part$linear.predictors
    # y in original, not flipped, units
    ylog.resi <- -1*flip.resi
    ylog.pred <- max(lnvar) - y.part$linear.predictors
    # y coefficients in original units
    y.coefficients <- y.part$coefficients * (-1)
    y.coefficients[1] <- max(lnvar) + y.part$coefficients[1]

    #  plot y resids vs x resids with smooth

    # option 1:  no smooth
    plot(x.part$residuals, ylog.resi, xlab = xname.part, ylab = yname.part, pch = pch.symb, col = col.symb)

    # option 2:  censored GAM
    if (smooth.method == "gam")
    { data.gam <- dat4cengam(ylog.resi, nonas[,2])
    ylog.gam <- gam(data.gam[,1] ~s(x.part$residuals), bs = gam.method, family = tobit1(link = "identity", left.threshold = data.gam[,2]))
    o <- order(x.part$residuals , ylog.resi)
    lines (x.part$residuals[o], ylog.gam$fitted.values[o], col = 'red')}

    if (multiplot == FALSE) {title(paste("Open circles are", yname, "nondetects"))}

    # compute and print one line (verbose = 1) of r2 and AIC for original and possible x transforms
    if(printstat==TRUE){
    cat(colnames(xnona[i]), "\n", "untransformed", "\n")
    notrans <- cencorreg(nonas[,1], nonas[,2], xnona, verbose = 1)
    AIC.none <- -2*notrans$loglik[2] + (2*notrans$df +1)

    cat("cube root", "\n")
    x.sign <- sign(xnona[i])
    x.cube <- abs(xnona[i])**(1/3) * x.sign
    cube.xnona <- cbind(x.cube, temp.x)
    cube.trans <- cencorreg(nonas[,1], nonas[,2], cube.xnona, verbose = 1)
    AIC.cube <- -2*cube.trans$loglik[2] + (2*cube.trans$df +1)

    cat("log transform", "\n")
    if (min(xnona[i]) > 0) {
      x.log <- log(xnona[i])
      log.xnona <- cbind(x.log, temp.x)
      log.trans <- cencorreg(nonas[,1], nonas[,2], log.xnona, verbose = 1)
      AIC.log <- -2*log.trans$loglik[2] + (2*log.trans$df +1) }
    else {AIC.log = AIC.cube
    cat("Cannot take logs of zero or negative values.", "\n")}
    AIC.diff <- AIC.none - min(AIC.log, AIC.cube)
    cat("Decrease in AIC from transformation of", colnames(xnona[i]), "=", max(0,AIC.diff), "\n", "\n")
    }  # end of cycle thru x variables
}
  }  # end of multiple x vars
  else (stop("For only one x variable partial plots not needed."))
}    # end of logs.

#  In original Y units
else {
  flip.y <- max(nonas[,1]) - nonas[,1]
  surv.y <- Surv(flip.y, as.logical(1-nonas[,2]) )

  if (is.data.frame(x.vars))  {
    yname.part <- paste(yname, ".partial", sep="")
    # multiple x variables
    for (i in 1:ncol(xnona))  {temp.x <- data.frame(xnona[, -i])
    xname.part <- paste(colnames(xnona[i]), ".partial", sep="")
    y.part <- survreg(surv.y ~ ., data = temp.x, dist = "gaussian")
    x.part <- lm(xnona[, i] ~ ., data = temp.x)
    flip.resi <- flip.y - y.part$linear.predictors
    # y in original, not flipped, units
    y.resi <- -1*flip.resi
    y.pred <- max(nonas[,1]) - y.part$linear.predictors
    # y coefficients in original units
    y.coefficients <- y.part$coefficients * (-1)
    y.coefficients[1] <- max(nonas[,1]) + y.part$coefficients[1]

    # option 1: no smooth
    plot(x.part$residuals, y.resi, xlab = xname.part, ylab = yname.part, pch = pch.symb, col = col.symb)

    # option 2:  censored GAM
    if (smooth.method == "gam")
    { data.gam <- dat4cengam(y.resi, nonas[,2])
    y.gam <- gam(data.gam[,1] ~s(x.part$residuals), bs = gam.method, family = tobit1(link = "identity", left.threshold = data.gam[,2]))

    o <- order(x.part$residuals , y.resi)
    lines (x.part$residuals[o], y.gam$fitted.values[o], col = 'red')}

    if (multiplot == FALSE) {title(paste("Open circles are", yname, "nondetects"))}

    # compute and print one line (verbose = 1) of r2 and AIC for original and possible x transforms
    if(printstat==TRUE){
    cat(colnames(xnona[i]), "\n", "untransformed", "\n")
    notrans <- cencorreg(nonas[,1], nonas[,2], xnona, verbose = 1)
    AIC.none <- -2*notrans$loglik[2] + (2*notrans$df +1)

    cat("cube root", "\n")
    x.sign <- sign(xnona[i])
    x.cube <- abs(xnona[i])**(1/3) * x.sign
    cube.xnona <- cbind(x.cube, temp.x)
    cube.trans <- cencorreg(nonas[,1], nonas[,2], cube.xnona, verbose = 1)
    AIC.cube <- -2*cube.trans$loglik[2] + (2*cube.trans$df +1)

    cat("log transform", "\n")
    if (min(xnona[i]) > 0) {
      x.log <- log(xnona[i])
      log.xnona <- cbind(x.log, temp.x)
      log.trans <- cencorreg(nonas[,1], nonas[,2], log.xnona, verbose = 1)
      AIC.log <- -2*log.trans$loglik[2] + (2*log.trans$df +1) }
    else {AIC.log = AIC.cube
    cat("Cannot take logs of zero or negative values.", "\n")}
    AIC.diff <- AIC.none - min(AIC.log, AIC.cube)
    cat("Decrease in AIC from transformation of", colnames(xnona[i]), "=", max(0,AIC.diff), "\n", "\n")
}
    }  # end of cycle thru multiple x variables
  }  # end of multiple variables
  else {stop("For only one x variable partial plots not needed.")}
} # end of original units
par(oldpar)
}

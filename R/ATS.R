#' Akritas-Theil-Sen line for censored data
#'
#' @description
#' Computes Kendall's tau and the Akritas-Theil-Sen (ATS) line for censored data, along with the test that the slope (and Kendall's tau) equal zero.  For one x variable regression.
#' @param y.var The column of y (response variable) values plus detection limits
#' @param y.cen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and `0` (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.var The column of x (explanatory variable) values plus detection limits
#' @param x.cen The x-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `x.var` column, and `0` (or `FALSE`) indicates a detected value in `x.var`.
#' @param LOG Indicator of whether to compute the ATS line in the original y units, or for their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @param retrans Indicator of whether to retransform the plot and line back to original Y-variable units.  Not needed when `LOG = FALSE`. When `retrans = FALSE` & `LOG = TRUE` the plot is drawn with logY units (default). When `retrans =  TRUE` & `LOG = TRUE` the plot is drawn with original Y units.
#' @param xlabel Custom label for the x axis of plots.  Default is x variable column name.
#' @param ylabel Custom label for the y axis of plots.  Default is y variable column name.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @param drawplot Logical `TRUE`/`FALSE` option of whether to draw plots or not. Default is `TRUE`
#'
#' @keywords trend censored
#' @export
#' @return  Coefficients (intercept and slope) for the ATS line are printed, along with Kendall's tau correlation coefficient, test statistic S, and the (single) p-value for the test that tau and slope both equal zero. A scatterplot with the fitted trend-line superimposed is also drawn.
#'
#' @importFrom graphics abline layout legend lines mtext par plot text
#' @importFrom utils data
#' @importFrom stats na.omit
#'
#' @references
#' Akritas, M.G., Murphy, S.A., LaValley, M.P., 1995. The Theil-Sen Estimator With Doubly Censored Data and Applications to Astronomy. Journal of the American Statistical Association 90, 170–177. https://doi.org/10.2307/2291140
#'
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @examples
#' \dontrun{
#' # Both y and x are censored
#' data(PbHeron)
#' with(PbHeron, ATS(Blood, BloodCen, Kidney, KidneyCen))
#'
#' # x is not censored
#' data(Brumbaugh)
#' with(Brumbaugh,ATS(Hg, HgCen, PctWetland))
#' }

ATS <- function(y.var, y.cen, x.var, x.cen = rep(0, length(x.var)),
                LOG = TRUE, retrans = FALSE, xlabel = NULL, ylabel = NULL,
                printstat = TRUE, drawplot = TRUE) {
  yname <- if (!is.null(ylabel)) ylabel else deparse(substitute(y.var))
  xname <- if (!is.null(xlabel)) xlabel else deparse(substitute(x.var))

  alldat <- na.omit(data.frame(y.var, y.cen = as.logical(y.cen),
                               x.var, x.cen = as.logical(x.cen)))
  nobs <- nrow(alldat)

  # Prepare y and shift if needed
  y.raw <- alldat$y.var
  y.shift <- 0
  if (LOG) {
    if (any(y.raw < 0)) stop('Cannot take logs of negative values. Use LOG = FALSE')
    y.raw <- log(y.raw)
    if (any(y.raw < 0)) {
      y.shift <- abs(min(y.raw)) + 1e-4
      y.raw <- y.raw + y.shift
    }
  } else if (any(y.raw < 0)) {
    y.shift <- abs(min(y.raw)) + 1e-4
    y.raw <- y.raw + y.shift
  }

  # Run Akritas-Theil-Sen
  ats <- cenken(y.raw, alldat$y.cen, alldat$x.var, alldat$x.cen)
  slope <- ats$slope
  intercept <- ats$intercept - y.shift
  tau <- ats$tau
  pval <- ats$p
  S <- tau * nobs * (nobs - 1) / 2

  # Subtitle logic
  eqn_lhs <- if (LOG) paste0("ln(", yname, ")") else yname
  slope_fmt <- if (slope >= 0) paste("+", round(slope, 4)) else round(slope, 4)
  subtitle <- paste(eqn_lhs, "=", round(intercept, 4), slope_fmt, "*", xname)

  # Print stats
  if (printstat) {
    cat("Akritas-Theil-Sen line for censored data\n\n")
    cat(subtitle, "\n")
    cat("Kendall's tau =", round(tau, 4), "  p-value =", round(pval, 5), "\n\n")
  }

  # Plot (log or original)
  if (drawplot) {
    x <- range(alldat$x.var)
    y.line <- x * slope + intercept
    if (LOG && !retrans) {
      kenplot(log(alldat$y.var), alldat$y.cen, alldat$x.var, alldat$x.cen,
              atsline = TRUE, xnam = xname, ynam = eqn_lhs)
    } else {
      if (LOG && retrans) {
        xo <- seq(min(alldat$x.var), max(alldat$x.var), length.out = 100)
        yo <- exp(xo * slope + intercept)
        kenplot(alldat$y.var, alldat$y.cen, alldat$x.var, alldat$x.cen,
                atsline = FALSE, xnam = xname, ynam = yname)
        lines(xo, yo, col = "purple")
        subtitle <- paste0(yname, " = e^(", round(intercept,4), ") * e^(",
                           round(slope,4), "*", xname, ")")
      } else {
        kenplot(alldat$y.var, alldat$y.cen, alldat$x.var, alldat$x.cen,
                xnam = xname, ynam = yname)
        lines(x, y.line, col = "purple")
      }
    }
    mtext(subtitle)
  }

  invisible(data.frame(intercept, slope, S, tau, pval))
}


# Older ATS code - saving just in case
# ATS <- function(y.var, y.cen, x.var, x.cen = rep(0, times=length(x.var)),
#                 LOG = TRUE, retrans = FALSE, xlabel = NULL, ylabel = NULL,
#                 printstat = TRUE,drawplot = TRUE){
#   yname <- deparse(substitute(y.var))
#   xname <- deparse(substitute(x.var))
#   y.cen <- as.logical(y.cen)
#   x.cen <- as.logical(x.cen)
#   alldat <-  data.frame(y.var, y.cen, x.var, x.cen)
#   alldat <- na.omit(alldat)
#   if (!is.null(xlabel)) {xname = xlabel}
#   if (!is.null(ylabel)) {yname = ylabel}
#   nobs <- length(alldat[,1])
#
#   if (LOG == TRUE)  {
#     if (min(alldat[,1]) < 0) stop('Cannot take logs of negative values.  Use LOG=FALSE')
#     y.log <- log(alldat[,1])
#     NoShift <- cenken(y.log, alldat[,2], alldat[,3], alldat[,4] )
#     slope <- NoShift$slope
#     intercept <- NoShift$intercept
#     pval <- NoShift$p
#     tau <- NoShift$tau
#     S <- tau*nobs*(nobs-1)*0.5
#
#     # if y.log goes negative, intercept = NA
#     if (min(y.log) < 0) {
#       shift.amt <- abs(min(y.log)) + 0.0001
#       y.shift <- y.log + shift.amt
#       Shifted<-cenken(y.shift, alldat[,2], alldat[,3], alldat[,4])
#       intercept <- Shifted$intercept - shift.amt
#     }
#     ylogname <- paste("ln(", yname, ")", sep = "")
#
#     if(slope>=0){
#       subtitle <- paste (ylogname, "=", round(intercept,4), "+", round(slope,4), "*", xname)
#     }
#     else{
#       subtitle <- paste(ylogname, "=", round(intercept,4), round(slope,4), "*", xname)
#     }
#
#     if(printstat==TRUE){
#     cat ("Akritas-Theil-Sen line for censored data", "\n", "\n")
#     if (slope >= 0) {cat (ylogname, "=", round(intercept,4), "+", round(slope,4), "*", xname, "\n")
#       subtitle
#       short.t <- paste(round(intercept,4), "+", round(slope,4), "*", xname)
#     }
#     else {cat (ylogname, "=", round(intercept,4), round(slope,4), "*", xname, "\n")  # neg slope
#       subtitle
#       short.t <- paste(round(intercept,4), round(slope,4), "*", xname)
#     }
#     cat("Kendall's tau =", round(tau,4), "  p-value =", round(pval, 5), "\n", "\n")
#     }
#     if(drawplot == TRUE){
#     # draw a scatterplot
#     x <- c(min(alldat[,3]), max(alldat[,3]))
#     y <- x*slope+intercept
#     z <- data.frame(x,y)
# #    print(z)
#     kenplot(log(alldat[,1]), alldat[,2], alldat[,3], alldat[,4], atsline = TRUE, xnam=xname, ynam=ylogname)
#  #   lines(z, col = "purple")
#     mtext(subtitle)
#     #
#     }
#     if (retrans == TRUE)
#     { # draw a 2nd scatterplot in original units
#       subtitle <- paste (yname, " = e^(", round(intercept,4), ")",  " * ", "e^(", round(slope,4), " * ", xname,")", sep="" )
#       delta <- seq(1:100)
#       diff<-max(alldat[,3]) - min(alldat[,3])
#       xo = delta *diff/100 + min(alldat[,3])
#       yo = exp(xo*slope)*exp(intercept)
#       z <- data.frame(xo,yo)
#       if(drawplot == TRUE){
#       kenplot(alldat[,1], alldat[,2], alldat[,3], alldat[,4], atsline = FALSE, xnam=xname, ynam=yname)
#       lines(xo, yo, col = "purple")
#       mtext(subtitle)
#       }
#       #
#     }
#     out <- data.frame(intercept, slope, S, tau, pval)
#   }
#   else  # no logarithms used for y
#   {
#     Noshift<-cenken(alldat[,1], alldat[,2], alldat[,3], alldat[,4])
#     slope <- Noshift$slope
#     intercept <- Noshift$intercept
#     pval <- Noshift$p
#     tau <- Noshift$tau
#     S <- tau*nobs*(nobs-1)*0.5
#
#     if (min(alldat[,1]) < 0) {  # y goes negative. Intercept is NA.
#       shift.amt <- abs(min(alldat[,1])) + 0.0001
#       y.shift <- alldat[,1] + shift.amt
#       Shifted<-cenken(y.shift, alldat[,2], alldat[,3], alldat[,4])
#       intercept <- Shifted$intercept - shift.amt
#     }
#
#     if(slope>=0){
#       subtitle <- paste (yname, "=", round(intercept,4), "+", round(slope,4), "*", xname)
#     }
#     else{
#       subtitle <- paste(yname, "=", round(intercept,4), round(slope,4), "*", xname)
#     }
#
#     if(printstat==TRUE){
#     cat ("Akritas-Theil-Sen line for censored data", "\n", "\n")
#     if (slope >= 0) {cat (yname, "=", round(intercept,4), "+", round(slope,4), "*", xname, "\n")
#       }
#     else {cat (yname, "=", round(intercept,4), round(slope,4), "*", xname, "\n")
#       }
#     cat("Kendall's tau =", round(tau,4), "  p-value =", round(pval, 5), "\n", "\n")
#     }
#
#     # draw a scatterplot with kenplot
#     x <- c(min(alldat[,3]),max(alldat[,3]))
#     y <- x*slope+intercept
#     z=data.frame(x,y)
#     if(drawplot == TRUE){
#     kenplot(alldat[,1], alldat[,2], alldat[,3], alldat[,4], xnam=xname, ynam=yname)
#     lines(z, col = "purple")
#     mtext(subtitle)
#     }
#     #
#     out <- data.frame(intercept, slope, S, tau, pval)
#   }
#   return(invisible(out))
# }

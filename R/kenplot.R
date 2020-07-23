#' Plot sensored trend for censored data
#'
#' @description
#' Computes Kendall's tau and the Akritas-Theil-Sen (ATS) line for censored data.  For one x variable regression.
#' @param y1 The column of y (response variable) values plus detection limits
#' @param ycen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and `0` (or `FALSE`) indicates a detected value in `y.var`.
#' @param x1 The column of x (explanatory variable) values plus detection limits
#' @param xcen The x-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `x.var` column, and `0` (or `FALSE`) indicates a detected value in `x.var`.
#' @param retrans Indicator of whether to retransform the plot and line back to original Y-variable units.  Not needed when `LOG = FALSE`. `retrans = FALSE` & `LOG = TRUE` draws the plot in logY units. `retrans =  TRUE` & `LOG = TRUE` draws the plot in original Y units.
#' @param xnam Custom label for the x axis of plots.  Default is x variable column name.
#' @param ynam Custom label for the y axis of plots.  Default is y variable column name.
#' @keywords trend kendall
#' @importFrom NADA cenken
#' @return
#'
#' #' @references
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' @seealso [NADA::cenken]
#'
#' @examples
#' library(NADA) #For example data
#'
#' # Both y and x are censored
#' data(Golden)
#' with(Golden, kenplot(Blood, BloodCen, Kidney, KidneyCen))
#'
#' # x is not censored
#' data(TCEReg)
#' with(TCEReg, kenplot(TCEConc, TCECen, PctIndLU,rep(0, times=length(PctIndLU))))

kenplot <- function(y1,ycen,x1,xcen,retrans = FALSE, xnam=NULL, ynam=NULL)  {
  alldat <- data.frame(y1, ycen, x1, xcen)

  if (!is.null(xnam)) {xnam = xnam}else(xnam=deparse(substitute(x1)))
  if (!is.null(ynam)) {ynam = ynam}else(ynam=deparse(substitute(y1)))

  xmin <- min(alldat[,3])
  xmax <- max(alldat[,3])
  ymin <- min(alldat[,1])
  ymax <- max(alldat[,1])

  if (ymin <= 0) {
    y1a <- y1+abs(min(y1))+1
    cka<- cenken(y1a, as.logical(ycen), x1, as.logical(xcen))
    int <- cka$intercept - (abs(min(y1))+1)
    slp <- cka$slope
    tau <- cka$tau
    pval <- cka$p
  }
  else { ck<- cenken(y1, as.logical(ycen), x1, as.logical(xcen))
  int <- ck$intercept;  slp <- ck$slope
  tau <-ck$tau;  pval<- ck$p
  }

  bothdetect <- as.integer(ycen)+as.integer(xcen)
  detected <- alldat[bothdetect == 0,]
  ynd <- alldat[as.integer(ycen) ==1 ,]  # set of censored y values
  nyc <- length(as.integer(ynd$ycen))    # number of censored y values
  xnd <- alldat[as.integer(xcen) ==1 ,]  # set of censored x values
  nxc <- length(as.integer(xnd$xcen))    # number of censored x values

  plot(detected$x1, detected$y1, ylim = c(ymin, ymax), xlim = c(xmin, xmax), ylab = ynam, xlab = xnam, pch=19, cex=0.7, main="Akritas - Theil - Sen line", xaxs="r", yaxs="r")
  if (retrans == FALSE) {
    abline(int, slp, col = "purple")}
    if (nyc != 0) {
      for (i in 1:nyc ){
        dashy <- min(0-0.5, ymin)
        dashx <- ynd[i,3]
        dash <- data.frame (dashx, dashy)
        dash[2,1] <- ynd[i,3]
        dash[2,2] <- ynd[i,1]
        lines(dash, lty="dashed", col = "red")
      }
    }
    # horizontal dashed lines for x censored
    if (nxc != 0) {
      for (i in 1:nxc ){
        dashy <- xnd[i,1]
        dashx <- min(0, xmin-0.5)
        dash <- data.frame (dashx, dashy)
        dash[2,1] <- xnd[i,3]
        dash[2,2] <- xnd[i,1]
        lines(dash, lty="dashed", col = "red")
      }
    }
  }


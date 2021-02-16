#' Plot robust median ATS line for censored data
#'
#' @description
#' Function used by other functions to plot the Akritas-Theil-Sen (ATS) line for censored data.  Only one x variable allowed. Both Y and X variables may be censored.
#' @param y1 The column of y (response variable) values plus detection limits
#' @param ycen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and `0` (or `FALSE`) indicates a detected value in `y.var`.
#' @param x1 The column of x (explanatory variable) values plus detection limits
#' @param xcen The x-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `x.var` column, and `0` (or `FALSE`) indicates a detected value in `x.var`.
#' @param atsline Indicator of whether to draw the ATS line or not. Default is FALSE.
#' @param xnam Custom label for the x axis of plots.  Default is x variable column name.
#' @param ynam Custom label for the y axis of plots.  Default is y variable column name.
#' @param Title Custom title for plots.  Default is "Akritas - Theil - Sen line".
#' @importFrom NADA cenken
#' @return
#' Scatterplot of data plus ATS line.  Censored values are drawn for both X and Y variables as dashed lines up to the detection limits.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @seealso [NADA::cenken]
#'
#' @examples
#' \dontrun{
#' # Both y and x are censored
#' data(PbHeron)
#' with(PbHeron, kenplot(Blood, BloodCen, Kidney, KidneyCen))
#'
#' # x is not censored
#' data(Brumbaugh)
#' with(Brumbaugh, kenplot(Hg, HgCen, PctWetland,rep(0, times=length(PctWetland))))
#' }

kenplot <- function(y1, ycen, x1, xcen, atsline = FALSE, xnam=NULL, ynam=NULL, Title="Akritas - Theil - Sen line")  {
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

#  oldpar <- par(no.readonly = TRUE)
#  on.exit(par(oldpar))

  plot(detected$x1, detected$y1, ylim = c(ymin, ymax), xlim = c(xmin, xmax), ylab = ynam, xlab = xnam, pch=19, cex=0.7, main=Title, xaxs="r", yaxs="r")
  if (atsline == TRUE) {
    abline(int, slp, col = "purple")}

  # vertial dashed lines for y censored
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


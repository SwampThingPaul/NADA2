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
#' @param ylim argument consistent with `base` plotting to adjust y-axis limits; see `par` for more details
#' @param xlim argument consistent with `base` plotting to adjust x-axis limits; see `par` for more details
#' @param pch argument consistent with `base` plotting to adjust point type, see `par` for more details
#' @param cex argument consistent with `base` plotting to adjust point size; see `par` for more details
#' @param xaxs argument consistent with `base` plotting to adjust x-axis type; see `par` for more details
#' @param yaxs argument consistent with `base` plotting to adjust y-axis type; see `par` for more details
#' @param ... argument to adjust other base plotting functions; see `par` for more details
#' @export
#' @return
#' Scatterplot of data plus ATS line.  Censored values are drawn for both X and Y variables as dashed lines up to the detection limits.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
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

kenplot <- function(y1, ycen, x1, xcen,
                    atsline = FALSE,
                    xnam = NULL, ynam = NULL, Title = "Akritas - Theil - Sen line",
                    ylim = NULL, xlim = NULL, pch = 19, cex = 0.7,
                    xaxs = "r", yaxs = "r", ...) {

  # Assemble data and check for consistency
  alldat <- data.frame(y1 = y1, ycen = ycen, x1 = x1, xcen = xcen)
  if (anyNA(alldat)) stop("Missing values are not supported in kenplot()")

  # Derive names for axes if not provided
  xlab <- if (!is.null(xnam)) xnam else deparse(substitute(x1))
  ylab <- if (!is.null(ynam)) ynam else deparse(substitute(y1))

  # Calculate plot limits
  xlim.val <- if (!is.null(xlim)) xlim else range(x1, na.rm = TRUE)
  ylim.val <- if (!is.null(ylim)) ylim else range(y1, na.rm = TRUE)

  # Fit censored Kendall model
  if (min(y1, na.rm = TRUE) <= 0) {
    offset <- abs(min(y1)) + 1
    y1a <- y1 + offset
    fit <- cenken(y1a, as.logical(ycen), x1, as.logical(xcen))
    int <- fit$intercept - offset
  } else {
    fit <- cenken(y1, as.logical(ycen), x1, as.logical(xcen))
    int <- fit$intercept
  }
  slope <- fit$slope

  # Separate subsets
  bothdetect <- !(as.logical(ycen) | as.logical(xcen))
  detected <- alldat[bothdetect, ]
  ycen_only <- alldat[as.logical(ycen), ]
  xcen_only <- alldat[as.logical(xcen), ]

  # Start plotting
  plot(detected$x1, detected$y1,
       xlim = xlim.val, ylim = ylim.val,
       xlab = xlab, ylab = ylab,
       main = Title, pch = pch, cex = cex,
       xaxs = xaxs, yaxs = yaxs, ...)

  # Optional regression line
  if (atsline) abline(int, slope, col = "purple")

  # Add vertical lines for censored y-values
  if (nrow(ycen_only) > 0) {
    for (i in seq_len(nrow(ycen_only))) {
      lines(c(ycen_only$x1[i], ycen_only$x1[i]),
            c(min(ylim.val[1] - 0.5, 0), ycen_only$y1[i]),
            lty = "dashed", col = "red")
    }
  }

  # Add horizontal lines for censored x-values
  if (nrow(xcen_only) > 0) {
    for (i in seq_len(nrow(xcen_only))) {
      lines(c(min(xlim.val[1] - 0.5, 0), xcen_only$x1[i]),
            c(xcen_only$y1[i], xcen_only$y1[i]),
            lty = "dashed", col = "red")
    }
  }

  invisible(list(slope = slope, intercept = int, tau = fit$tau, p = fit$p))
}


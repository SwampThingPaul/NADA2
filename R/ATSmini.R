#' Kendall's tau and ATS line for censored data
#'
#' @description Computes Kendall's tau and the Akritas-Theil-Sen (ATS) line for censored data.  Is called by censeaken because it is much faster than the ATS function.  It is not expected to be of much use to users on its own. The x variable (time) may not be censored.
#' @param y.var The column of y (response variable) values plus detection limits.
#' @param y.cen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.var The column of x (time for a trend test) values.
#'
#' @return Returns the intercept and slope (ATS line), tau (Kendall's tau), p-value and S-value (test statistic).
#' @export
#' @references
#' Akritas, M.G., Murphy, S.A., LaValley, M.P., 1995. The Theil-Sen Estimator With Doubly Censored Data and Applications to Astronomy. Journal of the American Statistical Association 90, 170–177. https://doi.org/10.2307/2291140
#'
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#' @examples
#' \dontrun{
#' # x may not be censored.  Use the ATS function when x is censored.
#' data(Brumbaugh)
#'
#' with(Brumbaugh, ATSmini(Hg, HgCen, SedLOI))
#' }

ATSmini <- function(y.var, y.cen, x.var) {
  y.cen <- as.logical(y.cen)
  nobs <- length(y.var)

  Noshift<-cenken(y.var, y.cen, x.var)
  slope <- Noshift$slope
  intercept <- Noshift$intercept
  pval <- Noshift$p
  tau <- Noshift$tau
  S <- tau*nobs*(nobs-1)*0.5

  if (min(y.var) < 0) {
    # y goes negative. Intercept is NA.
    shift.amt <- abs(min(y.var)) + 0.0001
    y.shift <- y.var + shift.amt
    Shifted<-cenken(y.shift, y.cen, x.var)
    intercept <- Shifted$intercept - shift.amt
  }

  out <- data.frame(intercept, slope, tau, pval, S)
  return(out)
}

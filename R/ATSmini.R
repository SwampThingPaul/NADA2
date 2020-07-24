#' Kendall's tau and ATS line for censored data
#'
#' @description Computes Kendall's tau and the Akritas-Theil-Sen (ATS) line for censored data.  For one x variable regression.
#' @param y.var The column of y (response variable) values plus detection limits.
#' @param y.cen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.var The column of x (explanatory variable) values plus detection limits.
#' @importFrom NADA cenken
#' @return Returns the intercept and slope (ATS line), tau (Kendall's tau), p-value and S-value (test statistic)
#' @export
#' @references
#' Akritas, M.G., Murphy, S.A., LaValley, M.P., 1995. The Theil-Sen Estimator With Doubly Censored Data and Applications to Astronomy. Journal of the American Statistical Association 90, 170–177. https://doi.org/10.2307/2291140
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1 edition. ed. John Wiley and Sons, USA, N.J.
#' @seealso [NADA::cenken]
#' @examples
#' library(NADA) #For example data
#'
#' # x is not censored
#' data(HgFish)
#'
#' with(HgFish, ATSmini(Hg, HgCen, SedLOI))

ATSmini <- function(y.var, y.cen, x.var) {
  y.cen <- as.logical(y.cen)
  nobs <- length(y.var)

  Noshift<-cenken(y.var, y.cen, x.var)
  slope <- Noshift$slope
  intercept <- Noshift$intercept
  pval <- Noshift$p
  tau <- Noshift$tau
  S <- tau*nobs*(nobs-1)*0.5

  if (min(y.var < 0)) {
    # y goes negative. Intercept is NA.
    shift.amt <- abs(min(y.var)) + 0.0001
    y.shift <- y.var + shift.amt
    Shifted<-cenken(y.shift, y.cen, x.var)
    intercept <- Shifted$intercept - shift.amt
  }

  out <- data.frame(intercept, slope, tau, pval, S)
  return(out)
}

#' Computes Kendall rank correlation S-value (score) for censored data
#'
#' @description A permutation test to compute Kendall rank correlation S-values on censored data.
#' @param x Column of the time variable, either a sequence of days or decimal times, etc.  Will be the scale used for time in the trend analysis.
#' @param y The column of y (response variable) values plus detection limits
#' @param ycen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y` column, and 0 (or `FALSE`) indicates a detected value in `y`.
#' @param seas Column of the season classifications. A factor in R, so usually though not necessarily a text variable.  If numeric, define as a factor before running the script.
#' @param R The number of repetitions in the permutation process.  R is often between 999 and 9999 (+ the 1 observed test statistic produces 1000 to 10000 repetitions).
#'
#' @return S-value for each permutation
#' @export
#' @seealso [Kendall::Kendall]
#' @examples
#' data(Brumbaugh)
#'
#' #Artifical time and season variables for demonstration purposes
#' Brumbaugh$time=1:nrow(Brumbaugh)
#' Brumbaugh$sea=as.factor(round(runif(nrow(Brumbaugh),1,4),0))
#'
#'
#' with(Brumbaugh,computeS(time,Hg,HgCen,sea,R=100))
#'
computeS<- function(x, y, ycen, seas = NULL, R=R) {
  #  called by censeaken.  Returns permutations of S for a single season
  nx = length(x); nj <- nx-1
  Sout <- matrix(0, nrow =R, ncol=1)
  delmax <- 0
  delmin <-0
  delsign<-0
  for (i in 1:R) {xperm <- sample(x)
  yperm <- y[order(xperm)]
  ycenperm <- ycen[order(xperm)]
  Ymin <- yperm*(1-ycenperm)   # 0 for NDs
  Ymax <- Ymin+(ycenperm*0.99*yperm)  # <1 becomes 0.99
  Sindiv <- rep(0, nx*(nj)/2)
  ns <- 1
  for (j in 1:nj) {
    for (k in (j+1):nx) {
      delmax <- sign(Ymax[k] - Ymax[j])
      delmin <- sign(Ymin[k] - Ymin[j])
      delsign <- abs(sign(delmax - delmin)) # if 0, delmax is correct.  if 1 Sindiv =0
      Sindiv[ns] <- (1-delsign)*delmax
      ns <- ns+1
    }
  }
  Sout[i,1] <- sum(Sindiv)
  }
  # print(Sout)
  return (Sout)
}

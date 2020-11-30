#' Kendall's S-statistic for permutations of censored data
#'
#' @description Computes a Kendall rank correlation S-statistic for permutations of censored data. Collectively these represent the variation in S expected when the null hypothesis is true.  Called by censeaken. computeS is not expected to be of much use to users on its own.
#' @param x Column of the time variable, either a sequence of days or decimal times, etc.  Time data for one season.
#' @param y The column of y (response variable) values plus detection limits for one season.
#' @param ycen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y` column, and 0 (or `FALSE`) indicates a detected value in `y`.
#' @param seas Name of a single season classification. Usually though not necessarily a text variable.
#' @param R The number of repetitions in the permutation process.  R is often between 999 and 9999 (+ the 1 observed test statistic produces 1000 to 10000 realizations).
#'
#' @return An Rx1 matrix containing an S-value for each of the R data permutations.
#' @export
#' @seealso [Kendall::Kendall]
#' @references
#' Helsel, D.R., Hirsch, R.M., Ryberg, K.R., Archfield, S.A., Gilroy, E.J., 2020. Statistical Methods in Water Resources. U.S. Geological Survey Techniques and Methods, book 4, chapter A3, 458p., https://doi.org/10.3133/tm4a3.

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

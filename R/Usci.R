#' Interval-censored U-Score
#'
#' @description Interval-censored computation of uscores and their ranks for 1 parameter.  Called by uscoresi. Usci is not expected to be of much use to users on its own.
#' @param ylo The lower end of the concentration interval
#' @param yhi The upper end of the concentration interval
#' @param rnk A `TRUE`/`FALSE` variable on whether to compute the multivariate pattern on the uscores, or the ranks of the uscores.  Default is rnk=`TRUE`, use the ranks. rnk = `FALSE` returns the uscores.
#' @return Returns a single column of uscores or the ranks of uscores for a single pair of (low, high) interval-censored data columns.
#'
#' @export
#'
#' @examples
#'
#' data(Brumbaugh)
#'
#' # for demonstration purposes create a lower end concentration interval
#' Brumbaugh$lowHg<-Brumbaugh$Hg*(1-Brumbaugh$HgCen)
#'
#' with(Brumbaugh,Usci(lowHg,Hg))


Usci <- function(ylo, yhi, rnk=TRUE){
  x <- na.omit(data.frame (ylo, yhi))
  n = length(x$ylo)
  yadj=x$yhi-(sign(x$yhi-x$ylo)*0.001*x$yhi) #sets a <1 to be <1
  overlap=x$yhi  # sets correct dimensions
  Score=overlap  # sets correct dimensions
  for (j in 1:n) {
    for (i in 1:n ){
      overlap[i]=sign(sign(yadj[i]-x$ylo[j])+sign(x$ylo[i]-yadj[j]))
    }
    Score[j] = -1*sum(overlap) # -1 so that low values = low scores
  }

  if (rnk) {Uscore=rank(Score)} else {Uscore = Score}
  # print(Score)
  # print(Uscore)
  return(Uscore)
}

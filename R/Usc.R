#' Interval-censored U-Score
#'
#' @description Interval-censored computation of uscores and their ranks for 1 parameter.
#' @param ylo the lower end of the concentration interval
#' @param yhi the upper end of the concentration interval
#' @param rnk A `TRUE`/`FALSE` variable on whether to compute the multivariate pattern on the uscores, or the ranks of the uscores.  Default is rnk=`TRUE`, use the ranks. rnk = `FALSE` returns the uscores.
#'
#' @export


Usc <- function(ylo, yhi, rnk=T){
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

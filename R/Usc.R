#' U-scores for (non-interval, sinle-column) Censored Data
#'
#' @description Computes the column of uscores from 2 columns of data in the indicator value format. Multiple detection limits allowed.  Called by the uscores function, Usc (this function) is not expected to be of much use to users on its own.
#' @param y The column of data values plus detection limits
#' @param ind The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y` column, and 0 (or `FALSE`) indicates a detected value in `y`.
#' @param rnk  A `TRUE`/`FALSE` variable on whether to compute the multivariate pattern on the uscores, or the ranks of the uscores.  Default is rnk=`TRUE`, use the ranks. rnk = `FALSE` returns the uscores.
#' @export
#' @return Returns a single column of uscores or the ranks of uscores for a single pair of (concentration, indicator) censored data columns.
#'
#' @examples
#' data(Brumbaugh)
#' uscore(Brumbaugh$Hg,Brumbaugh$HgCen)

Usc <- function(y, ind, rnk=TRUE){
  x <- na.omit(data.frame (y, ind))
  n=length(x$y)
  ylo=(1-as.integer(x$ind))*x$y
  # yadj=y
  yadj=x$y-(sign(x$y-ylo)*0.001*x$y)
  overlap=x$y
  Score=overlap
  for (j in 1:n) {
    for (i in 1:n ){
      overlap[i]=sign(sign(yadj[i]-ylo[j])+sign(ylo[i]-yadj[j]))
    }
    Score[j] = -1*sum(overlap) # -1 so that low values = low scores
  }

  if (rnk) {uscore=rank(Score)} else {uscore = Score}
  return(uscore)
}

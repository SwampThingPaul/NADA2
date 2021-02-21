
#' U-scores for interval-censored data (multiple columns)
#'
#' @description Computes uscores or the ranks off uscores within columns of interval-censored data (the "i").  Data may have one or more detection limits.
#' @param dat.frame A data frame. Default format is: paired = `TRUE`, the default input format is ylo1 yhi1 ylo2 yhi2 ylo3 yhi3, where ylo is the low end of the interval of possible values, and yhi is the high end. There is a pair of columns for each censored parameter.
#' @param paired An option to specify paired = `FALSE`, where the format would be ylo1 ylo2 ylo3 yhi1 yhi2 yhi3, low values for each parameter followed by the high values in the same order.
#' @param rnk rnk=`TRUE` returns the ranks of uscores.  rnk = `FALSE` returns the uscores themselves. Default is rnk = `TRUE` to return the ranks.
#' @param Cnames Cnames =`1` uses the "lo" column names to name the uscores columns (the default).  Cnames = `2` uses the "hi" column names.
#' @details Input is a data.frame of paired low and high possible range of values, in an interval-censored format. ylo = the lower end of the interval is the first (left) column in the pair.  yhi is the upper end of the interval, at the second (right) column in the pair. For a detected value, ylo=yhi.  For a ND,  ylo != yhi. The uscore is the number of observations known to be lower minus the number of observations known to be higher. Ties, including those such as <1 vs <3 or 4 vs 4 or <3 vs 2, are given a 0 value in the uscore computation. The ranks of uscores provides a scale that is often more manageable than the uscores themselves.
#' @return prints the uscore number of observations known to be lower - number of observations known to be higher, for each observation.
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @export
#'


uscoresi <- function(dat.frame, paired = TRUE, rnk=TRUE, Cnames = 1) {
  cols <- ncol(dat.frame)
  if (cols > 2) {half <- cols/2        # multiple pairs of columns
  dat.orig <- as.matrix(na.omit (dat.frame))
  j=0; data.0 <- dat.orig
  if (paired) { for (i in seq(1, to=(cols-1), by=2)) {j = j+1
  data.0[,j] <- dat.orig[,i]
  data.0[,j+half] <- dat.orig[,i+1]   # reorganizing paired columns.
  nvec <- colnames(data.0)
  if (Cnames == 1) {nvec <- nvec[seq(1, (cols-1), 2)] }
  else {nvec <- nvec[seq(2, cols, 2)] }
  nvec <- paste("usc.", nvec, sep="")
  } }
  else { data.0 <- dat.orig     # no reorg needed,  columns in blocks.
  nvec <- colnames(data.0)
  nvec <- nvec[1:half]
  nvec <- paste("usc.", nvec, sep="")
  }

  u.out <- data.0[,1:half]
  for (i in 1:half) {u <- Usci(data.0[,i], data.0[,i+half], rnk=rnk)
  if (i==1) {u.out <- u}
  else {u.out <- cbind(u.out, u)}
  }
  colnames(u.out) <- nvec
  return(u.out)
  }
  else {                   # only 1 pair of columns
    x <- na.omit(data.frame (dat.frame))
    names(x) <- c("ylo", "yhi")
    n = length(x$ylo)
    print (n)
    yadj=x$yhi-(sign(x$yhi-x$ylo)*0.001*x$yhi)   #  sets a <1 to be <1
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
}

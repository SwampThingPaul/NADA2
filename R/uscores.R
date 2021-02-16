#' Uscores for multiple columns of censored data
#'
#' @description Computes uscores or the ranks of uscores of censored data in the indicator format. Multiple DLs allowed.
#' @param dat.frame A data frame. Default format is paired = `TRUE`, where for 3 chemical parameters the input format is C1 I1 C2 I2 C3 I3, a concentration column followed by its corresponding indicator column.
#' @param paired When paired = `FALSE`, the input format is C1 C2 C3 I1 I2 I3 where the C columns contain concentrations or a detection limit, and the I columns are their associated indicators, in the same order as the concentration columns.
#' @param rnk A logical `TRUE`/`FALSE` variable on whether to compute the multivariate pattern on the uscores, or the ranks of the uscores.  Default is rnk=`TRUE`, use the ranks. rnk = `FALSE` returns the uscores.
#'
#' @return A matrix of uscores or ranks of uscores, one column for each chemical parameter.  If there is only one chemical parameter a vector of uscores or ranks of uscores is returned.
#' @export
#'
#' @examples
#' data(PbHeron)
#'
#' uscores(PbHeron[,4:15])

uscores <- function(dat.frame, paired = TRUE, rnk=TRUE) {
  cols <- ncol(dat.frame)
  if (cols >2) { half <- cols/2     # multiple pairs of columns
  dat.orig <- as.matrix(na.omit (dat.frame))
  j=0;  uscr=0; DLmax=0; data.1 <- dat.orig; data.0 <- dat.orig
  if (paired) { for (i in seq(1, to=(cols-1), by=2)) {j = j+1
  data.0[,j] <- dat.orig[,i]
  data.0[,j+half] <- dat.orig[,i+1]   # reorganizing paired columns
  nvec <- colnames(data.0)
  nvec <- nvec[seq(1, (cols-1), 2)]
  nvec <- paste("usc.", nvec, sep="")
  } }
  else {data.0 <- dat.orig   # no reorg needed.  columns in blocks.
  nvec <- colnames(data.0)
  nvec <- nvec[1:half]
  nvec <- paste("usc.", nvec, sep="")
  }

  u.out <- data.0[,1:half]
  for (i in 1:half) {u <- Usc(data.0[,i], data.0[,i+half], rnk=rnk)
  if (i==1) {u.out <- u}
  else {u.out <- cbind(u.out, u)}
  }
  colnames(u.out) <- nvec
  return(u.out)
  }

  else {    # only 1 pair of columns
    x <- na.omit(dat.frame)
    names(x) <- c("y", "ind")
    n=length(x$y)
    ylo=(1-as.integer(x$ind))*x$y
    yadj=x$y-(sign(x$y-ylo)*0.001*x$y)
    overlap=x$y  # sets correct dimensions
    Score=overlap  # sets correct dimensions
    for (j in 1:n) {
      for (i in 1:n ){
        overlap[i]=sign(sign(yadj[i]-ylo[j])+sign(ylo[i]-yadj[j]))
      }
      Score[j] = -1*sum(overlap)     # -1 so that low values = low scores
    }

    if (rnk) {uscore=rank(Score)} else {uscore = Score}
    # print(Score)
    # print(Urank)
    # print(uscore)
    return(uscore)
  }
}

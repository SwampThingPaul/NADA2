#' Computes ranks of data with one or more detection limits
#'
#' @description Computes ranks of data (within columns) with one or more DLs, re-censoring at the highest DL if multiple DLs in column are present.
#' @param dat.frame A data frame. Default format is paired = `TRUE`, where for 3 chemical parameters the input format is C1 I1 C2 I2 C3 I3, a concentration column followed by its corresponding indicator column.
#' @param paired An option to specify paired = `FALSE`, where the input format would be C1 C2 C3 I1 I2 I3 where the C columns contain concentrations or a detection limit, and the I columns are their associated indicators, in the same order as the concentration columns.
#' @export
#' @return Returns a rank for each chemical parameter.
#'
#' @examples
#' library(NADA) #For example data
#' data(Golden)
#'
#' ordranks(Golden[,4:15])
#'

ordranks <- function(dat.frame, paired = TRUE) {
  cols <- ncol(dat.frame)
  half <- cols/2
  dat.orig <- as.matrix(dat.frame)
  j=0;  ranks=0; DLmax=0; data.1 <- dat.orig; data.0 <- dat.orig
  if (paired) { for (i in seq(1, to=(cols-1), by=2)) {j = j+1
  data.0[,j] <- dat.orig[,i]
  data.0[,j+half] <- dat.orig[,i+1]
  nvec <- colnames(data.0)
  nvec <- nvec[seq(1, (cols-1), 2)]
  nvec <- paste("rnk.", nvec, sep="")
  } }
  else {data.0 <- dat.orig
  nvec <- colnames(data.0)
  nvec <- nvec[1:half]
  nvec <- paste("rnk.", nvec, sep="")
  }

  for (i in 1:half) {DLmax <- max(data.0[,i]*data.0[,i+half])
  data.1[,i] <- data.0[,i]*(1-data.0[,i+half])  # all <ND = 0 in data.1
  data.1[,i][data.1[,i] < DLmax] = 0   # all detects below DLmax = 0
  }
  data.2 <- data.1[,1:half]
  r.out <- data.2
  for (i in 1:half) {ranks <- rank(data.2[,i])
  if (i==1) {r.out <- ranks}
  else {r.out <- cbind(r.out, ranks)}
  }
  colnames(r.out) <- nvec
  return(r.out)
}

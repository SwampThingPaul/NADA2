#' Binary similarity coefficient matrix
#'
#' @description Computes a simple matching similarity coefficient
#' @param dat.frame A data frame containing only the 0/1 columns.
#' @export
#' @importFrom vegan designdist
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @return Returns a binary similarity matrix.
#' @seealso [vegan::designdist]
#' @examples
#' data(PbHeron)
#'
#' binarySim(PbHeron$LiverCen)
#'

binarySim <- function(dat.frame) {
  dat.symm <- designdist(dat.frame, method = "(a+d)/(a+b+c+d)", abcd = TRUE, terms = "binary", name = "simplematch")
  #  This returns a similarity matrix.  0 = disjoint.
  #  for a dissimilarity matrix, use binaryDiss
  return(dat.symm)
}

#' Binary similarity coefficient
#'
#' @description Computes a simple matching similarity coefficient
#' @param dat.frame A data frame containing only the 0/1 columns.
#' @export
#' @importFrom vegan designdist
#' @return Returns a binary similarity matrix.
#' @seealso [vegan::designdist]
#' @examples
#' library(NADA) #For example data
#' data(Golden)
#'
#' binarySim(Golden$LiverCen)
#'

binarySim <- function(dat.frame) {
  dat.symm <- designdist(dat.frame, method = "(a+d)/(a+b+c+d)", abcd = TRUE, terms = "binary", name = "simplematch")
  #  This returns a similarity matrix.  0 = disjoint.
  #  for a dissimilarity matrix, use binaryDiss
  return(dat.symm)
}

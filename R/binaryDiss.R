#' Binary dissimilarity coefficient matrix
#'
#' @description Computes a simple matching dissimilarity coefficient
#' @param dat.frame A data frame containing only the 0/1 columns.
#' @export
#' @importFrom vegan designdist
#' @return Returns a binary dissimilarity matrix.
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @seealso [vegan::designdist]
#' @examples
#' data(PbHeron)
#'
#' binaryDiss(PbHeron$LiverCen)

binaryDiss <- function(dat.frame) {
  dat.diss <- designdist(dat.frame, method = "1 - (a+d)/(a+b+c+d)", abcd = TRUE, terms = "binary", name = "simplematch")
  #  This returns a distance matrix.  0 = identical.
  #  for a similarity matrix, use binarySim
  return(dat.diss)
}

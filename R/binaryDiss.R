#' Binary dissimilarity coefficient
#'
#' @description Computes a simple matching dissimilarity coefficient
#' @param dat.frame A data frame containing only the 0/1 columns.
#' @export
#' @importFrom vegan designdist
#' @return Returns a binary dissimilarity matrix.
#' @seealso [vegan::designdist]
#' @examples
#' library(NADA) #For example data
#' data(Golden)
#'
#' binaryDiss(Golden$LiverCen)

binaryDiss <- function(dat.frame) {
  dat.diss <- designdist(dat.frame, method = "1 - (a+d)/(a+b+c+d)", abcd = TRUE, terms = "binary", name = "simplematch")
  #  This returns a distance matrix.  0 = identical.
  #  for a similarity matrix, use binarySymm
  return(dat.diss)
}

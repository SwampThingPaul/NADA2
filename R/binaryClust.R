#' Binary Cluster Matrix
#'
#' @description Performs clustering of a matrix of 0s and 1s, ie. the censoring indicator columns for multiple variables.
#' @param dat.frame A data frame containing only the 0/1 columns.
#' @param method Method of forming clusters.  The default is `"ward.D2"`, which is appropriate for a variety of types of data.  Another appropriate option is `“average”` – average distances between cluster centers.
#' @param group Optional grouping variable. Sites being cluster will be represented by their group name, rather than by the row number.
#' @param ncluster Optional number of clusters to be differentiated on the graph. Clusters are fenced off with rectagles.
#' @importFrom stats hclust cutree rect.hclust
#' @export
#' @return Prints a cluster dendrogram based on clustering method
#'
#' @examples
#' data(PbHeron)
#'
#' # without group specified
#' binaryClust(PbHeron[,4:15])
#'
#' # With Group argument
#' binaryClust(PbHeron[,4:15],group=PbHeron$DosageGroup)

binaryClust <- function(dat.frame, method = "ward.D2", group = NULL, ncluster = NULL) {

  simple_matching_coeff <- binaryDiss(dat.frame)
  dat.clust = hclust(simple_matching_coeff, method = method)
  if (is.null(group))  {plot(dat.clust)}
  else { plot(dat.clust, label = group)}

  if (!is.null(ncluster))   {c5=cutree(dat.clust, k=ncluster)
  rect.hclust(dat.clust, ncluster)
  }
}

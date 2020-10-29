#' Cluster Matrix of Binary Censored Data
#'
#' @description Performs clustering of a matrix of 0s and 1s, ie. the censoring indicator columns for multiple variables. If there are multiple detection limits within a column first convert the 0/1 designated to be Below vs Above the highest detection limit in the column.  Detection limits may differ among columns.
#' @param dat.frame A data frame containing only the 0/1 columns, one column per chemical parameter.
#' @param method Method of forming clusters.  The default is `"ward.D2"`, which is appropriate for a variety of types of data.  Another appropriate option is `“average”` – average distances between cluster centers.  See the vegan package for other possible clustering methods.
#' @param group Optional grouping variable. If used, sites being clustered will be represented by their group name, rather than by the row number.
#' @param ncluster Optional number of clusters to be differentiated on the graph. Clusters are fenced off with rectangles.
#' @importFrom stats hclust cutree rect.hclust
#' 
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.

#'@export
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

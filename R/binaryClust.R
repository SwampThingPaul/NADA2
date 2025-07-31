#' Cluster Matrix of Binary Censored Data
#'
#' @description Performs clustering of a matrix of 0s and 1s, ie. the censoring indicator columns for multiple variables. If there are multiple detection limits within a column first convert the 0/1 designated to be Below vs Above the highest detection limit in the column.  Detection limits may differ among columns.
#' @param data A data frame containing only the 0/1 columns, one column per chemical parameter.
#' @param method Method of forming clusters.  The default is `"ward.D2"`, which is appropriate for a variety of types of data.  Another appropriate option is `“average”` – average distances between cluster centers.  See the vegan package for other possible clustering methods.
#' @param group Optional grouping variable. If used, sites being clustered will be represented by their group name, rather than by the row number.
#' @param ncluster Optional number of clusters to be differentiated on the graph. Clusters are fenced off with rectangles.
#' @param plotncluster default is `TRUE` logical flags to add identification of clusters on dendrogram
#' @param clustIndex Optional, if not specified, potential number of clusters will be determined based on the mean best number of clusters across all indicies. For a specific index, see details
#' @importFrom stats hclust cutree rect.hclust
#' @importFrom graphics plot
#' @importFrom NbClust NbClust
#'
#' @details
#'
#' If a specific index is desired to determine the best number of clusters see `NbClust::NbClust` for index values.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.

#'@export
#' @return Prints a cluster dendrogram based on clustering method and outputs a list of clusters and hierarchical cluster analysis results
#'
#' @examples
#' data(PbHeron)
#'
#' # without group specified
#' binaryClust(PbHeron[,4:15])
#'
#' # With Group argument
#' binaryClust(PbHeron[,4:15],group=PbHeron$DosageGroup)

binaryClust <- function(data, method = "ward.D2", group = NULL, ncluster = NULL,plotncluster = TRUE,clustIndex = "all") {
  simple_matching_coeff <- binaryDiss(data)

  if (!inherits(simple_matching_coeff, "dist")) {
    stop("binaryDiss must return a 'dist' object")
  }

  dat.clust <- hclust(simple_matching_coeff, method = method)

  if (is.null(group)) {
    graphics::plot(dat.clust)
  } else {
    graphics::plot(dat.clust, labels = group)
  }

  if (is.null(ncluster)) {
    # to suppress plotting NbClust
    grDevices::pdf(file = NULL)
    invisible(capture.output({
    result <- NbClust(data, distance = "euclidean", min.nc = 2, max.nc = 10, method = method,index = clustIndex)
    }))
    grDevices::dev.off()

    ncluster <- round(mean(result$Best.nc[1,]))
    c5 <- cutree(dat.clust, k = ncluster)
    rect.hclust(dat.clust, k = ncluster)
    return(list(hclust = dat.clust,
                clusters = c5))

  }else{c5 <- cutree(dat.clust, k = ncluster)
  rect.hclust(dat.clust, k = ncluster)
  return(list(hclust = dat.clust,
              clusters = c5))
  }

}

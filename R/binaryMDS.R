#' Plot Nonmetric Multidimensional Scaling of censored data
#'
#' @description Plots an NMDS of a matrix of 0s and 1s, the censoring indicator columns for multiple variables.  Use the highest censoring limit within each column.  May have different censoring levels in different columns.
#' @param dat.frame A data frame containing only the 0/1 columns.
#' @param group Optional grouping variable. Sites will be represented by different colored symbols for each group.
#' @param title Optional title for the NMDS graph.
#' @param legend.pos For when group is specified, the location of the legend on the graph showing the colors representing each group’s data.  Default is “bottomleft”.  Alternatives are “topright” and “centerleft”, etc.
#' @export
#' @importFrom vegan metaMDS
#' @return Print NMDS plot of censored data groupings.
#' @details Prior to runnning this analysis it is suggested to consult best analysis practise when performing NMDS. As a rule of thumb, an NMDS ordination with a stress value around or above 0.2 is deemed suspect and a stress value approaching 0.3 indicates that the ordination is arbitrary. Stress values equal to or below 0.1 are considered fair, while values equal to or below 0.05 indicate good fit.
#' @seealso [vegan::metaMDS]
#'
#' @examples
#' \dontrun{
#' library(NADA) #For example data
#' data(Golden)
#'
#' # without group specified
#' binaryMDS(Golden[,4:15])
#'
#' # With Group argument
#' binaryMDS(Golden[,4:15],group=Golden$DosageGroup)
#' }



binaryMDS <- function(dat.frame, group = NULL, title=NULL, legend.pos = "bottomleft") {
  matname <- deparse(substitute(dat.frame))
  if (is.null(title)) {title <- paste("NMDS of", matname)}
  rownumb <- nrow(dat.frame)
  row.labels <- as.character(1:rownumb)
  dat.dist <- binaryDiss(dat.frame)  # distance matrix.  0 is identical.
  dat.mds <- metaMDS(dat.dist, zerodist = "add", autotransform = FALSE)
  plot(dat.mds, type = "n", main = title)
  if (is.null(group)) {text(dat.mds$points,labels=row.labels)}
  else {
    gp.factor <- as.factor(group)
    gp.nlevels <- c(1:nlevels(gp.factor))
    gp.col <- as.integer(gp.factor)
    gp.symb <- 18+gp.col
    plot.col=ordiplot(dat.mds, type="none", display="sites", main=title)
    #    print(gp.symb)
    #    points(plot.col, "sites", pch=gp.symb, col=gp.col)
    points(plot.col, "sites", pch=19, col=gp.col)
    text(plot.col, "sites", pos=4, cex = 0.8)
    #    legend(legend.pos, legend=levels(gp.factor), bty="n",col = gp.nlevels, pch = gp.symb)
    legend(legend.pos, legend=levels(gp.factor), bty="n",col = gp.nlevels, pch = 19)
  }

}

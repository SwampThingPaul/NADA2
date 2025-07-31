#' Plot Nonmetric Multidimensional Scaling of binary censored data
#'
#' @description Plots an NMDS of a matrix of 0s and 1s, the censoring indicator columns for multiple variables, to discern the pattern of data below vs. above the detection limit.  With multiple detection limits within a column, re-censoring to the highest limit in the column must be done prior to running this function.  May have different censoring levels in different columns.
#' @param dat.frame A data frame containing only the columns of 0/1 values.
#' @param group Optional grouping variable. Sites will be represented by different colored symbols for each group.
#' @param title Optional title for the NMDS graph.
#' @param legend.pos When group is specified, determines the location of the legend on the graph showing the colors representing each group’s data.  Default is “bottomleft”.  Alternatives are “topright” and “centerleft”, etc.
#' @importFrom vegan metaMDS ordiplot
#' @importFrom graphics points
#' @export
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @return Plots an NMDS of censored data represented as the binary Above vs Below a detection limit for each parameter.
#' @details Binary data may not provide sufficient information to discern differences in location on the plot if sample size is small.  Prior to running this analysis it is suggested to consult best analysis practice when performing NMDS. As a rule of thumb, an NMDS ordination with a stress value around or above 0.2 is deemed suspect and a stress value approaching 0.3 indicates that the ordination is arbitrary. Stress values equal to or below 0.1 are considered fair, while values equal to or below 0.05 indicate good fit.
#' @seealso [vegan::metaMDS]
#'
#' @examples
#' \dontrun{
#' data(PbHeron)
#'
#' # without group specified
#' binaryMDS(PbHeron[,4:15])
#'
#' # With Group argument
#' binaryMDS(PbHeron[,4:15],group=PbHeron$DosageGroup)
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

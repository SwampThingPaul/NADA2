
#' Plot U-score Nonmetric Multidimensional Scaling of censored data
#'
#' @description Plots an NMDS of uscores output from the `uscores` or `uscoresi` functions.
#' @param uscor A data frame of uscores or ranks of uscores produced by either the `uscores(...)` or `uscoresi(...)` functions
#' @param group Optional grouping variable. Sites will be represented by different colored symbols for each group.
#' @param title Optional title for the NMDS graph.
#' @param legend.pos For when group is specified, the location of the legend on the graph showing the colors representing each group’s data.  Default is “bottomleft”.  Alternatives are “topright” and “centerleft”, etc.
#' @importFrom vegan metaMDS
#' @importFrom stats dist
#' @return Prints an NMDS plot of censored data groupings based on U-scores
#' #' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @export
#'
#' @examples
#' data(PbHeron)
#'
#' PbHeron.u <- uscores(PbHeron[,4:15])
#' uMDS(PbHeron.u)
#'
#' # With group specific
#' uMDS(PbHeron.u,group=PbHeron$DosageGroup)

uMDS <- function(uscor, group = NULL, title=NULL, legend.pos = "bottomleft") {
  uname <- deparse(substitute(uscor))
  if (is.null(title)) {title <- paste("NMDS of", uname)}
  rownumb <- nrow(uscor)
  row.labels <- as.character(1:rownumb)
  uclid <- dist(uscor)  # euclidean distance matrix of uscores.
  u.mds <- metaMDS(uclid, zerodist = "add")
  if (is.null(group)) {plot(u.mds, type = "p", main = title)
    text(u.mds$points, pos=4, cex = 0.9, offset = 0.25)
  }
  else {
    gp.factor <- as.factor(group)
    gp.nlevels <- c(1:nlevels(gp.factor))
    gp.col <- as.integer(gp.factor)
    gp.symb <- 18+gp.col
    plot.col=ordiplot(u.mds, type="none", display="sites", pch = 19, main=title)
    points(plot.col, "sites", pch=19, col=gp.col)
    text(plot.col, "sites", pos=4, cex = 0.8)
    legend(legend.pos, legend=levels(gp.factor), bty="n",col = gp.nlevels, pch = 19)
  }
  return(invisible(u.mds))
}

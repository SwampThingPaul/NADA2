#' Permutation analysis of similarities between groups
#'
#' @param ano.out an `anosim` object. See details and example
#' @param hcol color of histogram
#' @param title title of histogram
#' @return prints a histogram of permutation tests
#' @export
#' @importFrom vegan anosim
#' @importFrom graphics hist
#' @seealso [vegan::anosim]
#' @references
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' Oksanen, J., Guillaume, F., 2018. Vegan: ecological diversity. CRAN R-Project.<https://cran.r-project.org/web/packages/vegan/index.html>
#'
#' @examples
#' data(PbHeron)
#'
#' # ROS model for each group
#' PbHeron.high <- with(subset(PbHeron,DosageGroup=="High"),NADA::ros(Blood,BloodCen))
#' PbHeron.high <- data.frame(PbHeron.high)
#' PbHeron.high$DosageGroup <- "High"
#'
#' PbHeron.low <- with(subset(PbHeron,DosageGroup=="Low"),NADA::ros(Blood,BloodCen))
#' PbHeron.low <- data.frame(PbHeron.low)
#' PbHeron.low$DosageGroup <- "Low"
#'
#' PbHeron.ros=rbind(PbHeron.high,PbHeron.low)
#'
#' # ANOSIM analysis
#' library(vegan)
#' PbHeron.anosim <- with(PbHeron.ros,anosim(modeled,DosageGroup))
#' summary(PbHeron.anosim)
#'
#' # Plot
#' anosimPlot(PbHeron.anosim)


anosimPlot <- function(ano.out, hcol = "light blue", title = "Histogram of anosim permutations") {
  min.p <- min(ano.out$perm)
  hset <- hist(ano.out$perm, plot = FALSE)
  yset <- max(hset$counts)
  label.txt <- paste("R =", signif(ano.out$statistic,2),"")
  xmax <- max(ano.out$perm, ano.out$statistic)

  hist(ano.out$perm, col = hcol, xlim = c(min.p, xmax), main = title, xlab = "Test statistics")
  abline (v = ano.out$statistic, lwd = 2)
  text (ano.out$statistic, yset, labels = label.txt, adj = c(1,1))
}

#' Permutation analysis of similarities between groups
#'
#' @param ano.out an `anosim` object. See details and example
#' @param hcol color of histogram
#' @param title title of histogram
#' @return prints a histogram of permutation tests
#' @export
#' @importFrom vegan anosim
#' @seealso [vegan::anosim]
#' @references
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' Oksanen, J., Guillaume, F., 2018. Vegan: ecological diversity. CRAN R-Project.<https://cran.r-project.org/web/packages/vegan/index.html>
#'
#' @examples
#' library(NADA) #For example data
#' data(Golden)
#'
#' # ROS model for each group
#' golden.high <- with(subset(Golden,DosageGroup=="High"),NADA::ros(Blood,BloodCen))
#' golden.high <- data.frame(golden.high)
#' golden.high$DosageGroup <- "High"
#'
#' golden.low <- with(subset(Golden,DosageGroup=="Low"),NADA::ros(Blood,BloodCen))
#' golden.low <- data.frame(golden.low)
#' golden.low$DosageGroup <- "Low"
#'
#' golden.ros=rbind(golden.high,golden.low)
#'
#' # ANOSIM analysis
#' library(vegan)
#' golden.anosim <- with(golden.ros,anosim(modeled,DosageGroup))
#' summary(golden.anosim)
#'
#' # Plot
#' anosimPlot(golden.anosim)


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

#' Peto-Peto nonparametric test
#'
#' @description
#' Performs a Peto-Peto nonparametric test of differences in cdfs between groups.  If more than two groups, the test is followed by a nonparametric multiple comparison test.  Uses the BH method of adjusting p-values.
#' @param y1 The column of data values plus detection limits
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the y1 column, and 0 (or `FALSE`) indicates a detected value in y1.
#' @param grp Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @keywords cdf
#' @export
#' @return  1-way test to test differences between groups
#'

#' @import survminer
#' @import survival
#'

#' @examples
#'
#' library(NADA)
#' data(Golden)
#' cen1way(Golden$Liver,Golden$LiverCen,Golden$Group)

cen1way <- function(y1,y2, grp) {
  yname <- deparse(substitute(y1))
  gname <- deparse(substitute(grp))
  rho=1
  fconst <- max(y1) + 1
  flip <- fconst - y1
  detect <- as.logical(1 - as.integer(y2))  # reverses TRUE/FALSE to fit survfit
  Factor <- as.factor(grp)
  df <- length(levels(Factor))-1
  CensData <- data.frame (flip, detect, Factor)
  Cen.stats <- matrix(0, nrow=df+1, ncol = 5)

  y.surv <- Surv(CensData$flip, CensData$detect, type="right")
  y.out<- survdiff(y.surv ~ CensData$Factor, rho=rho)
  pval = pchisq(y.out$chisq, df, lower.tail = FALSE)

  groupnames <- as.character(levels(Factor))
  for (i in 1:nlevels(Factor))    {
    y1gp <- y1[Factor==groupnames[i]]
    y2gp <- y2[Factor==groupnames[i]]
    Cstats <- suppressWarnings(cfit(y1gp, as.logical(y2gp), printstats=FALSE, Cdf = FALSE))
    Cstats <- Cstats[c(1:6)]
    Cstats <- Cstats[-3]
    Cstats <- data.frame(Cstats)
    if (i ==1) {Cen.stats <- Cstats
    cnames <- colnames(Cstats)}
    else {Cen.stats <- rbind(Cen.stats, Cstats)}
  }
  rownames(Cen.stats) <-  groupnames
  colnames(Cen.stats) <- cnames
  print.data.frame(Cen.stats, print.gap = 3)
  cat('\n',"     Oneway Peto-Peto test of CensData:", yname, "  by Factor:", gname, '\n', "     Chisq =", signif(y.out$chisq, 4), "  on", df, "degrees of freedom", "    p =", signif(pval,3), '\n')
  if (df >1) {mcomp <- pairwise_survdiff(Surv(flip, detect) ~ Factor, data=CensData, rho=rho)
  print(mcomp)}
}

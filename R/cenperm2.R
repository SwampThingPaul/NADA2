#' Censored two-group permutation test
#'
#' @description Performs a permutation test of differences in means between two groups of censored data.
#' @param y1 The column of data values plus detection limits
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y1` column, and 0 (or `FALSE`) indicates a detected value in `y1`.
#' @param grp Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param R The number of permutations used. Default is 9999
#' @param alternative indicates the alternative hypothesis and must be one of "`two.sided`", "`greater`" or "`less`". You may also specify just the initial letter. Default is "`two.sided`".
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @keywords permutation difference test
#' @export
#' @return Permutation test results with the number of permutations, range in group means and their difference, and range in `p-value`.
#' @details Because this is a permutation test it avoids the problem with MLE tests (`cen2means`) that assume a normal distribution.  No values are modeled as below zero and `p-values` are trustworthy. Ranges in means and p-values are due to interval-censoring of censored data means.
#'
#' @references
#' Good, P., 2000. Permutation Tests: A Practical Guide to Resampling Methods for Testing Hypotheses, 2nd ed, Springer Series in Statistics. Springer-Verlag, New York, NY. \doi{https://doi.org/10.1007/978-1-4757-3235-1}
#'
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#' data(PbHeron)
#' cenperm2(PbHeron$Liver,PbHeron$LiverCen,PbHeron$DosageGroup,alternative="t")

cenperm2 <- function(y1, y2, grp, R = 9999, alternative = "two.sided",printstat=TRUE) {
  yname <- deparse(substitute(y1))
  gname <- deparse(substitute(grp))
  xdat <- na.omit(data.frame(y1, y2, grp))
  x1 <- xdat[,1]
  x2 <- xdat[,2]
  Factor <- as.factor(xdat[,3])
  df <- length(levels(Factor))-1
  grpname <- levels(Factor)

  x1.lo <- x1*(1-x2)
  permdiff.lo <- vector(length = R)
  permdiff.hi <- vector(length = R)

  mu1hi <- mean(x1[Factor == grpname[1]])
  mu1lo <- mean(x1.lo[Factor == grpname[1]])
  mu2hi <- mean(x1[Factor == grpname[2]])
  mu2lo <- mean(x1.lo[Factor == grpname[2]])
  n1 <- length(x1[Factor == grpname[1]])
  n2 <- length(x1[Factor == grpname[2]])
  n <- n1+n2
  dbarhi <- mu1hi - mu2hi
  dbarlo <- mu1lo - mu2lo

  for (i in 1:R) {newFact <- sample(Factor)  # replace = FALSE
  permdiff.lo[i] <- mean(x1.lo[newFact == grpname[1]]) - mean(x1.lo[newFact == grpname[2]])
  permdiff.hi[i] <- mean(x1[newFact == grpname[1]]) - mean(x1[newFact == grpname[2]])
  }
  absdiff.hi = abs(dbarhi)
  absperm.hi = abs(permdiff.hi)
  if(alternative == 't'| alternative=="two.sided") pval.hi <- (1+sum(absdiff.hi<=absperm.hi))/(R+1)
  if(alternative == 'g'| alternative=="greater")  pval.hi <- (sum(permdiff.hi >= dbarhi)+1)/(R+1)
  if(alternative == 'l'| alternative=="less") pval.hi = (sum(permdiff.hi <= dbarhi)+1)/(R+1)

  absdiff.lo = abs(dbarlo)
  absperm.lo = abs(permdiff.lo)
  if(alternative == 't'| alternative=="two.sided") pval.lo <- (1+sum(absdiff.lo<=absperm.lo))/(R+1)
  if (alternative == 'g'| alternative=="greater")  pval.lo <- (sum(permdiff.lo >= dbarlo)+1)/(R+1)
  if (alternative == 'l'| alternative=="less") pval.lo = (sum(permdiff.lo <= dbarlo)+1)/(R+1)

  diff.lo <- min(dbarlo, dbarhi)
  diff.hi <- max(dbarlo, dbarhi)
  p.lo <- min(pval.lo, pval.hi)
  p.hi <- max(pval.lo, pval.hi)
  result <- data.frame(diff.lo, diff.hi, p.lo, p.hi)

  #added by PJ
  if(alternative%in%c('t','g','l')){alternative <- switch(alternative,"t"="two.sided","g"="greater","l"="less")}
  #  write results
  mean.txt <- paste("     Mean (", grpname[1], " - ", grpname[2], ") = ", sep="")

  if(printstat==TRUE){
  cat(" Permutation test of mean CensData:", yname, "  by Factor:", gname, '\n', "   ", R, "Permutations", "    alternative =", alternative, "\n")
  cat ("  mean of", grpname[1], "=", signif(mu1lo, 4), "to", signif(mu1hi,4), "    mean of", grpname[2], "=", signif(mu2lo,4), "to", signif(mu2hi,4), "\n")
  cat( mean.txt, signif(dbarlo, 4), "to", signif(dbarhi, 4), "      p =", pval.lo, "to", pval.hi, '\n', "\n")
  }
  return(invisible(result))
}

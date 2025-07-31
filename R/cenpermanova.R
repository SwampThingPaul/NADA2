#' Censored data one-factor permutation test
#'
#' @description Performs a permutation test of differences in means between groups of censored data.
#' @param y1 The column of data values plus detection limits.
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y1` column, and 0 (or `FALSE`) indicates a detected value in `y1`.
#' @param grp Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param R The number of permutations used. Default is 9999.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#'
#' @export
#' @return Permutation test results with the number of permutations, range in test statistics and `p-value` values through the various permutations. Group means are also listed.
#' @details Because this is a permutation test it avoids the problem with MLE tests (see `cenanova`) that assume a normal distribution.  No values are modeled as below zero and group means and `p-values` are trustworthy.
#'
#' @references
#' Good, P., 2000. Permutation Tests: A Practical Guide to Resampling Methods for Testing Hypotheses, 2nd ed, Springer Series in Statistics. Springer-Verlag, New York, NY. \doi{https://doi.org/10.1007/978-1-4757-3235-1}
#'
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @examples
#'
#' data(PbHeron)
#' cenpermanova(PbHeron$Liver,PbHeron$LiverCen,PbHeron$DosageGroup)
cenpermanova <- function(y1, y2, grp, R = 9999,printstat=TRUE) {
  yname <- deparse(substitute(y1))
  gname <- deparse(substitute(grp))
  xdat <- na.omit(data.frame(y1, y2, grp))
  x1 <- xdat[,1]
  x2 <- xdat[,2]
  Factor <- as.factor(xdat[,3])
  df <- length(levels(Factor))-1
  grpname <- levels(Factor)
  x1.lo <- x1*(1-x2)
  permW.lo <- vector(length = R)
  permW.hi <- vector(length = R)
  group.means <- vector(length = nlevels(Factor))
  mean.names <- vector(length = nlevels(Factor))
  gpname <- as.character(levels(Factor))

  # Test statistic Wprime is the between group sum of squares minus invariant terms.
  # see "Permutation Tests, Second Edition" by Good (2000), page 44.
  Wprime.groups <- function (x, Factor) {
    Wprime <- 0
    for (i in 1:nlevels(Factor)) {
      W.grp <- (sum(x[Factor == grpname[i]]))^2 / length(x[Factor == grpname[i]])
      Wprime <- Wprime + W.grp
    }
    return(invisible(Wprime))
  }

  # Test statistics on original data
  teststat.lo <- Wprime.groups(x1.lo, Factor)
  teststat.hi <- Wprime.groups(x1, Factor)
  absW.hi = abs(teststat.hi)
  absW.lo = abs(teststat.lo)
  # permutations
  for (i in 1:R) {newFact <- sample(Factor)  # replace = FALSE
  permW.lo[i] <- Wprime.groups(x1.lo, newFact)
  permW.hi[i] <- Wprime.groups(x1, newFact)
  }
  # p values always two-sided
  absperm.lo = abs(permW.lo)
  pval.lo <- (1+sum(absW.lo<=absperm.lo))/(R+1)
  absperm.hi = abs(permW.hi)
  pval.hi <- (1+sum(absW.hi<=absperm.hi))/(R+1)

  # group means
  for (i in 1:nlevels(Factor)) {
    ros.out <- suppressWarnings(cenros(x1[Factor == grpname[i]], as.logical(x2[Factor == grpname[i]])))
    group.means[i] <- signif(mean(ros.out$modeled),4)
    mean.names[i] <- paste("mean(", gpname[i], ")", sep="")
  }
  names(group.means) <- mean.names

  #  write test results
  p.lo <- min(pval.lo, pval.hi)
  p.hi <- max(pval.lo, pval.hi)
  result <- data.frame(teststat.lo, teststat.hi, p.lo, p.hi, group.means)
  if(printstat==TRUE){
    cat(" Permutation test of mean CensData:", yname, "  by Factor:", gname, '\n', "   ", R, "Permutations", "\n")
    cat( "Test Statistic =", signif(teststat.lo, 4), "to", signif(teststat.hi, 4), "      p =", pval.lo, "to", pval.hi, '\n', "\n")
    print(group.means, row.names = FALSE, print.gap = 3)
  }
  return(invisible(result))
}

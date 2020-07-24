#' Wilcoxcon Signed-Rank test for censored data
#'
#' @description Performs a nonparametric Wilcoxon signed-rank test of whether the median difference between two columns of paired censored data equals 0.  Uses the Pratt adjustment for pairs of equal values.
#' @param xd The first column of data values plus detection limits
#' @param xc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the xd column, and 0 (or `FALSE`) indicates a detected value in xd.
#' @param yd The second column of data values plus detection limits, or a single number representing a standard / guideline value.
#' @param yc The column of censoring indicators for yd, where 1 (or `TRUE`) indicates a detection limit in the yd column, and 0 (or `FALSE`) indicates a detected value in `yd`. Not needed if `yd` is a single standard number.
#' @param alternative The usual notation for the alternate hypothesis.  Default is `“two.sided”`.  Options are `“greater”` or `“less”`.
#'
#' @keywords Wilcoxon
#' @export
#' @return  A list of Wilcoxon Signed-Rank test with Pratt correction for ites statistics containing the following components:
#' \itemize{
#' \item `n` Number of samples
#' \item `Z` The value of the test statistic
#' \item `p.value` the p-value of the test
#' }
#' @references
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' Page, E.B., 1963. Ordered Hypotheses for Multiple Treatments: A Significance Test for Linear Ranks. Journal of the American Statistical Association 58, 216–230. <https://doi.org/10.2307/2282965>
#'
#' Pratt, J.W., 1959. Remarks on Zeros and Ties in the Wilcoxon Signed Rank Procedures. Journal of the American Statistical Association 54, 655–667. <https://doi.org/10.2307/2282543>
#'
#'
#' @examples
#'
#' library(NADA) #For example data
#' data(Atra)
#'
#' cen_signedranktest(Atra$June,Atra$JuneCen,Atra$Sept,Atra$SeptCen)

cen_signedranktest <- function(xd, xc, yd, yc, alternative="two.sided") {
  xname <- deparse(substitute(xd))
  yname <- deparse(substitute(yd))
  nonas <- na.omit(data.frame(xd, xc, yd, yc))
  if(length(nonas[,1]) != length(nonas[,3]))
    stop("Lengths of x and y must be the same for paired data.")
  if(alternative == "two.sided")
  {txt3 <- paste("alternative hypothesis: true difference", xname, "-", yname, "does not equal 0")}
  else if(alternative == "less")
  {txt3 <-  paste("alternative hypothesis: true difference", xname, "-", yname, "is less than 0")}
  else if (alternative == "greater")
  {txt3 <-  paste("alternative hypothesis: true difference", xname, "-", yname, "is greater than 0")}
  else {stop('Invalid choice for alternative, must match "two.sided", "greater", or "less"')}
  N = length(nonas[,1])
  for(i in seq(N)) {
    if(nonas[i,2] && nonas[i,4]) { # Both censored, make same
      xyMax <- max(nonas[i,1], nonas[i,3])
      nonas[i,1] <- nonas[i,3] <- xyMax
    }
    else if(nonas[i,2] && nonas[i,1] > nonas[i,3]) { # x censored > y
      nonas[i,3] <- nonas[i,1]
      nonas[i,4] <- TRUE
    }
    else if(nonas[i,4] && nonas[i,3] > nonas[i,1]) { # y censored > x
      nonas[i,1] <- nonas[i,3]
      nonas[i,2] <- TRUE
    }
  } # Done, no need to check uncensored observations
  # setting <s to slightly below detected obs at that value
  nonas[,1] <- nonas[,1] - nonas[,2] * (0.001) * nonas[,1]
  nonas[,3] <- nonas[,3] - nonas[,4] * (0.001) * nonas[,3]

  x <- nonas[,1]
  names(x) <- xname
  y <- nonas[,3]
  names(y) <- yname

  s.out <- wilcoxsign_test(x~y, alternative = alternative)
  txt <- paste("Censored signed-rank test for (x:", xname, " - ", "y:", yname, ") equals 0", "\n", txt3, "\n", sep = "")
  txt2 <- paste("n =", N, "  Z=", signif(statistic(s.out), 4), "  p-value =", signif(pvalue(s.out), 4))

  cat(txt, "\n","Pratt correction for ties", "\n", txt2, "\n")
}

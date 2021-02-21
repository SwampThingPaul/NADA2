#' Wilcoxcon Signed-Rank test for censored data
#'
#' @description Performs a nonparametric Wilcoxon signed-rank test of whether the median difference between two columns of paired censored data equals 0.  Uses the Pratt adjustment for pairs of equal or indistinguishable values.
#' @param xd The first column of data values plus detection limits
#' @param xc The column of censoring indicators for `xd`, where 1 (or `TRUE`) indicates a detection limit in the xd column, and 0 (or `FALSE`) indicates a detected value in xd.
#' @param yd The second column of data values plus detection limits, or a single number representing a standard / guideline value.
#' @param yc The column of censoring indicators for yd, where 1 (or `TRUE`) indicates a detection limit in the yd column, and 0 (or `FALSE`) indicates a detected value in `yd`. Not needed if `yd` is a single standard number.
#' @param alternative The usual notation for the alternate hypothesis.  Default is `“two.sided”`.  Options are `“greater”` or `“less”`.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#'
#' @importFrom coin wilcoxsign_test statistic pvalue
#' @export
#' @return  Prints a list of Wilcoxon Signed-Rank test with Pratt correction for ties statistics containing the following components:
#' \itemize{
#' \item `n` Number of samples
#' \item `Z` The value of the test statistic
#' \item `p.value` the p-value of the test
#' }
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Page, E.B., 1963. Ordered Hypotheses for Multiple Treatments: A Significance Test for Linear Ranks. Journal of the American Statistical Association 58, 216–230. \doi{https://doi.org/10.2307/2282965}
#'
#' Pratt, J.W., 1959. Remarks on Zeros and Ties in the Wilcoxon Signed Rank Procedures. Journal of the American Statistical Association 54, 655–667. \doi{https://doi.org/10.2307/2282543}
#'
#'
#' @examples
#'
#' data(atrazine)
#'
#' cen_signedranktest(atrazine$June,atrazine$JuneCen,atrazine$Sept,atrazine$SeptCen)

cen_signedranktest <- function(xd, xc, yd, yc, alternative="two.sided",printstat=TRUE) {
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

  names(N)<-"n"
  Z<-signif(statistic(s.out), 4)
  names(Z)<-"Z"
  pval<-signif(pvalue(s.out), 4)

  x <- list(n=N,statistic=Z,p.value=pval)
  out<-c(paste(names(x$n),"=",x$n),paste(names(x$statistic),"=",x$statistic),paste("p.value =",x$p.value))

  if(printstat==TRUE){
  cat(txt, "\n","Pratt correction for ties", "\n")
  cat(strwrap(paste(out,collapse=", ")),sep="\n")
  }
  invisible(x)

}

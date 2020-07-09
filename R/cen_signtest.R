#' Sign test for censored data
#'
#'@description Performs a nonparametric sign test of whether the median difference between two columns of paired censored data equals 0. Uses the Pratt adjustment for pairs of equal values.
#'@param xd The first column of data values plus detection limits
#'@param xc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the xd column, and 0 (or `FALSE`) indicates a detected value in xd.
#'@param yd The second column of data values plus detection limits.
#'@param yc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the `yd` column, and 0 (or `FALSE`) indicates a detected value in `yd`.
#'@param alternative The usual notation for the alternate hypothesis.  Default is `“two.sided”`.  Options are `“greater”` or `“less”`.
#'
#' @keywords median difference
#' @export
#' @return  Returns the number of `xd` and `yd` value greater than, less than and tied. Corrected and uncorrected `p-value` for ties also displayed.

#'
#' @examples
#'
#' library(NADA) #for example data
#' data(Atra)
#'
#' cen_signtest(Atra$June,Atra$JuneCen,Atra$Sept,Atra$SeptCen)


cen_signtest <- function(xd, xc, yd, yc, alternative="two.sided") {
  xname <- deparse(substitute(xd))
  yname <- deparse(substitute(yd))
  nonas <- na.omit(data.frame(xd, xc, yd, yc))
  if(length(nonas[,1]) != length(nonas[,3]))
    stop("Lengths of x and y must be the same for paired data.")
  if(alternative == "two.sided")
  {txt3 <- paste("  alternative hypothesis: true median difference is not = 0"); dir <- 0}
  else if(alternative == "less")
  {txt3 <-  paste("  alternative hypothesis: true median difference < 0"); dir <- (-1)}
  else if (alternative == "greater")
  {txt3 <-  paste("  alternative hypothesis: true median difference > 0"); dir <- 1}
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

  # test.dat <- data.frame(x, nonas[2], y, nonas[4])
  # print(test.dat)
  x.gt.y <- sum(x > y)
  x.lt.y <- sum (y > x)
  ties <- length(nonas[,1]) - x.gt.y - x.lt.y
  #  result.counts <- data.frame(x.gt.y, x.lt.y, ties)
  #  print(result.counts)
  s.binom <- binom.test(x.gt.y, length(nonas[,1])-ties, p=0.5, alternative = alternative)

  # LSA version in the coin package.  Not currently used.
  #  s.out <- sign_test(x~y, alternative = alternative)
  #  txt.s <- paste("n =", N, "  Z=", round(statistic(s.out), 4), "  p-value =", round(pvalue(s.out), 4))
  #  print(txt.s)

  # Fong's Modified Sign Test correction for ties
  dir <- dir*sign(x.gt.y - x.lt.y)  #dir=1: alt in same direction as data

  num.Fong <- 1-pbinom(max(x.gt.y, x.lt.y)-1,  N, p=0.5)
  denom.Fong <- 1-pbinom(floor((N-ties+1)/2)-1, N, p=0.5)
  p.Fong <- num.Fong/denom.Fong
  if (dir== 1) {num.Fong <- 1-pbinom((max(x.gt.y, x.lt.y)-1), N, p=0.5)
  denom.Fong <- 2*(1-pbinom(floor((N-ties+1)/2)-1, N, p=0.5))
  p.Fong <- num.Fong/denom.Fong}
  if (dir== -1) p.Fong <- 1-p.Fong/2+dbinom(min(x.gt.y, x.lt.y), N-ties, p=0.5)
  #{num.Fong <- pbinom((max(x.gt.y, x.lt.y)), N, p=0.5)
  # denom.Fong <- 2*(1-pbinom(floor((N-ties+1)/2)-1, N, p=0.5))
  # p.Fong <- num.Fong/denom.Fong}


  txt <- paste("Censored sign test for median(x:", xname, " - ", "y:", yname, ") equals 0", sep = "")
  txt2 <- paste("  n =", N, "  n+ =", x.gt.y, "  n- =", x.lt.y, "   ties:", ties, "\n")
  cat(txt, "\n",txt3, "\n", txt2, "\n", " No correction for ties:", "  p-value =", signif(s.binom$p.value, 4), "\n")
  if (ties != 0) cat("Fong correction for ties:", "  p-value =", signif(p.Fong, 4))
}

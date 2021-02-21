#' Test for difference in left-censored samples
#' @description Performs a nonparametric Paired Prentice-Wilcoxon test of whether the median difference between two columns of paired censored data equals 0 (O'Brien and Fleming, 1987)
#' @param xd The first column of data values plus detection limits
#' @param xc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the xd column, and 0 (or `FALSE`) indicates a detected value in `xd`.
#' @param yd The second column of data values plus detection limits
#' @param yc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the yd column, and 0 (or `FALSE`) indicates a detected value in `yd`
#' @param alternative The usual notation for the alternate hypothesis.  Default is `“two.sided”`.  Options are `“greater”` or `“less”`.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#'
#' @importFrom survival survfit Surv
#' @importFrom stats na.exclude pnorm
#' @return Paired Prentice-Wilcoxon test results including Z-statistic, n (sample size), p-value and median difference
#' @export
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' O’Brien, P.C., Fleming, T.R., 1987. A Paired Prentice-Wilcoxon Test for Censored Paired Data. Biometrics 43, 169–180. https://doi.org/10.2307/2531957
#'
#' @seealso [survival::survfit] [survival::Surv]
#'
#' @examples
#' data(PbHeron)
#' ppw.test(PbHeron$Liver,PbHeron$LiverCen,PbHeron$Bone,PbHeron$BoneCen)
#'

ppw.test <- function(xd, xc, yd, yc, alternative="two.sided",printstat=TRUE)
{
  xname <- deparse(substitute(xd))
  yname <- deparse(substitute(yd))
  OBrienFleming=TRUE
  ## Error checks
  if(length(xd) != length(yd))
    stop("Lengths of x and y must be the same for paired data.")
  keep <- !(is.na(xd) | is.na(yd) | is.na(xc) | is.na(yc))
  x <- xd[keep]
  y <- yd[keep]
  if(any(c(x, y) < 0))
    stop("Negative values in x or y")
  N <- length(x)
  if(OBrienFleming) {
    xd <- x
    xc <- xc[keep]
    yd <- y
    yc <- yc[keep]

    for(i in seq(N)) {
      if(xc[i] && yc[i]) { # Both censored, make same
        xyMax <- max(xd[i], yd[i])
        xd[i] <- yd[i] <- xyMax
      }
      else if(xc[i] && xd[i] > yd[i]) { # x censored > y
        yd[i] <- xd[i]
        yc[i] <- TRUE
      }
      else if(yc[i] && yd[i] > xd[i]) { # y censored > x
        xd[i] <- yd[i]
        xc[i] <- TRUE
      }
    } # Done, no need to check uncensored observations
  }
  ## Test requires stacked data
  group = factor(c(rep(xname, length(x)), rep(yname, length(y))),
                 levels=c(xname, yname))
  values <- c(xd, yd)
  cenflag <- as.logical(c(xc, yc))
  alternative <- pmatch(alternative, c("two.sided", "greater", "less"))
  if(is.na(alternative))
    stop('Invalid choice for alternative, must match "two.sided", "greater", or "less"')
  ## Define the test:
  PPW.test <- function(values, cenflag, group) {
    ## data must be in exactly 2 groups of equal size and stacked
    ## x first, then y (for alternative not two.sided)
    ## required adjustment to values
    adjust <- min(diff(sort(unique(values))))/10.
    values <- ifelse(cenflag, values-adjust, values)
    adjust <- adjust / length(values)
    dupes <- unique(values[!cenflag][duplicated(values[!cenflag])])
    in.values <- values
    if(length(dupes) > 0) { #there are dupes
      for(i in seq(along=dupes)) {
        sel <- which(values == dupes[i])
        mult <- 0L:(length(sel) - 1L)
        in.values[sel] <- dupes[i] + adjust * mult
      }
    }
    ## create data frame and add observed value at 0 to compute
    ## correct probs
    df <- data.frame(values=c(0, -in.values), cenflag=c(T, !cenflag))
    kmout <- survfit(Surv(values, cenflag) ~ 1, data=df, na.action=na.exclude)
    kmout$time <- - kmout$time # convert back to actual values
    ## compute mean survival for tied values
    St <- kmout$surv
    if(length(dupes) > 0L) { #there are dupes
      for(i in seq(along=dupes)) {
        ## take advantage of the fact that the output are sorted
        sel <- which(values == dupes[i])
        mult <- seq(0L, 1L - length(sel))
        sel <- which(kmout$time == dupes[i])
        St[sel] <- mean(St[sel + mult])
      }
    }

    ## Define the link between the observed data and the kaplan meier table
    ## and compute needed stats.
    link <- match(values, kmout$time)
    St <- St[link]
    Uncen <- values[!cenflag]
    UncenSt <- St[!cenflag]
    UncenSt <- UncenSt[order(Uncen)]
    Uncen <- sort(Uncen)
    Score <- 1. - 2*St
    ## This is not fast, but it works
    for(i in which(cenflag)) # fix each each censored score
      Score[i] <- 1. - UncenSt[which(Uncen > values[i])][1L]
    ## Compute d
    ## Reverse sense of d so that alternatives are in same direction--
    ## gives different sense of d from the original
    Score <- matrix(Score, ncol=2L)
    d <-  Score[, 2L] - Score[, 1L]
    return(list(Z=sum(d)/sqrt(sum(d^2)), Scores=Score, Diffs=d))
  } # end of PPW.test

  ret1 <- PPW.test(values, cenflag, group)
  stat <- ret1$Z
  names(stat) <- "Paired Prentice Z"
  meth <- "Paired Prentice-Wilcoxon test"
  param <- length(group)/2L
  names(param) <- "n"
  ## add finishing touches
  if(alternative == 1L)   # alternative is two-sided
  {  pvalue <- (1. - pnorm(abs(stat))) * 2.
     altern <- paste(xname, "not equal to", yname, sep = " ") }
  else if(alternative == 2L)  # alternative is greater than
  {  pvalue <- 1. - pnorm(stat)
     altern <- paste(xname, ">", yname, sep = " ") }
  else # alternaitve is less than
  {  pvalue <- pnorm(stat)
     altern <- paste(xname, "<", yname, sep = " ") }
  names(pvalue) <- "p value"
  mu <- 0
  names(mu) <- "difference"
  Scoremat <- cbind(ret1$Scores, ret1$Diffs)
  colnames(Scoremat) <- c("xScore", "yScore", "d")
  ## For diagnostic plot, create min and max differences
  d1 <- xd - ifelse(yc, 0, yd)
  d2 <- ifelse(xc, 0, xd) - yd
  mind <- pmin(d1, d2)
  maxd <- pmax(d1, d2)
  Scoremat <- cbind(Scoremat, minDiff=mind, maxDiff=maxd)

  sv.out <- survfit(Surv(Scoremat[,4], Scoremat[,5], type="interval2")~1)
  med.diff <- min(sv.out$time [sv.out$surv <= 0.50])

  # getting number of signif digits
  y.count <- vector(length=length(yd))
  for (i in 7:1) {y.count[yd == signif(yd, i)] <- i}
  # computing median difference
  median.diff <- paste ("Median difference equals", signif(med.diff, max(y.count)))

  retval <- list(statistic = stat, parameters = param,
                 p.value = pvalue, null.value = mu,
                 alternative = c("two.sided", "greater", "less")[alternative],
                 method = meth, data.name = paste(xname, "and", yname, sep = " "),
                 PPWmat=Scoremat)
  #oldClass(retval) <- c("htest", "ppw")

#  print(retval)
  txt <- paste("Paired Prentice Wilcoxon test for (x:", xname, " - ", "y:", yname, ") equals 0", "\n", "    alternative: ", altern, "\n", sep = "")
  txt2 <- paste("n =", param, "  Z =", signif(stat, 4), "  p-value =", signif(pvalue, 4))

  if(printstat==TRUE){
  cat(txt, "\n", txt2, "\n")
  cat(" ", median.diff, "\n")
  }
  return(invisible(retval))
}

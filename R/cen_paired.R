#' Censored data paired t-test
#'
#' @description
#'Performs a parametric test of whether the mean difference between two columns of paired censored data equals 0. Assumes that the paired differences follow a gaussian (normal) distribution.
#' @param xd The first column of data values plus detection limits
#' @param xc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the xd column, and 0 (or `FALSE`) indicates a detected value in xd.
#' @param yd The second column of data values plus detection limits, or a single number representing a standard / guideline value.
#' @param yc The column of censoring indicators for yd, where 1 (or `TRUE`) indicates a detection limit in the yd column, and 0 (or `FALSE`) indicates a detected value in `yd`. Not needed if `yd` is a single standard number.
#' @param alternative The usual notation for the alternate hypothesis.  Default is `“two.sided”`.  Options are `“greater”` or `“less”`.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#'
#' @importFrom fitdistrplus plotdistcens
#' @importFrom survival survreg Surv
#' @export
#' @return  A list of statistics containing the following components:
#' \itemize{
#' \item `n` Number of observations
#' \item `Z` The value of the test statistic
#' \item `p.value` the p-value of the test
#' \item `Mean difference` the mean difference between `xd` and `yd`
#' }
#' @details You may also test for whether the mean of the `xd` data exceeds a standard by entering the single number for the standard as `yd`.  In that case no `yc` is required.
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @seealso [survival::survreg]
#'
#' @examples
#'
#' data(atrazine)
#'
#' cen_paired(atrazine$June,atrazine$JuneCen,atrazine$Sept,atrazine$SeptCen)
#'
#' # Comparing standard/guieline value
#' cen_paired(atrazine$June, atrazine$JuneCen, 0.01, alternative = "greater")

cen_paired <- function(xd, xc, yd, yc, alternative="two.sided",printstat = TRUE) {
  xname <- deparse(substitute(xd))
  yname <- deparse(substitute(yd))
  mu <- 0
#  print(xd)
  if (length(yd)==1) {yc <- rep(FALSE, length(xd))
  mu <- yd[1]; yd <- rep(mu, length(xd))}
#  cat(xd, xc, yd, yc)
  nonas <- na.omit(data.frame(xd, xc, yd, yc))
  if(length(nonas[,1]) != length(nonas[,3]))
    stop("Lengths of x and y must be the same for paired data.")
  if(alternative == "two.sided")
  {txt3 <- paste("alternative hypothesis: true mean difference does not equal ", mu, ".", sep="")}
  else if(alternative == "less")
  {txt3 <-  paste("alternative hypothesis: true mean difference is less than ", mu, ".", sep="")}
  else if (alternative == "greater")
  {txt3 <-  paste("alternative hypothesis: true mean difference is greater than ", mu, ".", sep="")}
  else {stop('Invalid choice for alternative, must match "two.sided", "greater", or "less"')}
  N = length(nonas[,1])
  direction <- 0

  ## Create min and max differences
  d1 <- nonas[,1] - ifelse(nonas[,4], 0, nonas[,3])
  d2 <- ifelse(nonas[,2], 0, nonas[,1]) - nonas[,3]
  mind <- pmin(d1, d2)
  maxd <- pmax(d1, d2)
  #  print(data.frame(xd, as.logical(xc), yd, as.logical(yc), mind, maxd))

  ## Mean difference
  sv.out <- survreg(Surv(mind, maxd, type="interval2")~1, dist = "gaussian")
  mean.diff <- sv.out$coefficients

  #direction of mean vs expected
  if (alternative == "less") {direction <- -1}
  if (alternative == "greater") {direction <- 1}
  direction <- direction*sign(mean.diff)
  if (direction <0) {direction <- 0}

  test.diff <- summary(sv.out)
  Z <- test.diff$table[1,3]
  pvalue <- test.diff$table[1,4]
  if(alternative == "two.sided") {
    if (mu == 0) {txt3 <- paste("alternative hypothesis: true mean difference does not equal ", mu, ".", sep="")
    txt <- paste("Censored paired test for mean(", xname, " - ", yname, ") equals ", mu, ".", sep = "")}
    else {txt3 <-  paste("alternative hypothesis: true mean does not equal ", mu, ".", sep ="")
    txt <- paste("Censored paired test for mean(", xname, ") equals ", mu, sep = "")}
  }
  else if(alternative == "less")
  { pvalue <- ifelse (direction, pvalue/2, 1-(pvalue/2))
  if (mu == 0) {txt3 <-  paste("alternative hypothesis: true mean difference is less than ", mu, ".", sep="")
  txt <- paste("Censored paired test for mean(", xname, " - ", yname, ") equals ", mu, ".", sep = "")}
  else {txt3 <-  paste("alternative hypothesis: true mean ", xname, " is less than ", mu, ".", sep="")
  txt <- paste("Censored paired test for mean(", xname, ") equals ", mu, sep = "")}
  }
  else   # (alternative == "greater")
  { pvalue <- ifelse(direction, pvalue/2, 1-(pvalue/2) )
  if (mu == 0) {txt3 <-  paste("alternative hypothesis: true mean difference is greater than ", mu, ".", sep="")
  txt <- paste(" Censored paired test for mean(", xname, " - ", yname, ") equals ", mu, ".", sep = "")}
  else {txt3 <-  paste("alternative hypothesis: true mean ", xname, " exceeds ", mu, ".", sep = "")
  txt <- paste(" Censored paired test for mean(", xname, ") equals ", mu, sep = "")}
  }
  txt2 <- paste("n =", N, "  Z=", round(Z, 4), "  p-value =", signif(pvalue, 4))

  if(printstat==TRUE){
    if (mu == 0) {cat( txt, "\n", txt3, "\n", "\n", txt2, "\n", "Mean difference =", signif(mean.diff,4), "\n")}
    else {cat( txt, "\n", txt3, "\n", "\n", txt2, "\n",paste("Mean", xname), "=", signif(mean.diff+mu,4), "\n")}
  }

  names(N)<-"n"
  Z<-round(Z,4)
  names(Z)<-"Z"
  pval<-signif(as.numeric(pvalue), 4)
  diff.val<-as.numeric(signif(mean.diff+mu,4))
  names(diff.val)<-paste("Mean", xname)

  x <- list(n=N,statistic=Z,p.value=pval,MeanDifference=diff.val)


  #plot the cdf of differences
  sv.coefs <- data.frame(sv.out$icoef[1], exp(sv.out$icoef[2]))
  names(sv.coefs) <- c("mean", "sd")
  left <- mind;  right <- maxd
  minmax <- data.frame(left, right)

  oldpar<- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  plotdistcens(minmax, distr = "norm", para = sv.coefs, main = "Differences: CDF and Fitted Normal Dist")
  invisible(x)
}

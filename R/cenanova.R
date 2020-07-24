#' ANOVA for censored data
#'
#' @description Performs a parametric test of differences in means between groups of censored data, followed by a parametric Tukey's multiple comparison test.
#' @param y1 The column of data values plus detection limits
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y1` column, and 0 (or `FALSE`) indicates a detected value in `y1`.
#' @param grp Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param LOG Indicator of whether to compute tests in the original units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @keywords ANOVA
#' @export
#' @return
#' Returns the Maximum Likelihood Estimation (MLE) comparison results including Chi-Squared value, degrees of freedom and `p-value` of the test. Test assumes lognormal(`LOG=TRUE`) or nomal(`LOG=FALSE`) distribution
#'
#' In addition to comparison of censored data between groups, multiple comparisons of means are also printed based on the parametric survival model (`CensData~Factor`).
#' \itemize{
#' \item Group Names of groups (NOTE: `==` indicates `"two.sided"` comparisons`)
#' \item `Estimate` Slope of the parametric survival model
#' \item `Std. Error` Standard error of estimate
#' \item `z value` z value test statistic
#' \time `Pr(>|z|)` P-values of test
#' }
#'
#' @details When this happens the p-values may be unreal (often lower than they should be).  Because of this, testing in log units is preferable and the default.
#' @importFrom survival survreg Surv
#' @importFrom multcomp glht
#' @seealso [survival::survreg]
#'
#' @examples
#'
#' library(NADA) #for example data
#' data(Golden)
#' cenanova(Golden$Liver,Golden$LiverCen,Golden$DosageGroup)
#'
#' cenanova(Golden$Liver,Golden$LiverCen,Golden$DosageGroup,LOG=FALSE)
#'
#'
cenanova <- function(y1, y2, grp, LOG=TRUE) {
  yname <- deparse(substitute(y1))
  gname <- deparse(substitute(grp))

  # for both log and original units
  detect <- as.logical(1 - as.integer(y2))  # reverses TRUE/FALSE to fit survival functions
  Factor <- as.factor(grp)
  df <- length(levels(Factor))-1
  grpnames <- as.character(levels(Factor))

  # ln units for LOG = TRUE
  if (LOG == TRUE)  {
    lnvar <- log(y1)
    fconst <- max(lnvar)+1
    flip.log <- fconst - lnvar
    logCensData <- Surv(flip.log, detect, type="right")
    reg.out <- survreg(logCensData~Factor, dist = "gaussian")
    reg.chisq <- -2*(reg.out$loglik[1] - reg.out$loglik[2])
    reg.out$coefficients <- (-1)* reg.out$coefficients  # reversing the flip for all coeffs.
    reg.out$coefficients[1] <- fconst + reg.out$coefficients[1]  #reversing the flip for intercept
    dist.test <- "Assuming lognormal distribution of CensData"
    pval = pchisq(reg.chisq, df, lower.tail = FALSE)
    dist <- "Lognormal Dist";  statistic <- reg.chisq
    result <- data.frame(dist, statistic, df, pval)

    #  write test results
    cat('\n',"     MLE test of mean natural logs of CensData:", yname, "by Factor:", gname, '\n', "    ",dist.test,'\n', "     Chisq =", signif(reg.chisq, 4), " on", df, "degrees of freedom", "    p =", signif(pval,3), '\n', '\n')
    # multiple comparisons
    x.mc <- glht(reg.out, linfct = mcp(Factor = "Tukey"))

    # Q-Q plot of residuals
    reg.predict <- predict(reg.out)
    log.unit <- reg.predict - flip.log
    cenregQQ(log.unit, as.logical(y2), Factor, main = "Normal Q-Q Plot of Logscale Residuals", LOG = FALSE)
  }
  else  {
    # no logs used.  Original units
    CensData <- Surv(y1, detect, type="left")
    reg.out <- survreg(CensData~Factor, dist = "gaussian")   # officially a Tobit model
    reg.chisq <- -2*(reg.out$loglik[1] - reg.out$loglik[2])
    dist.test <- "Assuming normal distribution of CensData"
    pval = pchisq(reg.chisq, df, lower.tail = FALSE)
    dist <- "Normal Dist";  statistic <- reg.chisq
    result <- data.frame(dist, statistic, df, pval)

    #  write test results
    cat('\n',"     MLE test of mean CensData:", yname, "  by Factor:", gname, '\n', "    ",dist.test,'\n', "     Chisq =", signif(reg.chisq, 4), " on", df, "degrees of freedom", "    p =", signif(pval,3), '\n')
    # A warning
    cat("\n", "  NOTE: Data with nondetects may be projected below 0 with MLE normal distribution.", "\n", "  If so, p-values will be unreliable (often too small).  Use perm test instead.", "\n", '\n')

    # multiple comparisons
    x.mc <- glht(reg.out, linfct = mcp(Factor = "Tukey"))

    # Q-Q plot of residuals
    # reg.predict <- reg.out$linear.predictors
    from.groupmean <- residuals(reg.out, type = "response")
    cenregQQ(from.groupmean, as.logical(y2), Factor, LOG = FALSE)
  }

  # finish for both log and not log
  group.means <- as.vector(x.mc$coef)
  for (i in 2:length(levels(Factor))) {group.means[i] <- group.means[i] + group.means[1]}
  group.means <- t(group.means)
  mean.names <- paste("mean(", grpnames, ")", sep="")
  colnames(group.means) <- mean.names
  rownames(group.means) <- " "
  result <- cbind(result, group.means)

  #print group means and mult comparison results
  print(group.means, row.names = FALSE, print.gap = 3)
  print(summary(x.mc))

  return(invisible(result))
}

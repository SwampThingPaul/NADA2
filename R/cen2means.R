#' Censored data two-group test for difference in means
#'
#' @description
#' Performs a parametric test of differences in means between two groups of censored data, either in original or in log units (the latter becomes a test for difference in geometric means).
#' @param x1 The column of data values plus detection limits
#' @param x2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the y1 column, and 0 (or `FALSE`) indicates a detected value in y1.
#' @param group Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param LOG Indicator of whether to compute tests in the original units, or on their logarithms.  The default is to use the logarithms (LOG = `TRUE`).  To compute in original units, specify the option LOG = `FALSE` (or LOG = 0).
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @importFrom stats pchisq predict
#' @export
#' @return
#' Q-Q Plot with Shapiro-Francia test for normality W and p-values.
#' Returns the Maximum Likelihood Estimation (MLE) test results including Chi-Squared value, degrees of freedom and `p-value` of the test.
#'
#' @details Because this is an MLE procedure, when a normal distribution model is used (LOG=FALSE) values may be modeled as below zero.  When this happens the means may be too low and the p-values may be unreal (often lower than they should be).  Because of this, testing in log units is preferable and is the default.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @importFrom survival survreg Surv
#'
#' @examples
#'
#' data(PbHeron)
#' cen2means(PbHeron$Liver,PbHeron$LiverCen,PbHeron$DosageGroup)


cen2means <- function(x1, x2, group, LOG=TRUE,printstat=TRUE) {
  yname <- deparse(substitute(x1))
  gname <- deparse(substitute(group))

  ydat <- na.omit(data.frame(x1, x2, group))
  y1 <- ydat[,1];  y2 <- ydat[,2]; grp <- ydat[,3]

  # original units for LOG = FALSE
  fconst <- max(y1)
  flip <- fconst - y1
  # for both log and original units
  detect <- as.logical(1 - as.integer(y2))  # reverses TRUE/FALSE to fit survival functions
  Factor <- as.factor(grp)
  df <- length(levels(Factor))-1
  grpname <- as.character(levels(Factor))

  # ln units for LOG = 1
  if (LOG == TRUE)  {
    lnvar <- log(y1)
    fconst <- max(lnvar)
    flip.log <- max(lnvar) - lnvar

    logCensData <- Surv(flip.log, detect, type="right")
    reg.out <- survreg(logCensData~Factor, dist = "gaussian")
    reg.chisq <- -2*(reg.out$loglik[1] - reg.out$loglik[2])
    # print(reg.out$coefficients)   # before unflipping
    reg.out$coefficients <- (-1)* reg.out$coefficients
    reg.out$coefficients[1] <- fconst + reg.out$coefficients[1]  #reversing the flip
    #  print(reg.out$coefficients)
    dist.test <- "Assuming lognormal distribution of residuals around group geometric means"
    pval = pchisq(reg.chisq, df, lower.tail = FALSE)
    mean1 <- exp(reg.out$coefficients[1])
    mean2 <- exp(reg.out$coefficients[1] + reg.out$coefficients[2])
    dist <- "Lognormal Dist";  statistic <- reg.chisq
    result <- data.frame(dist, statistic, df, pval)

    #  write results
    if(printstat==TRUE){
      cat("     MLE 't-test' of mean natural logs of CensData:", yname, "by Factor:", gname, '\n', "   ",dist.test,'\n')
      cat("     geometric mean of", grpname[1], "=", signif(mean1, 4), "    geometric mean of", grpname[2], "=", signif(mean2,4), "\n")
      cat( "     Chisq =", signif(reg.chisq, 4), " on", df, "degrees of freedom", "    p =", signif(pval,3), '\n', "\n")
    }
    # Q-Q plot of residuals
    reg.predict <- predict(reg.out)
    two.group <- exp(reg.predict - flip.log)
    cenregQQ(two.group, as.logical(y2), Factor, LOG = TRUE)
  }
  else  # no logs used.  Better to use cenperm2 permutation test instead.
  { CensData <- Surv(flip, detect, type="right")
  reg.out <- survreg(CensData~Factor, dist = "gaussian")
  reg.out$coefficients <- (-1)* reg.out$coefficients  #reversing the flip
  reg.out$coefficients[1] <- fconst + reg.out$coefficients[1]  #reversing the flip for the intercept
  # print(reg.out$coefficients)
  reg.chisq <- -2*(reg.out$loglik[1] - reg.out$loglik[2])
  mean1 <- reg.out$coefficients[1]
  mean2 <- mean1 + reg.out$coefficients[2]

  dist.test <- "Assuming normal distribution of residuals around group means"
  pval = pchisq(reg.chisq, df, lower.tail = FALSE)
  dist <- "Normal Dist";  statistic <- reg.chisq
  result <- data.frame(dist, statistic, df, pval)

  #  write results
  if(printstat==TRUE){
    cat("     MLE 't-test' of mean CensData:", yname, "  by Factor:", gname, '\n', "   ",dist.test,'\n')
    cat("     mean of", grpname[1], "=", signif(mean1, 4), "    mean of", grpname[2], "=", signif(mean2,4), "\n")
    cat("     Chisq =", signif(reg.chisq, 4), " on", df, "degrees of freedom", "    p =", signif(pval,3), '\n')
    # A warning
    warning(paste("NOTE: Data with nondetects may be projected below 0 with MLE normal distribution.", "\n", "  If so, p-values will be unreliable (often too small).  Use perm test instead.", "\n"))
  }
  # Q-Q plot of residuals
  reg.predict <- predict(reg.out)
  two.group <- reg.predict - flip
  cenregQQ(two.group, as.logical(y2), Factor, LOG = FALSE)
  }
  return(invisible(result))
}

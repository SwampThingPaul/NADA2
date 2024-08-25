#' Compute an ECDF and Distribution Parameters for Censored Data
#'
#' @description Computes the empirical cumulative distribution function (ECDF) for censored data. Estimates parameters of the distribution, including the mean and quantiles.
#' @param y Concentrations plus detection limits for indicator formatted data.
#' @param cens Censoring indicators (logical. 1 or `TRUE` = censored, 0 or FALSE = detected) for indicator formatted data.
#' @param conf The confidence coefficient for confidence intervals around the Kaplan-Meier mean and median. Default = 0.95.
#' @param qtls Probabilities for the quantiles to be estimated.  Defaults are (0.10, 0.25, 0.50, 0.75, 0.90).  You may add and/or substitute probabilities -- all must be between and not including 0 to 1.
#' @param Cdf Logical `TRUE`/`FALSE` indicator of whether to plot the empirical cumulative distribution function (cdf).
#' @param ecdf.col Color for the ecdf plotted step function line.  Default is black.
#' @param km.orig  If `TRUE` (default), Kaplan-Meier results in the realm below detection limits reported as "NA".  If `FALSE`, information in the detection limits is used and results in the realm of detection limits reported as "< DL", where DL is the appropriate detection limit.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @param Ylab Optional input text in quotes to be used as the variable name on the ecdf plot.  The default is the name of the `y1` input variable.
#' @param q.type an integer between 1 and 9 selecting one of the nine quantile algorithms used only when km.orig = FALSE. See `stats::quantile` and Details below for more detail, default when km.orig = FLASE is set to 6.
#' @importFrom survival Surv survfit
#' @importFrom stats quantile
#' @return
#' If `printstat=TRUE`: Based on the provided `conf` value, Kaplan-Meier summary statistics (`mean`,`sd`,`median`), lower and upper confidence intervals around the mean and median value, sample size and percent of censored samples are returned. The specified quantile values are also printed and returned.
#'
#' If `Cdf=TRUE`: The ecdf of censored data is plotted.
#' @details Moment statistics are estimated using the enparCensored function of the EnvStats package. This avoids a small bias in the mean produced by the NADA package's cenfit function, which uses the reverse Kaplan-Meier procedure, converting left-censored to right-censored data prior to computing the ecdf and mean. See Gillespie et al.(2010) for more discussion on the bias of the estimate of the mean.
#'
#' Quantiles and their two-sided confidence limits are estimated using the quantile function of the survfit command. See ?quantiles or Helsel et al. (2020) for choosing the q.type; default q.type = 4 (Kaplan-Meier; prob = i/n) when km.orig = TRUE. This is standard procedure in the survival analysis discipline and in the survival package of R, and is also used by the cenfit function in the NADA package. While this is 'industry standard' in medical applications it is a poor choice for observed sample (rather than census) data because it means that the largest observation is assigned a probability equal to 1, the 100th percentile. This implies that this value is never expected to be exceeded when more sample data are collected. It also means the largest observation is not plotted on the ecdf in most software because it is "off the chart". For small datasets in particular it is unlikely that the current largest observation is the largest value in the population and so the resulting ecdf quantiles are likely not opotimal. The default q.type = 6 (Weibull; prob = i/(n+1)) when km.orig = FALSE, though that may be changed by the user. the largest observation plots at a probability less than 1 on the ecdf.  Differences in results when differing q.types are used will decrease as sample size increases.
#'
#' All printed values will also be output to an object if saved. Confidence intervals on the quantiles are also output when data include nondetects. Values are character because of the possibility of a `<1`, but if no `<` symbol can be converted to numeric using the `as.numeric(...)` command.  For data without censoring cfit will return numeric values. In that case it returns standard arithmetic mean, standard deviation and quantiles instead of K-M versions.
#'
#' @export
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Gillespie, B.W., et al., 2010.  Estimating Population Distributions When Some Data Are Below a Limit of Detection by Using a Reverse Kaplan-Meier Estimator. Epidemiology 21, 564-570.
#'
#' Helsel, D.R., Hirsch, R.M., Ryberg, K.R., Archfield, S.A., and Gilroy, E.J., 2020, Statistical Methods in Water Resources: U.S. Geological Survey Techniques and Methods, book 4, chapter A3, 458 p., https://doi.org/10.3133/tm4a3.
#'
#' Millard, S.P, 2013. EnvStats: An R Package for Environmental Statistics, 2nd ed. Springer Science+Business Media, USA, N.Y.  DOI 10.1007/978-1-4614-8456-1© Springer Science+Business Media New York 2013”
#'
#' Excerpt From: Steven P. Millard. “EnvStats.” Apple Books.
#'
#' @seealso [survival::survfit]
#' @examples
#'
#' data(Brumbaugh)
#'
#' cfit(Brumbaugh$Hg,Brumbaugh$HgCen)
#'

cfit <- function(y, cens, conf=0.95, qtls = c(0.10, 0.25, 0.50, 0.75, 0.90),
                 q.type = 6, Cdf = TRUE, ecdf.col = 1, km.orig = TRUE,
                 printstat = TRUE, Ylab = NULL){

  if (is.null(Ylab)) {Ylab <- deparse(substitute(y))}
  Conf <- 100*conf
  if (isFALSE(0.50 %in% qtls)) {qtls <- append(0.50, qtls)}     # median always needed
  k <- length(qtls)
  median.pos <- which(qtls == 0.50)
  ecdf.qtls <- vector("numeric", k)
  ecdf.ranks <- vector("integer", k)
  ecdf.yvals <- vector("numeric", k)
  max.tied <- vector("numeric", k)
  dl.tied <- vector("character", k)

  # remove lines with missing values.  Sort by increasing y
  nonas.1 <- na.omit(data.frame(y, cens))
  nonas  <- nonas.1[order(nonas.1[,1]),]
  y1 <- nonas[,1] ; y2 <- nonas[,2]
  y2 <- as.numeric(y2)
  N <- length(y2)
  names(N) <- " "
  delta = y2*min(y1)*(0.001)
  y1 <- y1-delta             # y1 nondetects slightly smaller than detects at same value

  #  if(is.logical(y2) | prod(as.numeric(y2)) == 0)   # Some zeros in y2 == indicator format for y2.
  #  reactivate previous line to allow a 2nd section for interval censoring in the future.

  overall.min <- min(nonas[,1])
  detect.min <- min(nonas[,1][y2==0])

  # if some values are nondetects
  if (sum(y2) != 0)   { DLs <- nonas[,1][y2==1]
  nd.max <- max(nonas[,1][y2==1])
  PctND <- signif(100*length(DLs)/N, 4)
  ycen2 <- 1-y2  # indicator for right-censored Surv
  flip.const <- max(y1)   # use the max for the constant
  yhi.flip <- flip.const - y1
  y.surv <- Surv(yhi.flip, ycen2)
  y.out<- survfit(y.surv ~ 1, conf.int = conf, conf.type = "plain")

  # left.censored.min="DL" is the Efron bias correction, generally standard for Kaplan-Meier
  if(is.logical(y2)==F){y2 <- as.logical(y2)}
  npar.out <- enparCensored(x=y1, censored = y2, ci = TRUE, pivot.statistic = "t", ci.sample.size = "Total", left.censored.min="DL")
  KMmean <- npar.out$parameters[1]
  std.err <- npar.out$parameters[3]
  KMsd <- signif(npar.out$parameters[2], 4)
  # could  use enparCensored output $limits [1,2] for conf limits on mean instead of computation below.  Get same numbers
  qt.ci <- qt(c((1-conf)/2, 1-(1-conf)/2), N-1)
  LCLmean <- signif(KMmean + qt.ci[1]*std.err, 4)
  UCLmean <- signif(KMmean + qt.ci[2]*std.err, 4)
  KMmean <- signif(KMmean, 4)

  # using EnvStats::ecdfPlotCensored produces the plot and computes quantiles for Kaplan-Meier.
  ecdf.out <- ecdfPlotCensored(y1, y2, plot.it = Cdf, plot.pos.con = 0.375, prob.method = "kaplan-meier", type = "s", ecdf.col = ecdf.col, xlab = Ylab, main = "Empirical CDF of Censored Data", ecdf.lwd = par("cex"))

  for (i in 1:k)  {ecdf.qtls[i] <- min(ecdf.out$Cumulative.Probabilities[ecdf.out$Cumulative.Probabilities >= qtls[i]])
  ecdf.ranks[i] <- min(which(ecdf.out$Cumulative.Probabilities == ecdf.qtls[i]))
  ecdf.yvals[i] <- nonas[,1][ecdf.ranks[i]]
  # these are standard Kaplan-Meier quantiles
  # gives the lowest DL for a series of censored data.
  }

  # get conf int on median by quantile function on a surv object (y.out).
  # requires unflipping the results using 1-qtls for probabilities.
  surv.median <- quantile(y.out, probs = 0.5, conf.int = TRUE)
  names(ecdf.yvals) <- paste (as.character(100*(qtls[1:k])), sep = "")
  surv.quantiles <- quantile(y.out, probs = 1-qtls, conf.int = TRUE, type=q.type)
  Qtiles <- flip.const - surv.quantiles$quantile    # numeric

  KMmedian <- signif(ecdf.yvals[median.pos], 4)
  Qtls <- as.character(signif(ecdf.yvals, 4))       # character

  LCLmedian <- signif(flip.const - surv.median$upper, 4)  # possibly == NA
  UCLmedian <- signif(flip.const - surv.median$lower, 4)   # possibly == NA
  Lower <- signif(flip.const - surv.quantiles$upper, 4)        # numeric
  names(Lower) <- paste (as.character(100*(qtls[1:k])), sep = "")
  Upper <- signif(flip.const - surv.quantiles$lower, 4)       # numeric
  names(Upper) <- paste (as.character(100*(qtls[1:k])), sep = "")

  # NAs in quantiles changed to the highest detection limit for a series of DLs with idential cumulative probabilites.  Provides more info than standard K-M.
  # the quantile function below gives the highest DL for a given cumulative probability. This is what should be used for a "<DL" result.
  if (km.orig == FALSE) {
    surv.quantiles <- quantile(y.out, probs = 1-qtls, conf.int = TRUE, type=q.type)
    Qtiles <- flip.const - surv.quantiles$quantile    # numeric
    names(Qtiles) <- paste (as.character(100*(qtls[1:k])), sep = "")
    Lower <- signif(flip.const - surv.quantiles$upper, 4)        # numeric
    names(Lower) <- paste (as.character(100*(qtls[1:k])), sep = "")
    Upper <- signif(flip.const - surv.quantiles$lower, 4)       # numeric
    names(Upper) <- paste (as.character(100*(qtls[1:k])), sep = "")
    Qtls <- as.character(signif(Qtiles, 4))    #character
    KMmedian <- signif(Qtiles[median.pos], 4)    #character
    LCLmedian <- as.character(signif(Lower[median.pos], 4))  # character, possibly == NA
    UCLmedian <- as.character(signif(Upper[median.pos], 4))   # character, possibly == NA

    for (i in 1:k) {max.tied[i] <- max(which(ecdf.out$Cumulative.Probabilities == ecdf.qtls[i]))

    if (is.na(Qtiles[i])) {Qtiles[i] <- nonas[,1][max.tied[i]]
    Qtls[i] <- paste("<", Qtiles[i], sep="")
    KMmedian <- Qtls[median.pos]}

    # getting <DL values for Lower CI endpoints
    if(is.na(Lower[i])) {Lower[i]<- max(DLs[DLs <= Qtiles[i]])
    Lower[i] <- as.character(paste("<", Lower[i], sep="") )
    LCLmedian <- Lower[median.pos] }

    # getting reality for upper CI endpoints
    if (is.na(Upper[i])) {
      Upper[i] <- "unknown"}
    }  # end of loop for establishing quantiles and CIs if <DL
    UCLmedian <- Upper[median.pos]
  }   # end of when km.orig == FALSE
  #  print(nonas[,1]); print(nonas[,2]); print(ecdf.ranks); print(max.tied)
  }  # end of data that include nondetects

  # when all data are detects
  else  {Mean <- mean(y1, na.rm = TRUE)
  PctND <- 0
  SD <- sd(y1, na.rm = TRUE)
  std.err <- SD/sqrt(length(na.omit(y1)))
  qt.ci <- qt(c((1-conf)/2, 1-(1-conf)/2), N-1)
  LCLmean <- signif(Mean + qt.ci[1]*std.err, 4)
  UCLmean <- signif(Mean + qt.ci[2]*std.err, 4)
  Mean <- signif(Mean, 4)
  SD <- signif(SD, 4)
  Qtls <- signif(quantile(y1, probs = qtls, na.rm = TRUE, type = q.type), 4)
  Median <-  signif(quantile(y1, probs = 0.50, na.rm = TRUE, type = q.type), 4)
  ecdf.out <- ecdfPlot(y1, plot.it = Cdf, plot.pos.con = 0.375, type = "s", xlab = Ylab, main = paste("Empirical CDF of", deparse(substitute(y))), ecdf.col = ecdf.col, ecdf.lwd = par("cex"))

  # all data are detects: getting CIs for median
  orderstat.CI <- quantile_CI (N, 0.50, alpha = 1-conf)
  order.y <- y1[sort.list(y1)]
  CI.median <- order.y[orderstat.CI$Interval]
  LCLmedian <- signif(CI.median[1], 4)
  UCLmedian <- signif(CI.median[2], 4)

  # end of when all data are detects
  }

  #print results
  quant.char <- paste ("Q", as.character(100*(qtls[1:k])), sep = "")
  Cpct <- paste(Conf, "%", sep = "")
  if (sum(y2) != 0) {stats <- data.frame (N, PctND, Conf, KMmean, KMsd, LCLmean, UCLmean, KMmedian, LCLmedian, UCLmedian)}
  else {stats <- data.frame (N, PctND, Conf, Mean, SD, LCLmean, UCLmean, Median, LCLmedian, UCLmedian)}
  if (printstat == TRUE) {cat("\n", "Output for", Ylab, "            ", Cpct, "Confidence Intervals", "\n", "Statistics:", "\n")
    print(stats)
    cat("\n")
    cat("Quantiles:", quant.char, "\n", sep="\t")
    cat("        ", Qtls, "\n", "\n", sep="\t")
  }

  names(Qtls) <- quant.char
  if (sum(y2) != 0) {Qtls.out <- cbind(Lower, Qtls, Upper)}
  else {Qtls.out <- Qtls }
  cstats <- list(stats, Qtls.out)
  return (invisible (cstats))
}



#' quantile confidence interval
#' @param n length
#' @param q quantile
#' @param alpha CI value
#' @keywords internal
#' @importFrom stats qbinom
#' @return Returns the order statistics and the actual coverage.


quantile_CI <- function(n, q, alpha=0.05) {
  #
  # Search over a small range of upper and lower order statistics for the
  # closest coverage to 1-alpha (but not less than it, if possible).
  #
  u <- qbinom(1-alpha/2, n, q) + (-2:2) + 1
  l <- qbinom(alpha/2, n, q) + (-2:2)
  u[u > n] <- Inf
  l[l < 0] <- -Inf
  coverage <- outer(l, u, function(a,b) pbinom(b-1,n,q) - pbinom(a-1,n,q))
  if (max(coverage) < 1-alpha) i <- which(coverage==max(coverage)) else
    i <- which(coverage == min(coverage[coverage >= 1-alpha]))
  i <- i[1]

  u <- rep(u, each=5)[i]
  l <- rep(l, 5)[i]
  return(list(Interval=c(l,u), Coverage=coverage[i]))

}

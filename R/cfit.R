#' Compute an ECDF for Censored Data
#'
#' @description Computes an estimate of an empirical cumulative distribution function (ECDF) for censored data
#' @param y1 lowest possible concentration (interval-censored) or concs+DLs for indicator format
#' @param y2 highest possible concentration (interval-censored) or censoring indicator (logical. 1 or `TRUE` = censored, 0 or FALSE = detected) for indicator format.
#' @param conf the confidence coefficient for confidence intervals around the Kaplan-Meier mean and median.
#' @param qtls Probabilities for the quantiles to be estimated.  The defaults are listed above.  You may add and/or substitute probabilities.
#' @param Cdf Logical `TRUE`/`FALSE` Indicator of whether to plot the data a cumulative distribution function (cdf) using Kaplan-Meier quantiles.
#' @param printstats Logical `TRUE`/`FALSE` Option of whether to print the resulting statisics in the console window, or not
#' @param Ylab Optional – input text in quotes to be used as the variable name on the cdf plot.  The default is the name of the `y1` input variable.
#'
#' @importFrom survival Surv survfit
#' @return
#'
#' If `printstats=TRUE` Based on the provided `conf` value, Kaplan-Meier summary statistics (`mean`,`sd`,`median`), lower and upper confidence intervals around the mean and median value, sample size and percent of censored samples. Additionally specified quartile values are returned.
#'
#' If `Cdf=TRUE` eCDF plot of censored data is printed.
#'
#' @export
#' @references
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#' @seealso [survival::survfit] [NADA::cenfit] [NADA::cendiff]
#' @examples
#' library(NADA) #For example data
#'
#' data(Cadmium)
#'
#' cfit(Cadmium$Cd,Cadmium$CdCen)
cfit <- function(y1, y2, conf=0.95, qtls = c(0.10, 0.25, 0.50, 0.75, 0.90), Cdf = TRUE, printstats = TRUE, Ylab = NULL) {
  N <- length(y2)
  if (is.null(Ylab)) Ylab <- deparse(substitute(y1))
  Conf <- 100*conf
  k <- length(qtls)
  names(N) <- " "

  if(is.logical(y2) | prod(y2) == 0)
  { y2 <- as.numeric(y2)
  DLs <- y1[y2==1]
  PctND <- signif(100*length(DLs)/N, 4)

  overall.min <- min(y1)
  yzero <- y1*(1-y2)
  flip.const <- max(y1)
  yhi.flip <- flip.const - y1
  ylo.flip <- flip.const - yzero
  y.surv <- Surv(yhi.flip, ylo.flip, type = "interval2")
  y.out<- survfit(y.surv ~ 1, conf.int = conf, conf.type = "plain")
  Rmean.flip <- survival:::survmean(y.out, rmean=flip.const) [[1]]["*rmean"]
  KMmean <- signif(flip.const - Rmean.flip, 4)
  std.err <- survival:::survmean(y.out,  rmean=flip.const) [[1]]["*se(rmean)"]
  KMsd <- signif(std.err*sqrt(N), 4)
  qt.ci <- qt(c((1-conf)/2, 1-(1-conf)/2), N-1)
  LCLmean <- signif(KMmean + qt.ci[1]*std.err, 4)
  UCLmean <- signif(KMmean + qt.ci[2]*std.err, 4)

  LCLmedian <- signif(flip.const - survival:::survmean(y.out,  rmean=flip.const) [[1]]["0.95UCL"], 4)
  UCLmedian <- signif(flip.const - survival:::survmean(y.out,  rmean=flip.const) [[1]]["0.95LCL"], 4)
  KMmedian <- signif(flip.const - survival:::survmean(y.out,  rmean=flip.const) [[1]]["median"], 4)

  KMmedian <- ifelse(KMmedian<overall.min, paste("<",overall.min, sep=""), KMmedian)
  flip.out <- y.out
  flip.out$time <- flip.const - y.out$time
  Qtls <- signif(quantile(flip.out, probs = (1-qtls), conf.int = FALSE), 4)

  Qtls <- ifelse(Qtls<overall.min, paste("<",overall.min, sep=""), Qtls)

  Qtls.char <- as.character(Qtls)
  quant.char <- paste ("Q", as.character(100*(qtls[1:k])), sep = "")
  Cpct <- paste(Conf, "%", sep = "")
  stats <- data.frame (N, PctND, Conf, KMmean, KMsd, KMmedian, LCLmean, UCLmean, LCLmedian, UCLmedian)

  if (printstats == TRUE) {cat("\n", "Output for", Ylab, "            ", Cpct, "Confidence Intervals", "\n", "Statistics:", "\n")
    print(stats)
    cat("\n")
    cat("Quantiles:", quant.char, "\n", sep="\t")
    cat("        ", Qtls.char, "\n", "\n", sep="\t")
  }

  names(Qtls) <- quant.char
  Qtls <- t(Qtls)

  cstats <- cbind(stats, Qtls)

  oldpar<- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (Cdf == TRUE) plot(flip.out, xlab = Ylab, ylab = "Cumulative Probability", ylim = c(0,1))
  return (invisible (cstats))
  }

  else {
  }
}


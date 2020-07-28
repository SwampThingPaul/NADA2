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
# @importFrom survival Surv survfit
#' @import survival
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
#'

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
  Rmean.flip <- NADA2.survmean(y.out, rmean=flip.const) [[1]]["*rmean"]
  KMmean <- signif(flip.const - Rmean.flip, 4)
  std.err <- NADA2.survmean(y.out,  rmean=flip.const) [[1]]["*se(rmean)"]
  KMsd <- signif(std.err*sqrt(N), 4)
  qt.ci <- qt(c((1-conf)/2, 1-(1-conf)/2), N-1)
  LCLmean <- signif(KMmean + qt.ci[1]*std.err, 4)
  UCLmean <- signif(KMmean + qt.ci[2]*std.err, 4)

  LCLmedian <- signif(flip.const - NADA2.survmean(y.out,  rmean=flip.const) [[1]]["0.95UCL"], 4)
  UCLmedian <- signif(flip.const - NADA2.survmean(y.out,  rmean=flip.const) [[1]]["0.95LCL"], 4)
  KMmedian <- signif(flip.const - NADA2.survmean(y.out,  rmean=flip.const) [[1]]["median"], 4)

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


#' Summary statistics of survival curve (from `survival:::survmean`)
#' @param x the result of a call to the survfit function.
#' @param scale a numeric value to rescale the survival time, e.g., if the input data to survfit were in days, scale=365 would scale the printout to years.
#' @param rmean restricited mean
#' @importFrom stats median
#' @keywords internal
#'
#' @export
NADA2.survmean=function(x, scale = 1, rmean)
{
  # Extracted from survival::survmean
  if (!is.null(x$start.time))
    start.time <- x$start.time
  else start.time <- min(0, x$time)
  pfun <- function(nused, time, surv, n.risk, n.event, lower,
                   upper, start.time, end.time) {
    minmin <- function(y, x) {
      tolerance <- .Machine$double.eps^0.5
      keep <- (!is.na(y) & y < (0.5 + tolerance))
      if (!any(keep))
        NA
      else {
        x <- x[keep]
        y <- y[keep]
        if (abs(y[1] - 0.5) < tolerance && any(y < y[1]))
          (x[1] + x[min(which(y < y[1]))])/2
        else x[1]
      }
    }
    if (!is.na(end.time)) {
      hh <- ifelse((n.risk - n.event) == 0, 0, n.event/(n.risk *
                                                          (n.risk - n.event)))
      keep <- which(time <= end.time)
      if (length(keep) == 0) {
        temptime <- end.time
        tempsurv <- 1
        hh <- 0
      }
      else {
        temptime <- c(time[keep], end.time)
        tempsurv <- c(surv[keep], surv[max(keep)])
        hh <- c(hh[keep], 0)
      }
      n <- length(temptime)
      delta <- diff(c(start.time, temptime))
      rectangles <- delta * c(1, tempsurv[-n])
      varmean <- sum(cumsum(rev(rectangles[-1]))^2 * rev(hh)[-1])
      mean <- sum(rectangles) + start.time
    }
    else {
      mean <- 0
      varmean <- 0
    }
    med <- minmin(surv, time)
    if (!is.null(upper)) {
      upper <- minmin(upper, time)
      lower <- minmin(lower, time)
      c(nused, max(n.risk), n.risk[1], sum(n.event), sum(mean),
        sqrt(varmean), med, lower, upper)
    }
    else c(nused, max(n.risk), n.risk[1], sum(n.event),
           sum(mean), sqrt(varmean), med, 0, 0)
  }
  stime <- x$time/scale
  if (is.numeric(rmean))
    rmean <- rmean/scale
  surv <- x$surv
  plab <- c("records", "n.max", "n.start", "events", "*rmean",
            "*se(rmean)", "median", paste(x$conf.int, c("LCL", "UCL"),
                                          sep = ""))
  ncols <- 9
  if (is.matrix(surv) && !is.matrix(x$n.event))
    x$n.event <- matrix(rep(x$n.event, ncol(surv)), ncol = ncol(surv))
  if (is.null(x$strata)) {
    if (rmean == "none")
      end.time <- NA
    else if (is.numeric(rmean))
      end.time <- rmean
    else end.time <- max(stime)
    if (is.matrix(surv)) {
      out <- matrix(0, ncol(surv), ncols)
      for (i in 1:ncol(surv)) {
        if (is.null(x$conf.int))
          out[i, ] <- pfun(x$n, stime, surv[, i], x$n.risk,
                           x$n.event[, i], NULL, NULL, start.time,
                           end.time)
        else out[i, ] <- pfun(x$n, stime, surv[, i],
                              x$n.risk, x$n.event[, i], x$lower[, i], x$upper[,
                                                                              i], start.time, end.time)
      }
      dimnames(out) <- list(dimnames(surv)[[2]], plab)
    }
    else {
      out <- matrix(pfun(x$n, stime, surv, x$n.risk, x$n.event,
                         x$lower, x$upper, start.time, end.time), nrow = 1)
      dimnames(out) <- list(NULL, plab)
    }
  }
  else {
    nstrat <- length(x$strata)
    stemp <- rep(1:nstrat, x$strata)
    last.time <- (rev(stime))[match(1:nstrat, rev(stemp))]
    if (rmean == "none")
      end.time <- rep(NA, nstrat)
    else if (is.numeric(rmean))
      end.time <- rep(rmean, nstrat)
    else if (rmean == "common")
      end.time <- rep(median(last.time), nstrat)
    else end.time <- last.time
    if (is.matrix(surv)) {
      ns <- ncol(surv)
      out <- matrix(0, nstrat * ns, ncols)
      if (is.null(dimnames(surv)[[2]]))
        dimnames(out) <- list(rep(names(x$strata), ns),
                              plab)
      else {
        cname <- outer(names(x$strata), dimnames(surv)[[2]],
                       paste, sep = ", ")
        dimnames(out) <- list(c(cname), plab)
      }
      k <- 0
      for (j in 1:ns) {
        for (i in 1:nstrat) {
          who <- (stemp == i)
          k <- k + 1
          if (is.null(x$lower))
            out[k, ] <- pfun(x$n[i], stime[who], surv[who,
                                                      j], x$n.risk[who], x$n.event[who, j],
                             NULL, NULL, start.time, end.time[i])
          else out[k, ] <- pfun(x$n[i], stime[who],
                                surv[who, j], x$n.risk[who], x$n.event[who,
                                                                       j], x$lower[who, j], x$upper[who, j],
                                start.time, end.time[i])
        }
      }
    }
    else {
      out <- matrix(0, nstrat, ncols)
      dimnames(out) <- list(names(x$strata), plab)
      for (i in 1:nstrat) {
        who <- (stemp == i)
        if (is.null(x$lower))
          out[i, ] <- pfun(x$n[i], stime[who], surv[who],
                           x$n.risk[who], x$n.event[who], NULL, NULL,
                           start.time, end.time[i])
        else out[i, ] <- pfun(x$n[i], stime[who], surv[who],
                              x$n.risk[who], x$n.event[who], x$lower[who],
                              x$upper[who], start.time, end.time[i])
      }
    }
  }
  if (is.null(x$lower))
    out <- out[, 1:7, drop = F]
  if (rmean == "none")
    out <- out[, -(5:6), drop = F]
  list(matrix = out[, , drop = T], end.time = end.time)
}

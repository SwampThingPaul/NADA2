#' Compute an ECDF and Distribution Parameters for Censored Data
#'
#' @description Computes the empirical cumulative distribution function (ECDF) for censored data. Estimates parameters of the distribution, including the mean and quantiles.
#' @param y1 Concentrations plus detection limits for indicator formatted data.
#' @param y2 Censoring indicators (logical. 1 or `TRUE` = censored, 0 or FALSE = detected) for indicator formatted data.
#' @param conf The confidence coefficient for confidence intervals around the Kaplan-Meier mean and median. Default = 0.95.
#' @param qtls Probabilities for the quantiles to be estimated.  Defaults are (0.10, 0.25, 0.50, 0.75, 0.90).  You may add and/or substitute probabilities -- all must be between and not including 0 to 1.
#' @param Cdf Logical `TRUE`/`FALSE` indicator of whether to plot the empirical cumulative distribution function (cdf).
#' @param printstats Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is TRUE.
#' @param Ylab Optional input text in quotes to be used as the variable name on the ecdf plot.  The default is the name of the `y1` input variable.
#' @param plot.pos numeric scalar between 0 and 1 containing the value of the plotting position constant. The default value is `plot.pos=0.375`, the Blom plotting position
#' @param q.type an integer between 1 and 9 selecting one of the nine quantile algorithms detailed below to be used. See `stats::quantile` for more detail, default is set to 7.
#' @importFrom survival Surv survfit
#' @importFrom stats quantile
#' @return
#' If `printstats=TRUE`: Based on the provided `conf` value, Kaplan-Meier summary statistics (`mean`,`sd`,`median`), lower and upper confidence intervals around the mean and median value, sample size and percent of censored samples are returned. The specified quantile values are also printed and returned.
#'
#' If `Cdf=TRUE`: The ecdf of censored data is plotted.
#' @details Quantiles and parameters are estimated using the enparCensored and ecdfPlotCensored functions of the EnvStats package. This avoids a small bias in the mean produced by the NADA package's cenfit function, which uses the reverse Kaplan-Meier procedure, converting left-censored to right-censored data prior to computing the ecdf and mean. See Gillespie et al. for more discussion on the bias.
#'
#' All printed values will also be output to an object if saved.  Values are character because of the possibility of a `<1`, but if no `<` symbol can be converted to numeric using the `as.numeric(...)` command.  For data without censoring cfit will also return values.  In that case the values labeled "KM" are not both Kaplan-Meier results and standard arithmetic mean, t-interval CIs on the mean, and quantiles.  See ?quantiles for choosing the q.type; default q.type = 7.
#'
#' @export
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Gillespie, B.W., et al., 2010.  Estimating Population Distributions When Some Data Are Below a Limit of Detection by Using a Reverse Kaplan-Meier Estimator. Epidemiology 21, 564-570.
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

cfit <- function(y1, y2, conf=0.95, qtls = c(0.10, 0.25, 0.50, 0.75, 0.90), plot.pos = 0.375, q.type = 7, Cdf = TRUE, printstats = TRUE, Ylab = NULL)
{ N <- length(y2)
if (is.null(Ylab)) Ylab <- deparse(substitute(y1))
Conf <- 100*conf
if (isFALSE(0.50 %in% qtls)) {qtls <- append(0.50, qtls)}     # median always needed
k <- length(qtls)
names(N) <- " "
ecdf.qtls <- vector("numeric", k)
ecdf.ranks <- vector("integer", k)
ecdf.yvals <- vector("numeric", k)
y2 <- as.numeric(y2)
#  if(is.logical(y2) | prod(as.numeric(y2)) == 0)   # Some zeros in y2 == indicator format for y2.
#  reactivate previous line to allow a 2nd section for interval censoring in the future.
DLs <- y1[y2==1]
PctND <- signif(100*length(DLs)/N, 4)
overall.min <- min(y1)
detect.min <- min(y1[y2==0])
nd.max <- ifelse(is.na(y1[y2==1]), 0, max(y1[y2==1]))
yzero <- y1*(1-y2)   # low end is zero for nondetects

# if some values are nondetects
if (sum(y2) != 0)
{ ycen2 <- 1-y2  # indicator for right-censored Surv
flip.const <- max(y1)   # use the max for the constant
yhi.flip <- flip.const - y1
y.surv <- Surv(yhi.flip, ycen2)
y.out<- survfit(y.surv ~ 1, conf.int = conf, conf.type = "plain")

npar.out <- enparCensored(y1, y2, ci = TRUE, pivot.statistic = "t", ci.sample.size = N)
KMmean <- npar.out$parameters[1]
std.err <- npar.out$parameters[3]
KMsd <- signif(npar.out$parameters[2], 4)
# could  use enparCensored output $limits [1,2] for conf limits on mean instead of computation below.  Get same numbers
qt.ci <- qt(c((1-conf)/2, 1-(1-conf)/2), N-1)
LCLmean <- signif(KMmean + qt.ci[1]*std.err, 4)
UCLmean <- signif(KMmean + qt.ci[2]*std.err, 4)
KMmean <- signif(KMmean, 4)

# using EnvStats::ecdfPlotCensored to get quantiles.  Blom plotting position is default
ecdf.out <- ecdfPlotCensored(y1, y2, plot.it = Cdf, plot.pos.con = plot.pos, type = "s", xlab = Ylab, main = "Empirical CDF of Censored Data", ecdf.lwd = par("cex"))
for (i in 1:k)  {ecdf.qtls[i] <- min(ecdf.out$Cumulative.Probabilities[ecdf.out$Cumulative.Probabilities >= qtls[i]])
ecdf.ranks[i] <- min(which(ecdf.out$Cumulative.Probabilities == ecdf.qtls[i]))
ecdf.yvals[i] <- ecdf.out$Order.Statistics[ecdf.ranks[i]] }

Qtls <- as.character(signif(ecdf.yvals, 4))
# change values below lowest detect or detection limit to < lowest limit
Qtls <- ifelse(ecdf.out$Censored[ecdf.ranks], paste("<", Qtls, sep=""), Qtls)
KMmedian <- Qtls[qtls == 0.50]

LCLmedian <- signif(flip.const - NADA2.survmean(y.out,  rmean=flip.const) [[1]]["0.95UCL"], 4)
LCLmedian <- ifelse(is.na(LCLmedian), KMmedian, LCLmedian)
UCLmedian <- signif(flip.const - NADA2.survmean(y.out,  rmean=flip.const) [[1]]["0.95LCL"], 4)
UCLmedian <- ifelse(is.na(UCLmedian), paste("<", UCLmedian, sep=""), UCLmedian)
}

# when all data are detects
else  {KMmean <- mean(y1, na.rm = TRUE)
KMsd <- sd(y1, na.rm = TRUE)
std.err <- KMsd/sqrt(length(na.omit(y1)))
qt.ci <- qt(c((1-conf)/2, 1-(1-conf)/2), N-1)
LCLmean <- signif(KMmean + qt.ci[1]*std.err, 4)
UCLmean <- signif(KMmean + qt.ci[2]*std.err, 4)
KMmean <- signif(KMmean, 4)
KMsd <- signif(KMsd, 4)
Qtls <- signif(quantile(y1, probs = qtls, na.rm = TRUE, type = q.type), 4)
KMmedian <-  signif(quantile(y1, probs = 0.50, na.rm = TRUE, type = q.type), 4)
ecdf.out <- ecdfPlot(y1, plot.it = Cdf, plot.pos.con = plot.pos, type = "s", xlab = Ylab, main = paste("Empirical CDF of", deparse(substitute(y1))), ecdf.lwd = par("cex"))

# getting CIs for median
orderstat.CI <- quantile_CI (length(na.omit(y1)), 0.50, alpha = 1-conf)
order.y <- y1[sort.list(y1)]
CI.median <- order.y[orderstat.CI$Interval]
LCLmedian <- signif(CI.median[1], 4)
UCLmedian <- signif(CI.median[2], 4)

# end of when all data are detects
}

#print results
quant.char <- paste ("Q", as.character(100*(qtls[1:k])), sep = "")
Cpct <- paste(Conf, "%", sep = "")
stats <- data.frame (N, PctND, Conf, KMmean, KMsd, KMmedian, LCLmean, UCLmean, LCLmedian, UCLmedian)

if (printstats == TRUE) {cat("\n", "Output for", Ylab, "            ", Cpct, "Confidence Intervals", "\n", "Statistics:", "\n")
  print(stats)
  cat("\n")
  cat("Quantiles:", quant.char, "\n", sep="\t")
  cat("        ", Qtls, "\n", "\n", sep="\t")
}

names(Qtls) <- quant.char
Qtls <- t(Qtls)
cstats <- cbind(stats, Qtls)
return (invisible (cstats))
}



#' quantile confidence interval
#'
#' @param n length
#' @param q quantile
#' @param alpha CI value
#' @keywords internal
#' @importFrom stats qbinom

#' @export
#'

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
  #
  # Return the order statistics and the actual coverage.
  #
  u <- rep(u, each=5)[i]
  l <- rep(l, 5)[i]
  return(list(Interval=c(l,u), Coverage=coverage[i]))
  # end of quantile_CI function
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

#' Seasonal Kendall permutation test on censored data
#'
#' @param time Column of the time variable, either a sequence of days or decimal times, etc.  Will be the scale used for time in the trend analysis.
#' @param y The column of y (response variable) values plus detection limits
#' @param y.cen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param group Column of the season classifications. A factor in R, so usually though not necessarily a text variable.  If numeric, define as a factor before running the script.
#' @param LOG Indicator of whether to compute the trend analysis in the original y units, or on their logarithms.  The default is to use the logarithms (LOG = `TRUE`).  To compute in original units, specify the option LOG = `FALSE` (or LOG = 0).
#' @param R The number of repetitions in the permutation process.  R is often between 999 and 9999 (adding +1 to represent the observed test statistic produces 1000 to 10000 repetitions). By default R=4999. Increasing R results in lower variation in the p-values produced between runs.
#' @param nmin The minimum number of observations needed for the entire time period to be tested, per season.  For example, with 1 sample per year per season over an 8-year period, you have 8 observations for each season.  You can increase this number if you want a higher minimum.  Don’t decrease it below the default of 4.  If there are fewer than nmin values that season is skipped and not included in the overall test and a note of this will be printed in the console.
#' @param seaplots In addition to the plot of the overall Seasonal Kendall trend line, plots of the trend in individual seasons can also be drawn.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param ... other inputs associated with modifying plots produced by this function.
#'
#' @return Prints the Kendall trend test results for each season individually. The overall Seasonal Kendall test and Theil-Sen line results are both printed and returned.
#' @details  For each season the ATS function is used to compute the season's Kendall-S statistic and p-value for the trend test.  The test is the usual ATS procedure, not a permutation test.  For the overall test there are R (default 4999) random rearrangements of data that are generated with no mixing of data from one season to another season within a permutation -- data over time within a season are randomized.  This retains any between-seasons differences while removing any trend from year to year to use as the permuted "representation of the null hypothesis" of no trend in the R Seasonal Kendall tests. For a 2-sided trend test the p-value is the number of the absolute value of the permutation S statistics that are equal to or greater than the absolute value of the observed S from the original data, plus 1 (for the observed S of the original data), divided by the total (R+1) S values.
#'
#' The censeaken function differs from other R packages that do not incorporate censored data in their trend tests (EnvStats, Kendall, rkt and others) in that censeaken uses all of the data across the years without an option to equalize each year's or season's influence by using the overall mean or median of the period's data rather than the original observations. Seasons with more data will have more influence on the outcome, just as years with more data will have more influence. If the numbers of observations differ enough between seasons or years that this is of concern, it is up to the user to perhaps eliminate some data in data-rich periods that are most unlike the conditions or times of data collected in sparse periods.  Or to compute summary statistics such as the median of each season/year combination and run the test on the medians, though this will result in a loss of power to detect trends and will require the user to use methods such as those in NADA2 to compute medians when there are censored data.
#'
#' If `seaplots=TRUE` each season's trend line will be plotted seperately. A plot of the overall Seasonal Kendall (Akritas-Theil-Sen) line is always plotted.
#' If `seaplots=FALSE` only the overall Seasonal Kendall (Akritas-Theil-Sen) line will be plotted on a data scatterplot.
#' @export
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Hirsch, R.M., Slack, J.R., Smith, R.A., 1982. Techniques of Trend Analysis for Monthly Water Quality Data, Water Res. Reseach 18, 107-121.

#'
#' @examples
#'\donttest{
#' data(Brumbaugh)
#'
#' # Artificial time and season variables for demonstration purposes
#' Brumbaugh$time=1:nrow(Brumbaugh)
#' Brumbaugh$sea=as.factor(round(runif(nrow(Brumbaugh),1,4),0))
#'
#' with(Brumbaugh,censeaken(time,Hg,HgCen,sea,seaplots = TRUE))
#'
#' censeaken(time,Hg,HgCen,sea,Brumbaugh,seaplots = TRUE)
#' }

censeaken <- function(time, y, y.cen, group, data = NULL,
                      LOG = FALSE, R = 4999, nmin = 4,
                      seaplots = FALSE, printstat = TRUE, ...) {
  # Capture variable names
  xname <- deparse(substitute(time))
  yname <- deparse(substitute(y))
  grpname <- deparse(substitute(group))

  # Evaluate variables in data if provided
  if (!is.null(data)) {
    time <- eval(substitute(time), data, parent.frame())
    y <- eval(substitute(y), data, parent.frame())
    y.cen <- eval(substitute(y.cen), data, parent.frame())
    group <- eval(substitute(group), data, parent.frame())
  }

  # Handle log transform if needed
  if (LOG) {
    yname <- paste0("ln(", yname, ")")
    y <- log(y)
  }

  # Clean and structure data
  dat <- na.omit(data.frame(time = time, y = y, y.cen = y.cen, group = group))

  # Summary information
  if (printstat) {
    cat("\nDATA ANALYZED:", yname, "vs", xname, "by", grpname, "\n")
  }

  # Prepare group splits
  yc <- split(dat$y, dat$group)
  xc <- split(dat$time, dat$group)
  cc <- split(dat$y.cen, dat$group)
  zc <- split(dat$group, dat$group)

  perm.sea <- matrix(0, nrow = R, ncol = length(yc))
  s_all <- 0
  denomall <- 0
  j <- 0  # valid group counter

  # Per-season computations
  for (i in seq_along(yc)) {
    sea <- zc[[i]][1]
    ntest <- length(xc[[i]])

    if (ntest < nmin) {
      warning(paste0(sea, "\nSeason dropped -- ", ntest, " are too few obs\n"))
      next
    }

    if (printstat) cat("----------\n")

    j <- j + 1
    perm.sea[1:R, j] <- computeS(xc[[i]], yc[[i]], cc[[i]], seas = sea, R = R)
    ats <- ATSmini(yc[[i]], cc[[i]], xc[[i]])
    medslop <- ats[2]
    int <- ats[1]
    s <- ats[[5]]
    tau <- signif(ats[3], 3)
    pval <- ats[4]
    denom <- ntest * (ntest - 1) / 2
    s_all <- s_all + s
    denomall <- denomall + denom

    if (seaplots) {
      x.seas <- range(xc[[i]])
      y.seas <- x.seas * as.numeric(medslop) + as.numeric(int)
      kenplot(yc[[i]], cc[[i]], xc[[i]], xcen = rep(0, ntest),
              xnam = xname, ynam = yname, atsline = TRUE, ...)
      lines(x.seas, y.seas, col = "purple")
      mtext(paste("Season =", sea))
    }

    if (printstat) {
      print(data.frame(
        Season = sea, N = ntest, S = s, Tau = tau, Pvalue = signif(pval, 5),
        Intercept = signif(int, 5), MedianSlope = signif(medslop, 4)
      ))
    }
  }

  if (printstat) cat("----------\n")

  # Overall stats
  Kendall_S <- rowSums(perm.sea, na.rm = TRUE)
  tau_all <- s_all / denomall
  ats_all <- ATSmini(dat$y, dat$y.cen, dat$time)
  medslope <- ats_all$slope
  intall <- ats_all$intercept

  K <- sum(abs(Kendall_S) >= abs(s_all))
  pval_all <- (1 + K) / (1 + length(Kendall_S))

  hist(Kendall_S, main = "Kendall's S statistic Permutation Test")
  abline(v = c(s_all, -s_all), col = "red", lty = 2)

  RESULTS <- data.frame(
    reps_R = R, N = nrow(dat), S_SK = s_all,
    tau_SK = signif(tau_all, 3), pval = signif(pval_all, 5),
    intercept = signif(intall, 5), slope = signif(medslope, 4)
  )

  if (printstat) {
    print(RESULTS)
    cat("----------\n")
  }

  kenplot(dat$y, dat$y.cen, dat$time, xcen = rep(0, nrow(dat)),
          xnam = xname, ynam = yname, Title = "Seasonal Kendall Test", ...)
  abline(intall, medslope, lwd = 2, col = "blue")
  mtext("Overall Trend Line", col = "blue")

  invisible(RESULTS)
}

#' Trend analysis of censored data with a covariate and seasonal blocks
#'
#' @description Computes the Seasonal Kendall trend test for censored data after adjustment of censored data for a covariate. The adjustment is by subtracting off a censored GAM smooth, removing the effect of the covariate.  Trend analysis is performed on the residuals from the GAM smooth.
#' @param y.var The column of y (response variable) values plus detection limits
#' @param y.cens The column of indicators, where 1 (or `TRUE`) indicates a detectionlimit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.var Column of a covariate (not time).  `y.var` will be smoothed versus `x.var` and residuals taken to subtract out the relationship between `y` and `x`.
#' @param time.var Column of the numerical time variable, either a sequence of numbered days or decimal times, etc.  Will be the scale used for time in ATS trend analysis.
#' @param R The number of repetitions in the permutation process.  R is often between 999 and 9999 (adding +1 to represent the observed test statistic produces 1000 to 10000 repetitions). By default R=4999. Increasing R results in lower variation in the p-values produced between runs.
#' @param nmin The minimum number of observations needed for the entire time period to be tested, per season.  For example, with 1 sample per year per season over an 8-year period, you have 8 observations for each season.  You can increase this number if you want a higher minimum.  Don’t decrease it below the default of 4.  If there are fewer than nmin values that season is skipped and not included in the overall test and a note of this will be printed in the console.
#' @param season  Column of the seasonal variable.  Usually a text variable but may be numeric.  Will be converted to a factor..
#' @param link Default = `"identity"` which means it uses data in the original units. See details.
#' @param Smooth Type of smoother used in the GAM. Default is `"cs"`, shrinkage cubic regression splines. See details for other options.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @param seaplots logical 'TRUE'/'FALSE' option to print plots with trend line for each season.  Default is 'TRUE'.
#' @keywords trend analysis GAM spline Seasonal Kendall
#' @export
#'
#' @importFrom mgcv gam
#' @importFrom cenGAM tobit1
#' @return
#'
#' Prints four plots: 1. Y data vs X covariate with GAM Smooth, 2. Residuals from GAM Smooth vs X covariate, 3. histogram of the SK test results illlustrating the p-value, and 4. Seasonal Kendall trend line of residuals vs time.  Plots of data and SK trend line for each season are produced when the seaplots option is TRUE.
#'
#' Returns the seasonal Kendall trend test results on residuals (intercept, slope, Kendall's tau, p-value for trend)
#'
#' @details
#'
#' Default `link` = identity. The y variables are then used in their original units. Other options are available see `cenGAM::tobit1` for more options.
#'
#' Default `Smooth` is `"cs"` for shrinkage cubic regression splines. See `mgcv::smooth.terms` for other types of smoothing algorithms.  '"ts"' is a thin-plate regression spline and is also commonly used.
#'
#' As with the censeaken function, observations are not edited when there are more data in some seasons than others. Seasons with more data will have more influence on the overall SK test than seasons with fewer data. To avoid this (as done with some Seasonal Kendall software) data in the seasons with more can be selectively deleted to better match the data frequency of the seasons with fewer data.

#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#'
#'@seealso [mgcv::gam]
#'
#' @examples
#'
#' data(Gales_Creek)
#' with(Gales_Creek,centrendsea(TCr,CrND,discharge,dectime,Season))

centrendsea <- function(y.var, y.cens, x.var, time.var, season, R = 4999, nmin = 4, link = "identity", Smooth = "cs", printstat=TRUE, seaplots = TRUE) {
  yname <- deparse(substitute(y.var))
  xname <- deparse(substitute(x.var))
  tname <- deparse(substitute(time.var))
  seasname <- deparse(substitute(season))
  resi.txt <- paste(yname, "residuals")

  y.cens <- as.logical(y.cens)
  y.var <- as.numeric(y.var)
  season <- as.factor(season)
  DL <- y.cens * y.var * 1.0001

  dat.all <- data.frame(y.var, x.var, time.var, DL, y.cens, season)
  dat.nonas <- na.omit(dat.all)
  colnames(dat.nonas) <- c(yname, xname, tname, "DL", "Cens", seasname)

  # cs, ts, ds  top 3, but ds gives bends that may be overfitting.
  gam.y <- gam(dat.nonas[,1] ~ s(dat.nonas[,2], bs = Smooth), family = tobit1(link = link, left.threshold = dat.nonas[,4]))
  # sorting fitted values
  o <- order(dat.nonas[,2] , dat.nonas [,1])

  # plot 1.  Y data vs covariate, with smooth
  x.cens = rep(0, times=length(dat.nonas[,3]) )
  cenxyplot(dat.nonas[,2], as.logical(x.cens), dat.nonas[,1], as.logical(dat.nonas[,5]), main = "1. Data and GAM Smooth", ylab = yname, xlab = xname, pch = 19, cex = 0.7)
  # draw the smooth
  lines (dat.nonas[o,2], gam.y$fitted.values[o], col = 'red')

  # plot 2. Y residuals from smooth vs. covariate
  plot(gam.y$residuals ~ dat.nonas[,2], main = "2. Residuals from GAM Smooth", xlab = xname, ylab = resi.txt)
  abline (h=0, col = "blue")

  # start censeaken on residuals
  if(printstat==TRUE){cat("\n", "Trend analysis by", seasname,"of:", yname, "adjusted for", xname, sep=" ","\n")}

  # bump up residuals so all are positive.  Changes intercept.
  bump.int <- min(gam.y$residuals)*(-1) + 1
  ttt<-dat.nonas[,3]
  yyy<-gam.y$residuals + bump.int
  ccc<- dat.nonas[,5]
  nall<-length(dat.nonas[,1])
  denom<-0;  denomall <- 0
  s_all <- 0;

  # Stats Not Needed for Trend Analysis
  # compute median of uncensored time (all data)
  # xmedian<-median(ttt)
  #  compute KM median of censored y.  Assumes all <ND go to 0 at lower end.
  ## y.dist <- cfit(yyy, ccc, Cdf = FALSE, printstat = FALSE, Ylab = resi.txt)
  ## ymedian <- y.dist$KMmedian
  # compute the Kaplan-Meier mean (all data)
  ## ymean <- y.dist$KMmean

  y.resi <- split(gam.y$residuals, dat.nonas[,6])
  yc<-split(yyy,dat.nonas[,6])
  xc<-split (ttt,dat.nonas[,6])
  cc<-split (ccc,dat.nonas[,6])
  zc<-split (dat.nonas[,6],dat.nonas[,6])
  # yc are the bumped residuals split into groups.  yc[i] are bumped residuals in the ith group.

  perm.sea <- matrix(0, nrow = R, ncol = length(yc))
  allslope<-rep(c(0),1)
  dsh=c("----------", "\n")

  #  EJG 7-1-2017 CHECK FOR NUMBER OF OBSERVATIONS IN EACH SEASON
  # splitting data season by season
  j=0        # j is number of seasons with sufficient data
  for(i in 1:length(yc))    # i is the number of groups/seasons
  {grp<-(zc[[i]])    # grp is the vector of group names for the ith group
  sea<-(grp[c(1)])   # sea is a single group name
  ntest<-length(xc[[i]])
  if (ntest < nmin) {
    cat(dsh)
    warning(paste(as.character(grp[1]),"\n","Season dropped --",ntest, "are too few obs","\n", sep=" "))
  }
  ## DONE CHECKING
  if (ntest < nmin) next

  if(printstat==TRUE){cat(dsh)}
  j=j+1
  #  Season by season computations
  perm.sea[1:R,j] <- computeS (xc[[i]], yc[[i]], cc[[i]], seas = sea, R = R)

  # ATS results for observed seasonal data
    ats.seas <- ATSmini(yc[[i]], cc[[i]], xc[[i]])

  medslop <- ats.seas[2]
  int <- ats.seas[1] - bump.int
  s <- ats.seas[[5]]
  tau <- signif(ats.seas[3], 3)
  pval <- ats.seas[4]
  denom <- (ntest*(ntest-1)/2)
  s_all <- s_all + s
  denomall <- denomall + denom

  # optional plots for each season
  if (seaplots == TRUE)  {
    x.seas <- c(min(xc[[i]]),max(xc[[i]]))
    y.seas <- 0
    slp=as.numeric(as.character(medslop[1]))
    int <- as.numeric(as.character(int[1]))
    y.seas <- x.seas*slp+int
    z <- data.frame(x.seas, y.seas)
    kenplot(y.resi[[i]], cc[[i]], xc[[i]], xcen = rep(0, times=ntest), xnam=tname, ynam=resi.txt, atsline = TRUE)
    lines(z, col = "purple")
    mtext(paste("Season =", sea))
  }

  # ATS results for the season, not pemutation tests
  RESULTS1<-data.frame(Season = sea, N = ntest, S=s, Tau=tau, Pvalue=signif(pval,5), Intercept = signif(int, 5), MedianSlope = signif(medslop, 4))
  if(printstat==TRUE){print(RESULTS1)}

  }         # end of seasonal computations
  if(printstat==TRUE){cat(dsh)}

  # Overall computations for all of the data
  if(printstat==TRUE){cat("Seasonal Kendall test and Akritas-Theil-Sen line on residuals", "\n");}
  Kendall_S <- as.vector(rowSums(perm.sea))   # S overall for each permutation
  tau_all<- (s_all)/denomall
  all.out <- ATSmini(yyy, ccc, ttt)
  medslope <- all.out$slope
  intall <- all.out$intercept - bump.int
  K <- 0

  # compute permutation two-sided p-value
  K <- length(Kendall_S[abs(Kendall_S) >= abs(s_all)])
  pval_all <- (1+K)/(1+length(Kendall_S))   # Adding 1 is due to the observed value from data

  hist(Kendall_S, main = "Kendall's S statistic Permutation Test")
  abline (v = s_all, col = "red", lty = 2)
  s_alt <- (-1)*s_all
  abline (v = s_alt, col = "red", lty = 2)

  RESULTS <- data.frame(reps_R = R, N=nall, S_SK=s_all, tau_SK =signif(tau_all,3), pval=signif(pval_all, 5), intercept = signif(intall,5), slope = signif(medslope,4))
  if(printstat==TRUE){print(RESULTS)}
  if(printstat==TRUE){cat(dsh)}

  kenplot(gam.y$residuals, ccc, ttt, xcen = rep(0, times=nall), xnam=tname, ynam=resi.txt, Title = "Overall Seasonal Kendall Test")
  abline(intall, medslope, lwd=2, col = "blue")
  mtext ("Overall Trend Line", col = "blue")

  return (invisible(RESULTS))
}


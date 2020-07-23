

#' Seasonal Kendall permutation test on censored data
#'
#' @param time Column of the time variable, either a sequence of days or decimal times, etc.  Will be the scale used for time in the trend analysis.
#' @param y The column of y (response variable) values plus detection limits
#' @param y.cen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param group Column of the season classifications. A factor in R, so usually though not necessarily a text variable.  If numeric, define as a factor before running the script.
#' @param LOG Indicator of whether to compute the regression in the original y units, or on their logarithms.  The default is to use the logarithms (LOG = `TRUE`).  To compute in original units, specify the option LOG = `FALSE` (or LOG = 0).
#' @param R The number of repetitions in the permutation process.  R is often between 999 and 9999 (+ the 1 observed test statistic produces 1000 to 10000 repetitions). By default R=4999. Increasing R simply results in lower variation in the pvalues produced between runs.
#' @param nmin The minimum number of observations needed for the entire time period to be tested, per season.  For example, with 1 sample per year per season over an 8-year period, you have 8 observations for each season.  You can increase this number if you want a higher minimum.  Don’t decrease it below 4.  If there are fewer than nmin values that season is skipped and not included in the overall test & a note will be printed.
#' @param seaplots In addition to the plot of the overall Seasonal Kendall line, plots for the individual seasons can be drawn.
#'
#' @return Prints Kendall Trend test results for each season (individually) and a Seasonal Kendall test and Theil-Sen line results.
#'
#' If `seaplots=TRUE` each season's trend line will be plotted along with overall Seasonal Kendall and Theil-Sen line.
#'
#' If `seaplots=FALSE` the entire dataset will be plotted with Thiel-Sen line overlayed.
#' @export
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' @seealso [NADA::cenken]
#'
#' @examples
#' library(NADA) #For example data
#'
#' data(HgFish)
#'
#' #Artifical time and season variables for demonstration purposes
#' HgFish$time=1:nrow(HgFish)
#' HgFish$sea=as.factor(round(runif(nrow(HgFish),1,4),0))
#'
#'
#' with(HgFish,censeaken(time,Hg,HgCen,sea,seaplots = T))


censeaken <- function(time, y, y.cen, group, LOG = FALSE,
                      R = 4999, nmin = 4, seaplots = FALSE)
{
  xname = deparse(substitute(time))
  yname = deparse(substitute(y))
  grpname = deparse(substitute(group))  # season column name

  if (LOG == TRUE)  { yname <- paste("ln(", yname, ")", sep = "")
  nonas <- na.omit(data.frame(time, log(y), y.cen, group)) }

  else {nonas <- na.omit(data.frame(time, y, y.cen, group)) }

  df = data.frame(TIME = xname, Y = yname, SEASON = grpname)
  cat("\n", "DATA ANALYZED:", yname, "vs", xname, "by", grpname, sep=" ","\n")

  xxx<-nonas$time
  yyy<-nonas[,2]
  ccc<- nonas$y.cen
  nall<-length(nonas[,1])
  denom<-0;  denomall <- 0
  s_all <- 0;

  # compute median of uncensored time (all data)
  xmedian<-median(xxx)
  #  compute KM median of censored y.  Assumes all <ND go to 0 at lower end.
  y.dist <- cfit(yyy, ccc, Cdf = FALSE, printstats = FALSE, Ylab = yname)
  ymedian <- y.dist$KMmedian

  # compute the Kaplan-Meier mean (all data)
  ymean <- y.dist$KMmean

  yc<-split(nonas[,2],nonas[,4])
  xc<-split (nonas[,1],nonas[,4])
  cc<-split (nonas[,3],nonas[,4])
  zc<-split (nonas[,4],nonas[,4])
  # yc are the y data split into groups.  yc[i] are data in the ith group.

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
    cat(as.character(grp[1]),"\n")
    cat("Note: Season dropped --",ntest, "are too few obs", sep=" ","\n")
  }
  ## DONE CHECKING
  if (ntest < nmin) next

  cat(dsh)
  j=j+1
  #  Season by season computations
  perm.sea[1:R,j] <- computeS (xc[[i]], yc[[i]], cc[[i]], seas = sea, R = R)

  # ATS results for observed seasonal data
  ats.seas <- ATSmini(yc[[i]], cc[[i]], xc[[i]])
  medslop <- ats.seas[2]
  int <- ats.seas[1]
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
    z=data.frame(x.seas, y.seas)
    kenplot(yc[[i]], cc[[i]], xc[[i]], xcen = rep(0, times=ntest), xnam=xname, ynam=yname)
    lines(z, col = "purple")
    mtext(paste("Season =", sea))
  }

  # ATS results for the season, not pemutation tests
  RESULTS1<-data.frame(Season = sea, N = ntest, S=s, Tau=tau, Pvalue=signif(pval,5), Intercept = signif(int, 5), MedianSlope = signif(medslop, 4))
  print(RESULTS1)

  }         # end of seasonal computations
  cat(dsh)

  # Overall computations for all of the data
  cat("Seasonal Kendall test and Theil-Sen line", "\n");
  Kendall_S <- as.vector(rowSums(perm.sea))   # S overall for each permutation
  tau_all<- (s_all)/denomall
  all.out <- ATSmini(yyy, ccc, xxx)
  medslope <- all.out$slope
  intall <- all.out$intercept
  K <- 0

  # compute permutation two-sided p-value
  K <- length(Kendall_S[abs(Kendall_S) >= abs(s_all)])
  pval_all <- (1+K)/(1+length(Kendall_S))   # Adding 1 is due to the observed value from data

  hist(Kendall_S, main = "Kendall's S statistic Permutation Test")
  abline (v = s_all, col = "red", lty = 2)
  s_alt <- (-1)*s_all
  abline (v = s_alt, col = "red", lty = 2)

  RESULTS <- data.frame(reps_R = R, N=nall, S_SK=s_all, tau_SK =signif(tau_all,3), pval=signif(pval_all, 5), intercept = signif(intall,5), slope = signif(medslope,4))
  print(RESULTS)
  cat(dsh)

  kenplot(yyy, ccc, xxx, xcen = rep(0, times=nall), xnam=xname, ynam=yname)
  abline(intall, medslope, lwd=2, col = "blue")
  mtext ("Overall Line")
}

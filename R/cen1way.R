#' Peto-Peto one-factor test
#'
#' @description
#' Performs a Peto-Peto nonparametric test of differences in cdfs between groups.  If more than two groups, the test is followed by a nonparametric multiple comparison test.  Uses the BH method of adjusting p-values.
#' @param x1 The column of data values plus detection limits
#' @param x2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the y1 column, and 0 (or `FALSE`) indicates a detected value in y1.
#' @param group Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param mcomp.method One of the standard methods for adjusting p-values for multiple comparisons.  Type ?p.adjust for the list of possible methods. Default is Benjamini-Hochberg "BH" false discover rate.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#'
#' @importFrom survminer pairwise_survdiff
#' @importFrom survival survdiff Surv
#' @importFrom stats pchisq
#' @importFrom EnvStats enparCensored
#' @export
#' @return  A list of summary statistics for each group evaluated containing the following components:
#' \itemize{
#' \item `N` Number of samples
#' \item `PctND` Percentage of non-detects
#' \item `KMmean` Kaplan-Meier estimate of the mean
#' \item `KMsd` Kaplan-Meier estimate of standard deviation
#' \item `KMmedian` Kaplan-Meier estmate of the median
#' }
#'
#' Peto-Peto test results including Chi-Squared value, degrees of freedom and `p-value` of the test.
#'
#' If more than two groups, `p-values` of the pairwise multiple comparisons, adjusted using the BH false-discovery rate, are reported.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Peto, R., Peto, J., 1972. Asymptotically Efficient Rank Invariant Test Procedures. Journal of the Royal Statistical Society. Series A (General) 135, 185. \doi{https://doi.org/10.2307/2344317}
#'
#' Benjamini, Y., Hochberg, Y., 1995. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.  Journal of the Royal Statistical Society. Series B (Methodological), 57, 289-300.
#'

#'
#' @examples
#' data(PbHeron)
#'
#' # Two Groups
#' cen1way(PbHeron$Liver,PbHeron$LiverCen,PbHeron$DosageGroup)
#'
#' # More than two groups
#' cen1way(PbHeron$Liver,PbHeron$LiverCen,PbHeron$Group)


cen1way <- function(x1, x2, group, mcomp.method = "BH", printstat=TRUE) {

  yname <- deparse(substitute(x1))
  gname <- deparse(substitute(group))

  ydat <- na.omit(data.frame(x1, x2, group))
  y1 <- ydat[,1];  y2 <- as.logical(ydat[,2]); grp <- ydat[,3]

  rho=1
  fconst <- max(y1) + 1
  flip <- fconst - y1
  detect <- as.logical(1 - as.integer(y2))  # reverses TRUE/FALSE to fit survfit
  Factor <- as.factor(grp)
  df <- length(levels(Factor))-1
  CensData <- data.frame (flip, detect, Factor)
  Cen.stats <- matrix(0, nrow=df+1, ncol = 5)

  y.surv <- Surv(CensData$flip, CensData$detect, type="right")
  y.out<- survdiff(y.surv ~ CensData$Factor, rho=rho)
  pval = pchisq(y.out$chisq, df, lower.tail = FALSE)

  groupnames <- as.character(levels(Factor))
  for (i in 1:nlevels(Factor))    {
    y1gp <- y1[Factor==groupnames[i]]
    y2gp <- y2[Factor==groupnames[i]]
    Cstats <- suppressWarnings(cfit(y1gp, as.logical(y2gp), printstat=FALSE, Cdf = FALSE))
    Cstats <- Cstats[[1]][c(1:6)] ## summary stats
    Cstats <- Cstats[-3] ## N, PcfND, Mean, SD and LCLmean only
    Cstats <- data.frame(Cstats)
    Cstats$grp <-  groupnames[i] # added to include group in summary data frame.
    rownames(Cstats) <-NULL # added to clean up row names.
    if (i ==1) {Cen.stats <- Cstats;cnames <- colnames(Cstats)}else {colnames(Cstats)=cnames;Cen.stats <- suppressWarnings(rbind(Cen.stats, Cstats))}
    # added for custom warning.
    if(sum(as.logical(y2gp))==0){warning("One or more group(s) do not have censored data.",call.=F)}
  }

  # rownames(Cen.stats) <-  groupnames
  Cen.stats <- Cen.stats[,c(6,1:5)] # reordered the data frame
  colnames(Cen.stats) <- cnames[c(6,1:5)] # reordered the data frame

  if(printstat==TRUE){
    print.data.frame(Cen.stats, print.gap = 3)
    cat('\n',"     Oneway Peto-Peto test of CensData:", yname, "  by Factor:", gname, '\n', "     Chisq =", signif(y.out$chisq, 4), "  on", df, "degrees of freedom", "    p =", signif(pval,3), '\n')
    if (df >1) {mcomp <- pairwise_survdiff(Surv(flip, detect) ~ Factor, data=CensData, p.adjust.method = mcomp.method, rho=rho)
    print(mcomp)}
  }

  PetoPeto<-list(data.name=paste(yname, "  by Factor:", gname,"\n"),
                 ChiSq=signif(y.out$chisq, 4),df=df,p.value=signif(pval,3))

  if(df>1){mcomp=mcomp}else{mcomp=NA}

  x<-list(SummaryStats=Cen.stats,
          PetoPeto=PetoPeto,
          mcomp=mcomp)
  invisible(x)
}

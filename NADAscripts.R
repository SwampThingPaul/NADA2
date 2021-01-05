#' Permutation Analysis of Similarity (anosim) for Censored Data
#'
#' @description Plots the permutation histogram and test statistic produced by an anosim (nonparametric multivariate) test of differences between groups.
#' @param ano.out an `anosim` object. See details and example
#' @param hcol color of histogram
#' @param title title of histogram
#' @return Plots a histogram of the permutation test statistics representing the null hypothesis along with the observed test statistic from the data.  The p-value is the proportion of test statistics equal to or more extreme than the observed test statistic.
#' @export
#' @importFrom vegan anosim
#' @importFrom graphics hist
#' @seealso [vegan::anosim]
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Oksanen, J., Guillaume, F., 2018. Vegan: ecological diversity. CRAN R-Project.<https://cran.r-project.org/web/packages/vegan/index.html>
#'
#' @examples
#' data(PbHeron)
#'
#' # ROS model for each group
#' PbHeron.high <- with(subset(PbHeron,DosageGroup=="High"),NADA::ros(Blood,BloodCen))
#' PbHeron.high <- data.frame(PbHeron.high)
#' PbHeron.high$DosageGroup <- "High"
#'
#' PbHeron.low <- with(subset(PbHeron,DosageGroup=="Low"),NADA::ros(Blood,BloodCen))
#' PbHeron.low <- data.frame(PbHeron.low)
#' PbHeron.low$DosageGroup <- "Low"
#'
#' PbHeron.ros=rbind(PbHeron.high,PbHeron.low)
#'
#' # ANOSIM analysis
#' library(vegan)
#' PbHeron.anosim <- with(PbHeron.ros,anosim(modeled,DosageGroup))
#' summary(PbHeron.anosim)
#'
#' # Plot
#' anosimPlot(PbHeron.anosim)


anosimPlot <- function(ano.out, hcol = "light blue", title = "Histogram of anosim permutations") {
  min.p <- min(ano.out$perm)
  hset <- hist(ano.out$perm, plot = FALSE)
  yset <- max(hset$counts)
  label.txt <- paste("R =", signif(ano.out$statistic,2),"")
  xmax <- max(ano.out$perm, ano.out$statistic)

  hist(ano.out$perm, col = hcol, xlim = c(min.p, xmax), main = title, xlab = "Test statistics")
  abline (v = ano.out$statistic, lwd = 2)
  text (ano.out$statistic, yset, labels = label.txt, adj = c(1,1))
}

#' Atrazine concentrations in Nebraska ground water
#'
#' @description
#' From the `NADA` R-Package.
#'
#' Atrazine concentrations in a series of Nebraska wells before (June) and after (September) the growing season.
#'
#' Objective is to determine if concentrations increase from June to September. There is one detection limit, at 0.01 ug/L. Used in Chapters 4, 5, and 9 of the NADA book.
#'
#' @usage data(atrazine)
#'
#' @docType data
#' @keywords dataset
#' @name atrazine
#' @source Junk et al., 1980. Areal, vertical, and temporal differences in ground water chemistry: II. Journal of Environmental Quality. 9(3) 479 - 483.
#' @references
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.

"atrazine"

#' Akritas-Theil-Sen line for censored data
#'
#' @description
#' Computes Kendall's tau and the Akritas-Theil-Sen (ATS) line for censored data, along with the test that the slope (and Kendall's tau) equal zero.  For one x variable regression.
#' @param y.var The column of y (response variable) values plus detection limits
#' @param y.cen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and `0` (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.var The column of x (explanatory variable) values plus detection limits
#' @param x.cen The x-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `x.var` column, and `0` (or `FALSE`) indicates a detected value in `x.var`.
#' @param LOG Indicator of whether to compute the ATS line in the original y units, or for their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @param retrans Indicator of whether to retransform the plot and line back to original Y-variable units.  Not needed when `LOG = FALSE`. When `retrans = FALSE` & `LOG = TRUE` the plot is drawn with logY units (default). When `retrans =  TRUE` & `LOG = TRUE` the plot is drawn with original Y units.
#' @param xlabel Custom label for the x axis of plots.  Default is x variable column name.
#' @param ylabel Custom label for the y axis of plots.  Default is y variable column name.
#' @keywords trend censored
#' @export
#' @return  Coefficients (intercept and slope) for the ATS line are printed, along with Kendall's tau correlation coefficient, test statistic S, and the (single) p-value for the test that tau and slope both equal zero. A scatterplot with the fitted trend-line superimposed is also drawn.
#'
#' @importFrom graphics abline layout legend lines mtext par plot text
#' @importFrom utils data
#' @importFrom NADA cenken
#' @importFrom stats na.omit
#'
#' @references
#' Akritas, M.G., Murphy, S.A., LaValley, M.P., 1995. The Theil-Sen Estimator With Doubly Censored Data and Applications to Astronomy. Journal of the American Statistical Association 90, 170–177. https://doi.org/10.2307/2291140
#'
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @examples
#' # Both y and x are censored
#' data(PbHeron)
#' with(PbHeron, ATS(Blood, BloodCen, Kidney, KidneyCen))
#'
#' # x is not censored
#' data(Brumbaugh)
#' with(Brumbaugh,ATS(Hg, HgCen, PctWetland))

ATS <- function(y.var, y.cen, x.var, x.cen = rep(0, times=length(x.var)), LOG = TRUE, retrans = FALSE, xlabel = NULL, ylabel = NULL, Title = NULL)
  {
  yname <- deparse(substitute(y.var))
  xname <- deparse(substitute(x.var))
  y.cen <- as.logical(y.cen)
  x.cen <- as.logical(x.cen)
  alldat <-  data.frame(y.var, y.cen, x.var, x.cen)
  alldat <- na.omit(alldat)
  if (!is.null(xlabel)) {xname = xlabel}
  if (!is.null(ylabel)) {yname = ylabel}
  nobs <- length(alldat[,1])

  if (LOG == TRUE)  {
    if (min(alldat[,1]) < 0) stop('Cannot take logs of negative values.  Use LOG=FALSE')
    y.log <- log(alldat[,1])
    NoShift <- cenken(y.log, alldat[,2], alldat[,3], alldat[,4] )
    slope <- NoShift$slope
    intercept <- NoShift$intercept
    pval <- NoShift$p
    tau <- NoShift$tau
    S <- tau*nobs*(nobs-1)*0.5

    # if y.log goes negative, intercept = NA
    if (min(y.log) < 0) {
      shift.amt <- abs(min(y.log)) + 0.0001
      y.shift <- y.log + shift.amt
      Shifted<-cenken(y.shift, alldat[,2], alldat[,3], alldat[,4])
      intercept <- Shifted$intercept - shift.amt
    }
    ylogname <- paste("ln(", yname, ")", sep = "")
    cat ("Akritas-Theil-Sen line for censored data", "\n", "\n")
    if (slope >= 0) {cat (ylogname, "=", round(intercept,4), "+", round(slope,4), "*", xname, "\n")
      subtitle <- paste (ylogname, "=", round(intercept,4), "+", round(slope,4), "*", xname)
      short.t <- paste(round(intercept,4), "+", round(slope,4), "*", xname)
    }
    else {cat (ylogname, "=", round(intercept,4), round(slope,4), "*", xname, "\n")  # neg slope
      subtitle <- paste(ylogname, "=", round(intercept,4), round(slope,4), "*", xname)
      short.t <- paste(round(intercept,4), round(slope,4), "*", xname)
    }
    cat("Kendall's tau =", round(tau,4), "  p-value =", round(pval, 5), "\n", "\n")
    # draw a scatterplot
    x <- c(min(alldat[,3]), max(alldat[,3]))
    y <- x*slope+intercept
    z <- data.frame(x,y)
#    print(z)
    kenplot(log(alldat[,1]), alldat[,2], alldat[,3], alldat[,4], atsline = TRUE, xnam=xname, ynam=ylogname)
 #   lines(z, col = "purple")
    mtext(subtitle)
    #
    if (retrans == TRUE)
    { # draw a 2nd scatterplot in original units
      subtitle <- paste (yname, " = e^(", round(intercept,4), ")",  " * ", "e^(", round(slope,4), " * ", xname,")", sep="" )
      delta <- seq(1:100)
      diff<-max(alldat[,3]) - min(alldat[,3])
      xo = delta *diff/100 + min(alldat[,3])
      yo = exp(xo*slope)*exp(intercept)
      z <- data.frame(xo,yo)
      kenplot(alldat[,1], alldat[,2], alldat[,3], alldat[,4], atsline = FALSE, xnam=xname, ynam=yname)
      lines(xo, yo, col = "purple")
      mtext(subtitle)
      #
    }
    out <- data.frame(intercept, slope, S, tau, pval)
  }
  else  # no logarithms used for y
  {
    Noshift<-cenken(alldat[,1], alldat[,2], alldat[,3], alldat[,4])
    slope <- Noshift$slope
    intercept <- Noshift$intercept
    pval <- Noshift$p
    tau <- Noshift$tau
    S <- tau*nobs*(nobs-1)*0.5

    if (min(alldat[,1]) < 0) {  # y goes negative. Intercept is NA.
      shift.amt <- abs(min(alldat[,1])) + 0.0001
      y.shift <- alldat[,1] + shift.amt
      Shifted<-cenken(y.shift, alldat[,2], alldat[,3], alldat[,4])
      intercept <- Shifted$intercept - shift.amt
    }

    cat ("Akritas-Theil-Sen line for censored data", "\n", "\n")
    if (slope >= 0) {cat (yname, "=", round(intercept,4), "+", round(slope,4), "*", xname, "\n")
      subtitle <- paste (yname, "=", round(intercept,4), "+", round(slope,4), "*", xname)}
    else {cat (yname, "=", round(intercept,4), round(slope,4), "*", xname, "\n")
      subtitle <- paste (yname, "=", round(intercept,4), round(slope,4), "*", xname) }
    cat("Kendall's tau =", round(tau,4), "  p-value =", round(pval, 5), "\n", "\n")

    # draw a scatterplot with kenplot
    x <- c(min(alldat[,3]),max(alldat[,3]))
    y <- x*slope+intercept
    z=data.frame(x,y)
    kenplot(alldat[,1], alldat[,2], alldat[,3], alldat[,4], xnam=xname, ynam=yname)
    lines(z, col = "purple")
    mtext(subtitle)
    #
    out <- data.frame(intercept, slope, S, tau, pval)
  }
  return(invisible(out))
}

#' Kendall's tau and ATS line for censored data
#'
#' @description Computes Kendall's tau and the Akritas-Theil-Sen (ATS) line for censored data.  Is called by censeaken because it is much faster than the ATS function.  It is not expected to be of much use to users on its own. The x variable (time) may not be censored.
#' @param y.var The column of y (response variable) values plus detection limits.
#' @param y.cen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.var The column of x (time for a trend test) values.
#' @importFrom NADA cenken
#' @return Returns the intercept and slope (ATS line), tau (Kendall's tau), p-value and S-value (test statistic).
#' @export
#' @references
#' Akritas, M.G., Murphy, S.A., LaValley, M.P., 1995. The Theil-Sen Estimator With Doubly Censored Data and Applications to Astronomy. Journal of the American Statistical Association 90, 170–177. https://doi.org/10.2307/2291140
#'
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#' @seealso [NADA::cenken]
#' @examples
#' # x may not be censored.  Use the ATS function when x is censored.
#' data(Brumbaugh)
#'
#' with(Brumbaugh, ATSmini(Hg, HgCen, SedLOI))

ATSmini <- function(y.var, y.cen, x.var) {
  y.cen <- as.logical(y.cen)
  nobs <- length(y.var)

  Noshift<-cenken(y.var, y.cen, x.var)
  slope <- Noshift$slope
  intercept <- Noshift$intercept
  pval <- Noshift$p
  tau <- Noshift$tau
  S <- tau*nobs*(nobs-1)*0.5

  if (min(y.var < 0)) {
    # y goes negative. Intercept is NA.
    shift.amt <- abs(min(y.var)) + 0.0001
    y.shift <- y.var + shift.amt
    Shifted<-cenken(y.shift, y.cen, x.var)
    intercept <- Shifted$intercept - shift.amt
  }

  out <- data.frame(intercept, slope, tau, pval, S)
  return(out)
}

#' Cluster Matrix of Binary Censored Data
#'
#' @description Performs clustering of a matrix of 0s and 1s, ie. the censoring indicator columns for multiple variables. If there are multiple detection limits within a column first convert the 0/1 designated to be Below vs Above the highest detection limit in the column.  Detection limits may differ among columns.
#' @param dat.frame A data frame containing only the 0/1 columns, one column per chemical parameter.
#' @param method Method of forming clusters.  The default is `"ward.D2"`, which is appropriate for a variety of types of data.  Another appropriate option is `“average”` – average distances between cluster centers.  See the vegan package for other possible clustering methods.
#' @param group Optional grouping variable. If used, sites being clustered will be represented by their group name, rather than by the row number.
#' @param ncluster Optional number of clusters to be differentiated on the graph. Clusters are fenced off with rectangles.
#' @importFrom stats hclust cutree rect.hclust
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.

#'@export
#' @return Prints a cluster dendrogram based on clustering method
#'
#' @examples
#' data(PbHeron)
#'
#' # without group specified
#' binaryClust(PbHeron[,4:15])
#'
#' # With Group argument
#' binaryClust(PbHeron[,4:15],group=PbHeron$DosageGroup)

binaryClust <- function(dat.frame, method = "ward.D2", group = NULL, ncluster = NULL) {

  simple_matching_coeff <- binaryDiss(dat.frame)
  dat.clust = hclust(simple_matching_coeff, method = method)
  if (is.null(group))  {plot(dat.clust)}
  else { plot(dat.clust, label = group)}

  if (!is.null(ncluster))   {c5=cutree(dat.clust, k=ncluster)
  rect.hclust(dat.clust, ncluster)
  }
}

#' Binary dissimilarity coefficient matrix
#'
#' @description Computes a simple matching dissimilarity coefficient
#' @param dat.frame A data frame containing only the 0/1 columns.
#' @export
#' @importFrom vegan designdist
#' @return Returns a binary dissimilarity matrix.
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @seealso [vegan::designdist]
#' @examples
#' data(PbHeron)
#'
#' binaryDiss(PbHeron$LiverCen)

binaryDiss <- function(dat.frame) {
  dat.diss <- designdist(dat.frame, method = "1 - (a+d)/(a+b+c+d)", abcd = TRUE, terms = "binary", name = "simplematch")
  #  This returns a distance matrix.  0 = identical.
  #  for a similarity matrix, use binarySim
  return(dat.diss)
}

#' Plot Nonmetric Multidimensional Scaling of binary censored data
#'
#' @description Plots an NMDS of a matrix of 0s and 1s, the censoring indicator columns for multiple variables, to discern the pattern of data below vs. above the detection limit.  With multiple detection limits within a column, re-censoring to the highest limit in the column must be done prior to running this function.  May have different censoring levels in different columns.
#' @param dat.frame A data frame containing only the columns of 0/1 values.
#' @param group Optional grouping variable. Sites will be represented by different colored symbols for each group.
#' @param title Optional title for the NMDS graph.
#' @param legend.pos When group is specified, determines the location of the legend on the graph showing the colors representing each group’s data.  Default is “bottomleft”.  Alternatives are “topright” and “centerleft”, etc.
#' @importFrom vegan metaMDS ordiplot
#' @importFrom graphics points
#' @export
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @return Plots an NMDS of censored data represented as the binary Above vs Below a detection limit for each parameter.
#' @details Binary data may not provide sufficient information to discern differences in location on the plot if sample size is small.  Prior to runnning this analysis it is suggested to consult best analysis practise when performing NMDS. As a rule of thumb, an NMDS ordination with a stress value around or above 0.2 is deemed suspect and a stress value approaching 0.3 indicates that the ordination is arbitrary. Stress values equal to or below 0.1 are considered fair, while values equal to or below 0.05 indicate good fit.
#' @seealso [vegan::metaMDS]
#'
#' @examples
#' \dontrun{
#' data(PbHeron)
#'
#' # without group specified
#' binaryMDS(PbHeron[,4:15])
#'
#' # With Group argument
#' binaryMDS(PbHeron[,4:15],group=PbHeron$DosageGroup)
#' }



binaryMDS <- function(dat.frame, group = NULL, title=NULL, legend.pos = "bottomleft") {
  matname <- deparse(substitute(dat.frame))
  if (is.null(title)) {title <- paste("NMDS of", matname)}
  rownumb <- nrow(dat.frame)
  row.labels <- as.character(1:rownumb)
  dat.dist <- binaryDiss(dat.frame)  # distance matrix.  0 is identical.
  dat.mds <- metaMDS(dat.dist, zerodist = "add", autotransform = FALSE)
  plot(dat.mds, type = "n", main = title)
  if (is.null(group)) {text(dat.mds$points,labels=row.labels)}
  else {
    gp.factor <- as.factor(group)
    gp.nlevels <- c(1:nlevels(gp.factor))
    gp.col <- as.integer(gp.factor)
    gp.symb <- 18+gp.col
    plot.col=ordiplot(dat.mds, type="none", display="sites", main=title)
    #    print(gp.symb)
    #    points(plot.col, "sites", pch=gp.symb, col=gp.col)
    points(plot.col, "sites", pch=19, col=gp.col)
    text(plot.col, "sites", pos=4, cex = 0.8)
    #    legend(legend.pos, legend=levels(gp.factor), bty="n",col = gp.nlevels, pch = gp.symb)
    legend(legend.pos, legend=levels(gp.factor), bty="n",col = gp.nlevels, pch = 19)
  }

}

#' Binary similarity coefficient matrix
#'
#' @description Computes a simple matching similarity coefficient
#' @param dat.frame A data frame containing only the 0/1 columns.
#' @export
#' @importFrom vegan designdist
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @return Returns a binary similarity matrix.
#' @seealso [vegan::designdist]
#' @examples
#' library(NADA) #For example data
#' data(Golden)
#'
#' binarySim(Golden$LiverCen)
#'

binarySim <- function(dat.frame) {
  dat.symm <- designdist(dat.frame, method = "(a+d)/(a+b+c+d)", abcd = TRUE, terms = "binary", name = "simplematch")
  #  This returns a similarity matrix.  0 = disjoint.
  #  for a dissimilarity matrix, use binaryDiss
  return(dat.symm)
}

#' Brumbaugh
#'
#' @description
#' From the `NADA` R-Package.
#'
#' Mercury concentrations in fish across the United States.
#'
#' Objective is to determine if mercury concentrations differ by watershed land use. Can concentrations be related to water and sediment characteristics of the streams?
#'
#' There are three detection limits, at 0.03, 0.05, and 0.10 ug/g wet weight. Used in Chapters 10, 11 and 12 of the NADA book.
#'
#' @usage data(Brumbaugh)
#'
#' @docType data
#' @keywords dataset
#' @name Brumbaugh
#' @source Brumbaugh et al., 2001, USGS Biological Science Report BSR-2001-0009.
#' @references
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.

"Brumbaugh"

#' Draws censored boxplots
#'
#' @description
#' Portions below the maximum detection limit(s) are not shown by default, as their percentiles are not known.
#' @param y1 The column of y (response variable) values plus detection limits.
#' @param y2 The y-variable censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the y1 column, and 0 (or `FALSE`) indicates a detected value in y1.
#' @param group An optional column of a grouping variable.  Draws side-by-side boxplots if this variable is present.
#' @param LOG `TRUE`/`FALSE` indicator of whether to plot the Y axis data on the original scale (`FALSE`) or log scale (`TRUE`).
#' @param show `TRUE`\/`FALSE` indicator of whether to show estimated values computed using ROS for the portion of the box below the maximum DL (`TRUE`), or just leave the lower portion blank (`FALSE`).
#' @param minmax `NULL`/`FALSE` indicator of whether to draw outliers individually. Default is to show outliers. Setting `minmax = FALSE` (or any text) will draw the whiskers out to the max and min of the dataset.
#' @param ordr A vector indicating the order of boxes to be drawn on the boxplot, if not in alphabetical order (the default).  Example: for 4 boxplots for groups A, B, C, D, to change the order to the reverse type ordr = c(4, 3, 2, 1).  Example 2: To change the order to A, C, D, B, type ordr = c(1, 3, 4, 2)
#' @param Ylab Y axis label text, if something is wanted other than the Y variable name in the dataset.
#' @param Xlab X axis label text, if something is wanted other than the group variable name in the dataset.
#' @param Title Text to show as the graph title.  Default is blank.
#' @param dl.loc Location indicator of where to plot the "MaxDL=" text on some versions of the plot.  Possible entries are “topleft”, “topright”, “topcenter”, and the corresponding “bottom” text.
#' @param dl.col Color of the max detection limit line(s), and the legend text stating the max DL.  Default is “red”, but all standard R colors may be used.
#' @param bxcol Color for interior of boxplots. Specify just one color if all boxes are to be the same color.  If a different color is desired for each of three boxplots, as one example, use bxcol = c(“red”, “white”, “blue”) etc.
#' @param Ymax Maximum Y value to be shown on the plot.  Used to cut off high outliers on plot and better show the bulk of the boxplots.
#' @details If maximum detection limits vary among groups, separate maxDL lines will be drawn for each group's boxplot. If one group has fewer than 3 detected observations its boxplot will not be drawn.  Its detection limits will not count when computing the maximum limit.  However, if only one boxplot is drawn for the entire dataset by not specifying a group variable, the detection limits from the portion that is the mostly ND group will be used when computing the maximum limit.
#' @export
#' @importFrom graphics boxplot lines plot polygon
#' @importFrom grDevices adjustcolor
#' @importFrom NADA ros
#' @import utils
#'
#' @return Prints a boxplot with detection limit identified and a concatenated list of the maximum detection limit for each group.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @examples
#' data(PbHeron)
#' cboxplot(PbHeron$Liver,PbHeron$LiverCen,PbHeron$Group)



cboxplot <- function(y1, y2, group=NULL, LOG =FALSE, show=FALSE, ordr = NULL, Ylab=yname, Xlab = gname, Title = NULL, dl.loc = "topright", dl.col = "red", bxcol = "white", Ymax = NULL, minmax = NULL) {

  oldpar<- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  box.fill <- adjustcolor( "white", alpha.f = 0.6)
  if (show==TRUE) {bdl.col <- box.fill}
  else {bdl.col <- "white"}
  yname <- deparse(substitute(y1))
  if (is.null(Ylab))  Ylab <- yname
  y1 <- as.numeric(y1)
  y2 <- as.integer(y2)
  xmin = 0
  Ylim <- NULL
  grp.all = "1"
  gname = NULL

  if (sum(y2) > 0)    # not all data are detects
  { dlmax <- max(y1[y2 == 1])
  dltxt <- paste("Max DL=", signif(dlmax, 5), sep="")
  nonas <- na.omit(data.frame(y1, y2))         # omits NAs for the sake of the ros function
  y.nona <- nonas[,1]
  nd <- (nonas[,2] == 1)   # nd is vector of logical T/F
  if(sum(1-nonas[,2]) <= 2 )  {
    cat("Note: Data had fewer than 3 detects and so cannot be plotted", "\n")
    stop }
  try(if(sum(as.integer(nonas[,2])) == 0) stop("All data are detected.  Use the boxplot function."))

  if (LOG == TRUE)  {
    if (is.null(group) == TRUE) {    # log scale, no group
      xx = c(0.5, 1.5, 1.5, 0.5)
      y.ros <- suppressWarnings(ros(y.nona, nd))
      y.min <- log(min(y.ros$modeled))
      yy = c( y.min, y.min, log(dlmax), log(dlmax))
      Ylab = c(Ylab, "(natural logs)" )
      if (is.null(Ymax) == FALSE) { Ylim = c(y.min, log(Ymax))}
      if (is.null(minmax)) {
        # boxplot with separate outliers beyond 1.5*IQR
        boxplot(log(y.ros$modeled), ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplecol = "white", outcex = 0.8, main =Title, col=bxcol, ylim = Ylim)}
      else {
        #  boxplot drawn to max and min
        boxplot(log(y.ros$modeled), ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplewex = 0.2, range=0,  main =Title, col=bxcol, ylim = Ylim)}

      polygon (xx, yy, col = bdl.col, border=bdl.col)
      abline (h=y.min, col=bdl.col, lwd=8)
      abline(h=log(dlmax), col = dl.col, lty="longdash", lwd=2)
      text(1.35, log(dlmax), labels = dltxt, pos=1, col=dl.col, cex=0.8)
    }
    else {              # log scale, with groups
      nonas <- na.omit(data.frame(y1, y2, group))  # omits NAs for the sake of the ros function
      y.nona <- nonas[,1]
      nd <- (nonas[,2] == 1)   # nd is vector of logical T/F
      y.grp <- nonas[,3]
      grp <- as.factor(y.grp)
      gname <- deparse(substitute(group))
      if (is.null(ordr) == FALSE)  {
        grp = factor(grp, levels(grp)[ordr])
      }
      glevels <- levels(grp)
      gpnum <- length(levels(grp))

      xx <- c(0.5, gpnum+0.5, gpnum+0.5, 0.5)
      xgrp <- 0
      ygrp <- 0
      maxDL = 0
      y.rosmin <- c(1000000)
      j=0
      for (i in 1:gpnum) {
        if(sum(1-as.integer(nd[y.grp==levels(grp)[i]])) <=2 )  {
          cat("Note: One group had fewer than 3 detects and so cannot be plotted.", "\n")
          #       maxDL[i] <- 0
          maxDL[i] <- max(y.nona[y.grp == levels(grp)[i]] * as.integer(nd[y.grp == levels(grp)[i]]) )
          k=length(y.nona[y.grp == levels(grp)[i]])
          y.spacer <- rep(NA, k)
          if (j==0) {y.all <- y.spacer
          grp.all <- rep(levels(grp)[i],k)   }
          else {y.all <- c(y.all, y.spacer)
          grp.all <- c(grp.all, rep(levels(grp)[i],k))   }
          j <- k
          next }

        maxDL[i] <- max(y.nona[y.grp == levels(grp)[i]] * as.integer(nd[y.grp == levels(grp)[i]]) )
        y.ros <- suppressWarnings(ros(y.nona[y.grp == levels(grp)[i]], nd[y.grp == levels(grp)[i]]))
        y.rosmin[i] <- min(y.ros$modeled)
        k=length(y.ros$modeled)
        if (j==0) {y.all <- y.ros$modeled
        grp.all <- rep(levels(grp)[i],k)   }
        else {y.all <- c(y.all, y.ros$modeled)
        grp.all <- c(grp.all, rep(levels(grp)[i],k))   }
        j <- k
      }
      ymin.all <- min(y.rosmin[is.na(y.rosmin) == FALSE])

      grp.all <- as.factor (grp.all)
      if (is.null(ordr) == FALSE)  {
        grp.all = factor(grp.all, levels(grp.all)[ordr])
      }
      glevels <- levels(grp.all)
      gpnum <- length(levels(grp.all))

      yy = c(log(ymin.all), log(ymin.all), log(dlmax), log(dlmax))
      Ylab = c(Ylab, "(natural logs)" )
      if (is.null(Ymax) == FALSE) { Ylim = c(log(ymin.all), log(Ymax))}
      if (is.null(minmax)) {
        # boxplot with separate outliers beyond 1.5*IQR
         boxplot(log(y.all)~grp.all, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplecol = "white", outcex = 0.8, main =Title, col=bxcol, ylim = Ylim)}
      else {
        #  boxplot drawn to max and min
         boxplot(log(y.all)~grp.all, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplewex = 0.2, range=0,  main =Title, col=bxcol, ylim = Ylim)}

      multiDL <- ifelse (sd(maxDL) != 0, TRUE, FALSE)
      if (multiDL == FALSE) {     # plot all with same max DL
        # Everything below the max DL is gray, not black
         polygon(xx, yy, col = bdl.col, border=bdl.col)
         abline(h=log(dlmax), col = dl.col, lty="longdash",lwd=2)
         legend(dl.loc, legend = dltxt, bty="n", text.col=dl.col, cex=0.8)
      }
      else {  cat("Maximum DL","Group", "\n", sep = "   ")
        for (i in 1:gpnum) {xgrp <- c(i-0.5, i+0.5, i+0.5, i-0.5)   # different max DL per group
        ygrp = c(log(ymin.all), log(ymin.all), log(maxDL[i]), log(maxDL[i]))
        polygon(xgrp, ygrp, col = bdl.col, border = bdl.col)
        lx = c(i-0.45, i+0.45)
        ly = c(log(maxDL[i]), log(maxDL[i]))
        lines(lx, ly, col = dl.col, lty="longdash", lwd=1.5)
        cat("  ", maxDL[i], levels(grp.all)[i], "\n", sep = "   ")
        }
      }  }
  }

  else  {     #  LOG = FALSE
    if (is.null(group) == TRUE) {     # no groups
      xx = c(0.5, 1.5, 1.5, 0.5)
      y.ros <- suppressWarnings(ros(y.nona, nd))
      y.min <- min(y.ros$modeled)
      yy = c(y.min, y.min, dlmax, dlmax)
      if (is.null(Ymax) == FALSE) { Ylim = c(y.min, Ymax)}
      if (is.null(minmax)) {
        # boxplot with separate outliers beyond 1.5*IQR
        boxplot(y.ros$modeled, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplecol = "white", outcex = 0.8, main =Title, col=bxcol, ylim = Ylim)}
      else {
        #  boxplot drawn to max and min
        boxplot(y.ros$modeled, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplewex = 0.2, range=0,  main =Title, col=bxcol, ylim = Ylim)}

      # Everything below the max DL is gray, not black
      polygon(xx, yy, col = bdl.col, border=bdl.col)
      abline (h=y.min, col=bdl.col, lwd=8)
      abline(h=dlmax, col = dl.col, lty="longdash", lwd=2)
      text(1.35, dlmax, labels = dltxt, pos=3, col=dl.col, cex=0.8)
    }

    else {            #  LOG = FALSE, with groups
      nonas <- na.omit(data.frame(y1, y2, group))  # omits NAs for the sake of the ros function
      y.nona <- nonas[,1]
      nd <- (nonas[,2] == 1)   # nd is vector of logical T/F
      y.grp <- nonas[,3]

      grp <- as.factor(y.grp)
      gname <- deparse(substitute(group))
      gpnum <- length(levels(grp))
      glevels <- levels(grp)
      if (is.null(ordr) == FALSE)  {
        grp = factor(grp, levels(grp)[ordr])
      }
      levels(grp.all) <- glevels
      xx <- c(0.5, gpnum+0.5, gpnum+0.5, 0.5)
      xgrp <- 0
      ygrp <- 0
      maxDL = 0
      y.rosmin <- c(1000000)
      j=0
      for (i in 1:gpnum) {
        if(sum(1-as.integer(nd[y.grp==levels(grp)[i]])) <=2 )  {
          cat("Note: One group had fewer than 3 detects and so cannot be plotted.", "\n")
          #       maxDL[i] <- 0
          maxDL[i] <- max(y.nona[y.grp == levels(grp)[i]] * as.integer(nd[y.grp == levels(grp)[i]]) )
          k=length(y.nona[y.grp == levels(grp)[i]])
          y.spacer <- rep(NA, k)
          if (j==0) {y.all <- y.spacer
          grp.all <- rep(levels(grp)[i],k)   }
          else {y.all <- c(y.all, y.spacer)
          grp.all <- c(grp.all, rep(levels(grp)[i],k))   }
          j <- k
          next }

        maxDL[i] <- max(y.nona[y.grp == levels(grp)[i]] * as.integer(nd[y.grp == levels(grp)[i]]) )
        y.ros <- suppressWarnings(ros(y.nona[y.grp == levels(grp)[i]], nd[y.grp == levels(grp)[i]]))
        y.rosmin[i] <- min(y.ros$modeled)
        k=length(y.ros$modeled)
        if (j==0) {y.all <- y.ros$modeled
        grp.all <- rep(levels(grp)[i],k)   }
        else {y.all <- c(y.all, y.ros$modeled)
        grp.all <- c(grp.all, rep(levels(grp)[i],k))   }
        j <- k
      }
      ymin.all <- min(y.rosmin[is.na(y.rosmin) == FALSE])
      yy = c( ymin.all, ymin.all, dlmax, dlmax)
      grp.all <- as.factor (grp.all)
      if (is.null(ordr) == FALSE)  {
        grp.all = factor(grp.all, levels(grp.all)[ordr]) }
      else {levels(grp.all) <- glevels}

      if (is.null(Ymax) == FALSE) { Ylim = c(ymin.all, Ymax)}
      if (is.null(minmax)) {
        # boxplot with separate outliers beyond 1.5*IQR
        boxplot(y.all~grp.all, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplecol = "white", outcex = 0.8, main =Title, col=bxcol, ylim = Ylim)}
      else {
        #  boxplot drawn to max and min
        boxplot(y.all~grp.all, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplewex = 0.2, range=0,  main =Title, col=bxcol, ylim = Ylim)}

      multiDL <- ifelse (sd(maxDL) != 0, TRUE, FALSE)
      if (multiDL == FALSE) {      # all groups use same max DL
        polygon(xx, yy, col = bdl.col, border=bdl.col)
        abline (h=0, col=bdl.col, lwd=8)
        abline(h=dlmax, col = dl.col, lty="longdash", lwd=2)
        legend(dl.loc, legend = dltxt, bty="n", text.col=dl.col, cex=0.8)
      }
      else {cat("Maximum DL","Group", "\n", sep = "   ")
        for (i in 1:gpnum) {xgrp <- c(i-0.5, i+0.5, i+0.5, i-0.5)   # different max DL per group
        ygrp = c(ymin.all, ymin.all, maxDL[i], maxDL[i])
        polygon(xgrp, ygrp, col = bdl.col, border = bdl.col)
        lx = c(i-0.45, i+0.45)
        ly = c(maxDL[i], maxDL[i])
        lines(lx, ly, col = dl.col, lty="longdash", lwd=1.5)
        cat("  ", maxDL[i], levels(grp.all)[i], "\n", sep = "   ")
        }
      }
    }
  }
 }   # end of when there are nondetects

  else   # when there are no nondetects
  { LOG <- ifelse (LOG, "y", "")
    if (is.null(group) == TRUE)    # no group
    {  boxplot(y1, na.action = na.omit, ylab = Ylab, col = bxcol, main = Title, log=LOG)}
    else        # with groups
    {gname <- deparse(substitute(group))
        if (is.null(ordr) == FALSE)  {
          group = factor(group, levels(group)[ordr]) }
        glabs <- levels(group)
        boxplot(y1~group, na.action = na.omit, ylab = Ylab, xlab = gname, names = glabs, col = bxcol, main = Title, log=LOG)
        }
  }
}

#' Censored data paired t-test
#'
#' @description
#'Performs a parametric test of whether the mean difference between two columns of paired censored data equals 0. Assumes that the paired differences follow a gaussian (normal) distribution.
#' @param xd The first column of data values plus detection limits
#' @param xc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the xd column, and 0 (or `FALSE`) indicates a detected value in xd.
#' @param yd The second column of data values plus detection limits, or a single number representing a standard / guideline value.
#' @param yc The column of censoring indicators for yd, where 1 (or `TRUE`) indicates a detection limit in the yd column, and 0 (or `FALSE`) indicates a detected value in `yd`. Not needed if `yd` is a single standard number.
#' @param alternative The usual notation for the alternate hypothesis.  Default is `“two.sided”`.  Options are `“greater”` or `“less”`.
#'
#' @importFrom fitdistrplus plotdistcens
#' @importFrom survival survreg Surv
#' @export
#' @return  A list of statistics containing the following components:
#' \itemize{
#' \item `n` Number of samples
#' \item `Z` The value of the test statistic
#' \item `p.value` the p-value of the test
#' \item `Mean difference` the mean difference between `xd` and `yd`
#' }
#' @details You may also test for whether the x data exceed a standard by entering the single number for the standard as `yd`.  In that case no `yc` is required.
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @seealso [survival::survreg]
#'
#' @examples
#'
#' data(atrazine)
#'
#' cen_paired(atrazine$June,atrazine$JuneCen,atrazine$Sept,atrazine$SeptCen)
#'
#' # Comparing standard/guieline value
#' cen_paired(atrazine$June, atrazine$JuneCen, 0.01, alternative = "greater")

cen_paired <- function(xd, xc, yd, yc, alternative="two.sided") {
  xname <- deparse(substitute(xd))
  yname <- deparse(substitute(yd))
  mu <- 0
  if (length(yd)==1) {yc <- rep(FALSE, length(xd))
  mu <- yd[1]; yd <- rep(mu, length(xd))}
  nonas <- na.omit(data.frame(xd, xc, yd, yc))
  if(length(nonas[,1]) != length(nonas[,3]))
    stop("Lengths of x and y must be the same for paired data.")
  if(alternative == "two.sided")
  {txt3 <- paste("alternative hypothesis: true mean difference does not equal ", mu, ".", sep="")}
  else if(alternative == "less")
  {txt3 <-  paste("alternative hypothesis: true mean difference is less than ", mu, ".", sep="")}
  else if (alternative == "greater")
  {txt3 <-  paste("alternative hypothesis: true mean difference is greater than ", mu, ".", sep="")}
  else {stop('Invalid choice for alternative, must match "two.sided", "greater", or "less"')}
  N = length(nonas[,1])
  direction <- 0

  ## Create min and max differences
  d1 <- nonas[,1] - ifelse(nonas[,4], 0, nonas[,3])
  d2 <- ifelse(nonas[,2], 0, nonas[,1]) - nonas[,3]
  mind <- pmin(d1, d2)
  maxd <- pmax(d1, d2)
  #  print(data.frame(xd, as.logical(xc), yd, as.logical(yc), mind, maxd))

  ## Mean difference
  sv.out <- survreg(Surv(mind, maxd, type="interval2")~1, dist = "gaussian")
  mean.diff <- sv.out$coefficients

  #direction of mean vs expected
  if (alternative == "less") {direction <- -1}
  if (alternative == "greater") {direction <- 1}
  direction <- direction*sign(mean.diff)
  if (direction <0) {direction <- 0}

  test.diff <- summary(sv.out)
  Z <- test.diff$table[1,3]
  pvalue <- test.diff$table[1,4]
  if(alternative == "two.sided") {
    if (mu == 0) {txt3 <- paste("alternative hypothesis: true mean difference does not equal ", mu, ".", sep="")
    txt <- paste("Censored paired test for mean(", xname, " - ", yname, ") equals ", mu, ".", sep = "")}
    else {txt3 <-  paste("alternative hypothesis: true mean does not equal ", mu, ".", sep ="")
    txt <- paste("Censored paired test for mean(", xname, ") equals ", mu, sep = "")}
  }
  else if(alternative == "less")
  { pvalue <- ifelse (direction, pvalue/2, 1-(pvalue/2))
  if (mu == 0) {txt3 <-  paste("alternative hypothesis: true mean difference is less than ", mu, ".", sep="")
  txt <- paste("Censored paired test for mean(", xname, " - ", yname, ") equals ", mu, ".", sep = "")}
  else {txt3 <-  paste("alternative hypothesis: true mean ", xname, " is less than ", mu, ".", sep="")
  txt <- paste("Censored paired test for mean(", xname, ") equals ", mu, sep = "")}
  }
  else   # (alternative == "greater")
  { pvalue <- ifelse(direction, pvalue/2, 1-(pvalue/2) )
  if (mu == 0) {txt3 <-  paste("alternative hypothesis: true mean difference is greater than ", mu, ".", sep="")
  txt <- paste(" Censored paired test for mean(", xname, " - ", yname, ") equals ", mu, ".", sep = "")}
  else {txt3 <-  paste("alternative hypothesis: true mean ", xname, " exceeds ", mu, ".", sep = "")
  txt <- paste(" Censored paired test for mean(", xname, ") equals ", mu, sep = "")}
  }
  txt2 <- paste("n =", N, "  Z=", round(Z, 4), "  p-value =", signif(pvalue, 4))
  if (mu == 0) {cat( txt, "\n", txt3, "\n", "\n", txt2, "\n", "Mean difference =", signif(mean.diff,4), "\n")}
  else {cat( txt, "\n", txt3, "\n", "\n", txt2, "\n",paste("Mean", xname), "=", signif(mean.diff+mu,4), "\n")}

  #plot the cdf of differences
  sv.coefs <- data.frame(sv.out$icoef[1], exp(sv.out$icoef[2]))
  names(sv.coefs) <- c("mean", "sd")
  left <- mind;  right <- maxd
  minmax <- data.frame(left, right)

  oldpar<- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  plotdistcens(minmax, distr = "norm", para = sv.coefs, main = "Differences: CDF and Fitted Normal Dist")

}

#' Wilcoxcon Signed-Rank test for censored data
#'
#' @description Performs a nonparametric Wilcoxon signed-rank test of whether the median difference between two columns of paired censored data equals 0.  Uses the Pratt adjustment for pairs of equal values.
#' @param xd The first column of data values plus detection limits
#' @param xc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the xd column, and 0 (or `FALSE`) indicates a detected value in xd.
#' @param yd The second column of data values plus detection limits, or a single number representing a standard / guideline value.
#' @param yc The column of censoring indicators for yd, where 1 (or `TRUE`) indicates a detection limit in the yd column, and 0 (or `FALSE`) indicates a detected value in `yd`. Not needed if `yd` is a single standard number.
#' @param alternative The usual notation for the alternate hypothesis.  Default is `“two.sided”`.  Options are `“greater”` or `“less”`.
#'
#' @importFrom coin wilcoxsign_test statistic pvalue
#' @export
#' @return  Prints a list of Wilcoxon Signed-Rank test with Pratt correction for ties statistics containing the following components:
#' \itemize{
#' \item `n` Number of samples
#' \item `Z` The value of the test statistic
#' \item `p.value` the p-value of the test
#' }
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Page, E.B., 1963. Ordered Hypotheses for Multiple Treatments: A Significance Test for Linear Ranks. Journal of the American Statistical Association 58, 216–230. <https://doi.org/10.2307/2282965>
#'
#' Pratt, J.W., 1959. Remarks on Zeros and Ties in the Wilcoxon Signed Rank Procedures. Journal of the American Statistical Association 54, 655–667. <https://doi.org/10.2307/2282543>
#'
#'
#' @examples
#'
#' data(atrazine)
#'
#' cen_signedranktest(atrazine$June,atrazine$JuneCen,atrazine$Sept,atrazine$SeptCen)

cen_signedranktest <- function(xd, xc, yd, yc, alternative="two.sided") {
  xname <- deparse(substitute(xd))
  yname <- deparse(substitute(yd))
  nonas <- na.omit(data.frame(xd, xc, yd, yc))
  if(length(nonas[,1]) != length(nonas[,3]))
    stop("Lengths of x and y must be the same for paired data.")
  if(alternative == "two.sided")
  {txt3 <- paste("alternative hypothesis: true difference", xname, "-", yname, "does not equal 0")}
  else if(alternative == "less")
  {txt3 <-  paste("alternative hypothesis: true difference", xname, "-", yname, "is less than 0")}
  else if (alternative == "greater")
  {txt3 <-  paste("alternative hypothesis: true difference", xname, "-", yname, "is greater than 0")}
  else {stop('Invalid choice for alternative, must match "two.sided", "greater", or "less"')}
  N = length(nonas[,1])
  for(i in seq(N)) {
    if(nonas[i,2] && nonas[i,4]) { # Both censored, make same
      xyMax <- max(nonas[i,1], nonas[i,3])
      nonas[i,1] <- nonas[i,3] <- xyMax
    }
    else if(nonas[i,2] && nonas[i,1] > nonas[i,3]) { # x censored > y
      nonas[i,3] <- nonas[i,1]
      nonas[i,4] <- TRUE
    }
    else if(nonas[i,4] && nonas[i,3] > nonas[i,1]) { # y censored > x
      nonas[i,1] <- nonas[i,3]
      nonas[i,2] <- TRUE
    }
  } # Done, no need to check uncensored observations
  # setting <s to slightly below detected obs at that value
  nonas[,1] <- nonas[,1] - nonas[,2] * (0.001) * nonas[,1]
  nonas[,3] <- nonas[,3] - nonas[,4] * (0.001) * nonas[,3]

  x <- nonas[,1]
  names(x) <- xname
  y <- nonas[,3]
  names(y) <- yname

  s.out <- wilcoxsign_test(x~y, alternative = alternative)
  txt <- paste("Censored signed-rank test for (x:", xname, " - ", "y:", yname, ") equals 0", "\n", txt3, "\n", sep = "")
  txt2 <- paste("n =", N, "  Z=", signif(statistic(s.out), 4), "  p-value =", signif(pvalue(s.out), 4))

  cat(txt, "\n","Pratt correction for ties", "\n", txt2, "\n")
}

#' Sign test for censored data
#'
#' @description Performs a nonparametric sign test of whether the median difference between two columns of paired censored data equals 0. Uses the Fong adjustment for pairs of equal values.
#' @param xd The first column of data values plus detection limits
#' @param xc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the xd column, and 0 (or `FALSE`) indicates a detected value in xd.
#' @param yd The second column of data values plus detection limits.
#' @param yc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the `yd` column, and 0 (or `FALSE`) indicates a detected value in `yd`.
#' @param alternative The usual notation for the alternate hypothesis.  Default is `“two.sided”`.  Options are `“greater”` or `“less”`.
#'
#' @importFrom stats dbinom pbinom binom.test
#' @export
#' @return  Returns the number of `xd` and `yd` values greater than, less than and tied. Corrected and uncorrected `p-value` for ties also displayed.
#' @references
#' Fong, D.Y.T., Kwan, C.W., Lam, K.F., Lam, K.S.L., 2003. Use of the Sign Test for the Median in the Presence of Ties. The American Statistician 57, 237–240.
#'
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'#'
#' @examples
#'
#' data(atrazine)
#'
#' cen_signtest(atrazine$June,atrazine$JuneCen,atrazine$Sept,atrazine$SeptCen)


cen_signtest <- function(xd, xc, yd, yc, alternative="two.sided") {
  xname <- deparse(substitute(xd))
  yname <- deparse(substitute(yd))
  nonas <- na.omit(data.frame(xd, xc, yd, yc))
  if(length(nonas[,1]) != length(nonas[,3]))
    stop("Lengths of x and y must be the same for paired data.")
  if(alternative == "two.sided")
  {txt3 <- paste("  alternative hypothesis: true median difference is not = 0"); dir <- 0}
  else if(alternative == "less")
  {txt3 <-  paste("  alternative hypothesis: true median difference < 0"); dir <- (-1)}
  else if (alternative == "greater")
  {txt3 <-  paste("  alternative hypothesis: true median difference > 0"); dir <- 1}
  else {stop('Invalid choice for alternative, must match "two.sided", "greater", or "less"')}
  N = length(nonas[,1])
  for(i in seq(N)) {
    if(nonas[i,2] && nonas[i,4]) { # Both censored, make same
      xyMax <- max(nonas[i,1], nonas[i,3])
      nonas[i,1] <- nonas[i,3] <- xyMax
    }
    else if(nonas[i,2] && nonas[i,1] > nonas[i,3]) { # x censored > y
      nonas[i,3] <- nonas[i,1]
      nonas[i,4] <- TRUE
    }
    else if(nonas[i,4] && nonas[i,3] > nonas[i,1]) { # y censored > x
      nonas[i,1] <- nonas[i,3]
      nonas[i,2] <- TRUE
    }
  } # Done, no need to check uncensored observations
  # setting <s to slightly below detected obs at that value
  nonas[,1] <- nonas[,1] - nonas[,2] * (0.001) * nonas[,1]
  nonas[,3] <- nonas[,3] - nonas[,4] * (0.001) * nonas[,3]

  x <- nonas[,1]
  names(x) <- xname
  y <- nonas[,3]
  names(y) <- yname

  # test.dat <- data.frame(x, nonas[2], y, nonas[4])
  # print(test.dat)
  x.gt.y <- sum(x > y)
  x.lt.y <- sum (y > x)
  ties <- length(nonas[,1]) - x.gt.y - x.lt.y
  #  result.counts <- data.frame(x.gt.y, x.lt.y, ties)
  #  print(result.counts)
  s.binom <- binom.test(x.gt.y, length(nonas[,1])-ties, p=0.5, alternative = alternative)

  # LSA version in the coin package.  Not currently used.
  #  s.out <- sign_test(x~y, alternative = alternative)
  #  txt.s <- paste("n =", N, "  Z=", round(statistic(s.out), 4), "  p-value =", round(pvalue(s.out), 4))
  #  print(txt.s)

  # Fong's Modified Sign Test correction for ties
  dir <- dir*sign(x.gt.y - x.lt.y)  #dir=1: alt in same direction as data

  num.Fong <- 1-pbinom(max(x.gt.y, x.lt.y)-1,  N, prob=0.5)
  denom.Fong <- 1-pbinom(floor((N-ties+1)/2)-1, N, prob=0.5)
  p.Fong <- num.Fong/denom.Fong
  if (dir== 1) {num.Fong <- 1-pbinom((max(x.gt.y, x.lt.y)-1), N, prob=0.5)
  denom.Fong <- 2*(1-pbinom(floor((N-ties+1)/2)-1, N, prob=0.5))
  p.Fong <- num.Fong/denom.Fong}
  if (dir== -1) p.Fong <- 1-p.Fong/2+dbinom(min(x.gt.y, x.lt.y), N-ties, prob=0.5)
  #{num.Fong <- pbinom((max(x.gt.y, x.lt.y)), N, prob=0.5)
  # denom.Fong <- 2*(1-pbinom(floor((N-ties+1)/2)-1, N, prob=0.5))
  # p.Fong <- num.Fong/denom.Fong}


  txt <- paste("Censored sign test for median(x:", xname, " - ", "y:", yname, ") equals 0", sep = "")
  txt2 <- paste("  n =", N, "  n+ =", x.gt.y, "  n- =", x.lt.y, "   ties:", ties, "\n")
  cat(txt, "\n",txt3, "\n", txt2, "\n", " No correction for ties:", "  p-value =", signif(s.binom$p.value, 4), "\n")
  if (ties != 0) cat("Fong correction for ties:", "  p-value =", signif(p.Fong, 4))
}

#' Peto-Peto one-factor test
#'
#' @description
#' Performs a Peto-Peto nonparametric test of differences in cdfs between groups.  If more than two groups, the test is followed by a nonparametric multiple comparison test.  Uses the BH method of adjusting p-values.
#' @param y1 The column of data values plus detection limits
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the y1 column, and 0 (or `FALSE`) indicates a detected value in y1.
#' @param grp Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @importFrom survminer pairwise_survdiff
#' @importFrom stats pchisq
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
#' Peto, R., Peto, J., 1972. Asymptotically Efficient Rank Invariant Test Procedures. Journal of the Royal Statistical Society. Series A (General) 135, 185. <https://doi.org/10.2307/2344317>
#'
#'Benjamini, Y., Hochberg, Y., 1995. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.  Journal of the Royal Statistical Society. Series B (Methodological), 57, 289-300.
#'
#' @importFrom survival survdiff Surv
#'
#' @examples
#' data(PbHeron)
#'
#' # Two Groups
#' cen1way(PbHeron$Liver,PbHeron$LiverCen,PbHeron$DosageGroup)
#'
#' # More than two groups
#' cen1way(PbHeron$Liver,PbHeron$LiverCen,PbHeron$Group)


cen1way <- function(y1,y2, grp) {
  yname <- deparse(substitute(y1))
  gname <- deparse(substitute(grp))
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
    Cstats <- suppressWarnings(cfit(y1gp, as.logical(y2gp), printstats=FALSE, Cdf = FALSE))
    Cstats <- Cstats[c(1:6)]
    Cstats <- Cstats[-3]
    Cstats <- data.frame(Cstats)
    if (i ==1) {Cen.stats <- Cstats
    cnames <- colnames(Cstats)}
    else {Cen.stats <- rbind(Cen.stats, Cstats)}
  }
  rownames(Cen.stats) <-  groupnames
  colnames(Cen.stats) <- cnames
  print.data.frame(Cen.stats, print.gap = 3)
  cat('\n',"     Oneway Peto-Peto test of CensData:", yname, "  by Factor:", gname, '\n', "     Chisq =", signif(y.out$chisq, 4), "  on", df, "degrees of freedom", "    p =", signif(pval,3), '\n')
  if (df >1) {mcomp <- pairwise_survdiff(Surv(flip, detect) ~ Factor, data=CensData, rho=rho)
  print(mcomp)}
}

#' Censored data two-group test for difference in means
#'
#' @description
#' Performs a parametric test of differences in means between two groups of censored data, either in original or in log units (the latter becomes a test for difference in geometric means).
#' @param y1 The column of data values plus detection limits
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the y1 column, and 0 (or `FALSE`) indicates a detected value in y1.
#' @param grp Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param LOG Indicator of whether to compute tests in the original units, or on their logarithms.  The default is to use the logarithms (LOG = `TRUE`).  To compute in original units, specify the option LOG = `FALSE` (or LOG = 0).
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


cen2means <- function(y1, y2, grp, LOG=TRUE) {
  yname <- deparse(substitute(y1))
  gname <- deparse(substitute(grp))
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
    cat("     MLE 't-test' of mean natural logs of CensData:", yname, "by Factor:", gname, '\n', "   ",dist.test,'\n')
    cat("     geometric mean of", grpname[1], "=", signif(mean1, 4), "    geometric mean of", grpname[2], "=", signif(mean2,4), "\n")
    cat( "     Chisq =", signif(reg.chisq, 4), " on", df, "degrees of freedom", "    p =", signif(pval,3), '\n', "\n")
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
  cat("     MLE 't-test' of mean CensData:", yname, "  by Factor:", gname, '\n', "   ",dist.test,'\n')
  cat("     mean of", grpname[1], "=", signif(mean1, 4), "    mean of", grpname[2], "=", signif(mean2,4), "\n")
  cat("     Chisq =", signif(reg.chisq, 4), " on", df, "degrees of freedom", "    p =", signif(pval,3), '\n')
  # A warning
  cat("\n", "  NOTE: Data with nondetects may be projected below 0 with MLE normal distribution.", "\n", "  If so, p-values will be unreliable (often too small).  Use perm test instead.", "\n")
  # Q-Q plot of residuals
  reg.predict <- predict(reg.out)
  two.group <- reg.predict - flip
  cenregQQ(two.group, as.logical(y2), Factor, LOG = FALSE)
  }
  return(invisible(result))
}

#' ANOVA for censored data
#'
#' @description Performs a parametric test of differences in means between groups of censored data, followed by a parametric Tukey's multiple comparison test.
#' @param y1 The column of data values plus detection limits
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y1` column, and 0 (or `FALSE`) indicates a detected value in `y1`.
#' @param grp Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param LOG Indicator of whether to compute tests in the original units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @importFrom survival survreg Surv
#' @importFrom multcomp glht mcp
#' @importFrom stats residuals
#' @export
#' @return
#' Returns the Maximum Likelihood Estimation (MLE) comparison results including Chi-Squared value, degrees of freedom and `p-value` of the test. Test assumes lognormal(`LOG=TRUE`) or nomal(`LOG=FALSE`) distribution of residuals from group means.
#'
#' Tukey's multiple comparison p-values of pairwise differenes in group means are also printed.
#' \itemize{
#' \item Group Names of groups (NOTE: `== 0` indicates null hypothesis of "equals zero").
#' \item `Estimate` Estimated difference between group means.
#' \item `Std. Error` Standard error of estimate.
#' \item `z value` Test statistic.
#' \time `Pr(>|z|)` P-values for test that difference in means equals zero.
#' }
#'
#' @details Test is computed using Maximum Likelihood Estimation. When a gaussian	distribution model is	used (LOG=FALSE) modeled values may fall below zero, producing unreal p-values (often lower than they should be).  Because of this, testing in log units is preferable and is the default.
#' @seealso [survival::survreg]
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#'#' @examples
#' data(PbHeron)
#' cenanova(PbHeron$Liver,PbHeron$LiverCen,PbHeron$DosageGroup)
#'
#' cenanova(PbHeron$Liver,PbHeron$LiverCen,PbHeron$DosageGroup,LOG=FALSE)
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

#' Comparison of empirical cdf of censored data
#'
#' @description Plots the empirical cdf and cdfs of three theoretical distributions, fit by maximum likelihood estimation (MLE).
#' @param y.var The column of y (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param dist3 Name of the third distribution to be plotted, default is `norm` (normal distrubtion). Alternate third distribution is `weibull`(for Weibull).  Lognormal and gamma distributions are always used.
#' @param Yname Optional – input text in quotes to be used as the variable name.  The default is the name of the `y.var` input variable.
#' @export
#' @return prints a plot of the empirial CDFs with BIC value for each distribution.
#' @importFrom fitdistrplus cdfcompcens fitdistcens
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Delignette-Muller, M., Dutang, C., 2015. fitdistrplus : An R Package for Fitting Distributions. Journal of Statistical Software, 64, 1-34. http://www.jstatsoft.org/v64/i04/.
#'
#' @examples
#'
#' library(NADA) #For example data
#'
#' data(Brumbaugh)
#' cenCompareCdfs(Brumbaugh$Hg,Brumbaugh$HgCen)
#'
#' # With Weibull distribution
#' cenCompareCdfs(Brumbaugh$Hg,Brumbaugh$HgCen,dist3="weibull")
#'
#' # Using an distribution not supported by this function (yet)
#' # you will get an error message
#' \dontrun{cenCompareCdfs(Brumbaugh$Hg,Brumbaugh$HgCen,dist3="beta")}
#'
#' # With Yname specified
#' cenCompareCdfs(Brumbaugh$Hg,Brumbaugh$HgCen,Yname="TCE Conc (ug/L)\nLong Island, NY USA")


cenCompareCdfs <- function(y.var, cen.var, dist3="norm", Yname = yname)  {
  #added to stop if dist3 is not from the list
  if(!(dist3%in%c("norm","lnorm","gamma","weibull"))){stop(paste0(dist3," distribution is not supported with this function, try again."))}

  dist.vals <- c("norm","lnorm","gamma","weibull")
  dist.vals.text <- c("Normal","Lognormal","Gamma","Weibull")

  yname <- deparse(substitute(y.var))
  if (sum(as.integer(cen.var)) > 0)    # not all data are detects

  {left <- y.var*(1-as.integer(cen.var))
  right <- y.var
  var.frame <- data.frame(left, right)

  y.dist1 <- fitdistcens(var.frame, "lnorm")
  y.dist2 <- fitdistcens(var.frame, "gamma")
  y.dist3 <- fitdistcens(var.frame, dist3)

  bic.dist1 <- paste("Lognormal BIC =", signif(y.dist1$bic, 3) )
  bic.dist2 <- paste("Gamma BIC =", signif(y.dist2$bic, 3) )
  bic.dist3 <- paste(dist.vals.text[match(dist3,dist.vals)],"BIC =", signif(y.dist3$bic, 3) )

  cdfcompcens(list(y.dist1, y.dist2, y.dist3), legendtext=c(bic.dist1, bic.dist2, bic.dist3), xlab = Yname, fitlty = c(1, 5, 3), lwd = c(1, 1, 2))

}
    else            # all data are detects
  {
  y.dist1 <- fitdist(var.frame, "lnorm", "mle")
  y.dist2  <- fitdist(var.frame, "gamma", "mle")
  y.dist3  <- fitdist(var.frame, dist3, "mle")

  bic.dist1 <- paste("Lognormal BIC =", signif(y.dist1$bic, 3) )
  bic.dist2 <- paste("Gamma BIC =", signif(y.dist2$bic,3) )
  bic.dist3 <- paste(dist.vals.text[match(dist3,dist.vals)],"BIC =", signif(y.dist3$bic, 3) )

    cdfcomp(list(y.dist1, y.dist2, y.dist3), legendtext=c(bic.dist1, bic.dist2, bic.dist3), do.points = FALSE, verticals = TRUE, xlab = Yname, fitlty = c(1, 5, 3), lwd = c(1, 1, 2))
  }

  # prior version edited by PJ
  # y.lnorm <- fitdistcens(var.frame, "lnorm")
  # y.gamma <- fitdistcens(var.frame, "gamma")
  # y.norm <- fitdistcens(var.frame, "norm")
  # bic.lnorm <- paste("Lognormal BIC =", signif(y.lnorm$bic, 3) )
  # bic.gamma <- paste("Gamma BIC =", signif(y.gamma$bic,3) )
  # bic.norm <- paste("Normal BIC =",signif(y.norm$bic, 3) )
  #
  # if (dist3 == "normal") {
  #   cdfcompcens(list(y.lnorm, y.gamma, y.norm), legendtext=c(bic.lnorm, bic.gamma, bic.norm), xlab = Yname, fitlty = c(1, 5, 3), lwd = c(1, 1, 2))
  # }   else {
  #   y.weib <- fitdistcens(var.frame, "weibull")
  #   bic.weib <- paste("Weibull BIC =",signif(y.weib$bic, 3) )
  #   cdfcompcens(list(y.lnorm, y.gamma, y.weib), legendtext=c(bic.lnorm, bic.gamma, bic.weib), xlab = Yname, fitlty = c(1, 5, 3), lwd = c(1, 1, 2))
  # }
  }

#' Censored Q-Q Plot comparison
#'
#' @description Produces three quantile-quantile (Q-Q) plots, also called probability plots, based on three distributions (normal, lognormal and gamma distributions).
#' @param y.var The column of y (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param Yname Optional – input text in quotes to be used as the variable name on all plots.  The default is the name of the `y.var` input variable.
#' @export
#' @return Plots three Q-Q plots based on normal, lognormal and gamma distributions and prints the best-fit distribution.
#' @details Produces three Q-Q plots and reports which has the highest Shapiro-Francia test statistic (W).  The distribution with the highest W is the best fit of the three.
#'
#' @importFrom EnvStats distChooseCensored qqPlotCensored
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Millard, S.P., 2013. EnvStats: An R Package for Environmental Statistics. Springer-Verlag, New York.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#' data(Brumbaugh)
#'
#' \dontrun{cenCompareQQ(Brumbaugh$Hg,Brumbaugh$HgCen)}

cenCompareQQ <- function(y.var, cen.var, Yname = yname)  {
  yname <- deparse(substitute(y.var))
  if (sum(as.integer(cen.var)) > 0)    # not all data are detects

  { cen.logical <- as.logical(cen.var)
  var.choose <- distChooseCensored(y.var, cen.logical)
  norm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$norm$statistic, 3) )
  lnorm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$lnorm$statistic, 3) )
  gamma.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$gamma$statistic, 3) )
  all.W <- c(var.choose$test.results$norm$statistic, var.choose$test.results$lnorm$statistic, var.choose$test.results$gamma$statistic)
  all.W <- all.W/ max(all.W)
  best.text <- c("normal", "lognormal", "gamma")
  max.distrib <- best.text[all.W==1.0]

  if (var.choose$decision != "Nonparametric") {
    best.dist <- paste (var.choose$decision, "is a good fit")}
  else { best.dist <- paste ("Best of the three distributions is the", max.distrib)
  }
  cat(best.dist, "\n")
  par(mfrow=c(2,2))
  qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
  mtext(norm.text)
  #  legend("bottomright", legend = norm.text)

  qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = Yname, main = "Lognormal Q-Q Plot")
  mtext(lnorm.text)
  #  legend("bottomright", legend = lnorm.text)

  qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
  mtext(gamma.text)
  #   legend("bottomright", legend = gamma.text)

  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.47, y = 0.6, best.dist, pos = 1, cex = 1.2, col = "black", family="sans", font=1, adj=1)
  par(mfrow=c(1,1))
  }

else  # all data are detects
{ cen.logical <- as.logical(cen.var)
var.choose <- distChoose(y.var, method = "sf", alpha = 0.05, choices = c("norm", "gamma", "lnorm"))
norm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$norm$statistic, 3) )
lnorm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$lnorm$statistic, 3) )
gamma.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$gamma$statistic, 3) )
all.W <- c(var.choose$test.results$norm$statistic, var.choose$test.results$lnorm$statistic, var.choose$test.results$gamma$statistic)
all.W <- all.W/ max(all.W)
best.text <- c("normal", "lognormal", "gamma")
max.distrib <- best.text[all.W==1.0]

if (var.choose$decision != "Nonparametric") {
  best.dist <- paste (var.choose$decision, "is a good fit")}
else { best.dist <- paste ("Best of the three distributions is the", max.distrib)
}
cat(best.dist, "\n")
par(mfrow=c(2,2))
EnvStats::qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
mtext(norm.text)

EnvStats::qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = Yname, main = "Lognormal Q-Q Plot")
mtext(lnorm.text)

EnvStats::qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
mtext(gamma.text)

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.47, y = 0.6, best.dist, pos = 1, cex = 1.2, col = "black", family="sans", font=1, adj=1)
par(mfrow=c(1,1))
  }
}

#' Correlation and Regression with censored data
#'
#' @description Computes three parametric correlation coefficients for one X variable and the corresponding R2 for multiple X variables, and a regression equation for censored data.
#' @param y.var The column of y (response variable) values plus detection limits.
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value is in `y.var`.
#' @param x.vars One or more uncensored explanatory variable(s). See Details
#' @param LOG Indicator of whether to compute the regression in the original y units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @export
#' @return
#' When `x.vars` is one variable, likelihood, rescaled likelihood and McFaddens correlation coefficient (`R`) are printed.
#' When `x.vars` is more than one variable, likelihood, rescaled likelihood and McFaddens coefficent of determination (`R2`) are printed.
#'
#' Model coefficients (intercept and slopes), Chi-Squared statistic and p-value for test that all slope coefficients equal zero (overall test), and model AIC and BIC are provided.
#'
#' Q-Q plot of model residuals with corresponding Shapiro-Francia W and p-value are plotted for evaluation of model distributional assumptions.
#'
#' @importFrom survival survreg Surv
#' @importFrom EnvStats gofTestCensored qqPlotCensored
#'
#' @details
#'
#' `x.vars`: If 1 x variable only, enter its name.  If multiple x variables, enter the name of a data frame of columns of the x variables. No extra columns unused in the regression allowed. Create this by `x.frame <- data.frame (Temp, Flow, Time)` for 3 variables (temperature, flow and time).
#'
#' AIC and BIC are printed to help evaluate the ‘best’ regression model.
#'
#' The default is that the Y variable will be log transformed.
#'
#' @seealso [survival::survreg]
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' @examples
#'
#' data(Brumbaugh)
#'
#' # One variable
#' cencorreg(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh$SedMeHg)
#'
#' # More than one variable for demostration purposes
#'cencorreg(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh[,c("SedMeHg","PctWetland")])


cencorreg <- function(y.var, cen.var, x.vars, LOG = TRUE, verbose = TRUE) {
  yname <- deparse(substitute(y.var))
  nonas <- na.omit(cbind(y.var, cen.var, x.vars))
  xnona <- nonas[,-(1:2)]

  if (LOG == TRUE)  {lnvar <- log(nonas[,1])    # take logs of Y (default)
  flip.log <- max(lnvar) +1 - lnvar
  #  print(max(lnvar)+1)
  surv.log <- Surv(flip.log, as.logical(1-nonas[,2]) )

  if (is.data.frame(x.vars))  {             # multiple x variables
    reg.out <- survreg(surv.log ~ ., data = xnona, dist = "gaussian")
    cn <- names(reg.out$coefficients[-1])
    xvars.txt <- cn[1]
    for (i in 1:length(cn))  {j <-(i+1)
    if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
    }
    reg.out$call[3] <- xvars.txt
  }

  else { xname <- deparse(substitute(x.vars))          # 1 x variable
  x.df <- as.data.frame(xnona)
  names(x.df) <- xname
  reg.out <- survreg(surv.log~., data = x.df, dist = "gaussian")
  cn <- names(reg.out$coefficients[-1])
  reg.out$call[3] <- cn
  }

  ylog.pred <- max(lnvar) +1 - reg.out$linear.predictors
  reg.out$call[2] <- paste("log(", yname, ")", sep = "")
  reg.out$coefficients <- reg.out$coefficients * (-1)
  reg.out$coefficients[1] <- max(lnvar)+1 + reg.out$coefficients[1]   # coeffs in cenreg lognormal
  ylog.resi <- lnvar - ylog.pred
  reg.out$linear.predictors <- ylog.pred
  reg.out$resids <- ylog.resi
  vtext<- paste("Quantiles of", yname, "residuals (log units)")

  if (verbose == TRUE) {
  testnorm <- gofTestCensored(ylog.resi,nonas[,2])
  ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
  qqPlotCensored(ylog.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Lognormal Q-Q Plot of residuals")
  mtext(ptext)
  } }  # end of taking logs of Y

  else                                          #  Y in original units, normal Q-Q plot
  { if(min(nonas[,1] >= 0)) { y.low <- nonas[,1]*(1-nonas[,2])}   #  0 for low end of all NDs
  surv.norm <- Surv(y.low, nonas[,1], type="interval2")
  if (is.data.frame(x.vars))  {       # multiple x variables
    reg.out <- survreg(surv.norm ~ ., data = xnona, dist = "gaussian")
    cn <- names(reg.out$coefficients[-1])
    xvars.txt <- cn[1]
    for (i in 1:length(cn))  {j <-(i+1)
    if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
    }
    reg.out$call[3] <- xvars.txt
  }

  else { xname <- deparse(substitute(x.vars))          # 1 x variable
  x.df <- as.data.frame(xnona)
  names(x.df) <- xname
  reg.out <- survreg(surv.norm~., data = x.df, dist = "gaussian")
  cn <- names(reg.out$coefficients[-1])
  reg.out$call[3] <- cn
  }

  reg.out$call[2] <- yname
  ynorm.pred <- reg.out$linear.predictors
  ynorm.resi <- y.low - ynorm.pred
  reg.out$linear.predictors <- ynorm.pred
  reg.out$resids <- ynorm.resi
  vtext<- paste("Quantiles of", yname, "residuals")

  if (verbose == TRUE) {
  testnorm <- gofTestCensored(ynorm.resi,nonas[,2])
  ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
  qqPlotCensored(ynorm.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Normal Q-Q Plot of residuals")
  mtext(ptext)
  } }  # end of low censored = 0

 #  negative y values.  Use flip variable to set -Inf as low end for censored values
 else{ flip.norm <- max(nonas[,1]) +1 - nonas[,1]
 surv.norm <- Surv(flip.norm, as.logical(1-nonas[,2]) )

 if (is.data.frame(x.vars))  {       # multiple x variables
   reg.out <- survreg(surv.norm ~ ., data = xnona, dist = "gaussian")
   cn <- names(reg.out$coefficients[-1])
   xvars.txt <- cn[1]
   for (i in 1:length(cn))  {j <-(i+1)
   if (i != length(cn)) xvars.txt <- paste(xvars.txt, cn[j], sep = "+")
   }
   reg.out$call[3] <- xvars.txt
 }

 else { xname <- deparse(substitute(x.vars))          # 1 x variable
 x.df <- as.data.frame(xnona)
 names(x.df) <- xname
 reg.out <- survreg(surv.norm~., data = x.df, dist = "gaussian")
 cn <- names(reg.out$coefficients[-1])
 reg.out$call[3] <- cn
 }

 reg.out$call[2] <- yname
 ynorm.pred <- max(nonas[,1]) +1 - reg.out$linear.predictors
 reg.out$coefficients <- reg.out$coefficients * (-1)
 reg.out$coefficients[1] <- max(nonas[,1]) +1 + reg.out$coefficients[1]   # coeffs in orig scale
 ynorm.resi <- nonas[,1] - ynorm.pred
 reg.out$linear.predictors <- ynorm.pred
 reg.out$resids <- ynorm.resi

 vtext<- paste("Quantiles of", yname, "residuals")

 if (verbose == TRUE) {
   testnorm <- gofTestCensored(ynorm.resi,nonas[,2])
   ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
   qqPlotCensored(ynorm.resi, nonas[,2], add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = "Normal Q-Q Plot of residuals")
   mtext(ptext)
 }
   }   # end of flipping
     } # end of Y in original units

  if (is.data.frame(x.vars))  {             # multiple x variables.  Print r-squared.
    LRr2 <- signif(1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y)),4)
    McFr2 <- signif((1-reg.out$loglik[2]/reg.out$loglik[1]),4)
    Nag.r2 <- signif((1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))) /(1-exp(2*reg.out$loglik[1]/length(reg.out$y))),4)
    AIC <- -2*reg.out$loglik[2] + (2*reg.out$df +1)
    BIC <- -2*reg.out$loglik[2] + log(reg.out$df+reg.out$df.residual+1)*reg.out$df
    cat(" Likelihood R2 =", LRr2, "                   ", "AIC =", AIC,"\n")
    if (verbose == TRUE) cat(" Rescaled Likelihood R2 =", Nag.r2, "          ", "BIC =", BIC, "\n", "McFaddens R2 =", McFr2, "\n","\n")
  }
  else {                                    # 1 x variable.  Print correlation coefficients
    LRcorr <- signif(sign(reg.out$coefficients[2])*sqrt(1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))),4)
    McFcorr <- signif(sign(reg.out$coefficients[2])*sqrt(1-reg.out$loglik[2]/reg.out$loglik[1]),4)
    Nag.cor <- signif(sign(reg.out$coefficients[2])*sqrt((1-exp(-2*(reg.out$loglik[2]-reg.out$loglik[1])/length(reg.out$y))) /(1-exp(2*reg.out$loglik[1]/length(reg.out$y)))),4)
    AIC <- -2*reg.out$loglik[2] + (2*reg.out$df +1)
    BIC <- -2*reg.out$loglik[2] + log(reg.out$df+reg.out$df.residual+1)*reg.out$df
    cat(" Likelihood R =", LRcorr, "                   ", "AIC =", AIC,"\n", "Rescaled Likelihood R =", Nag.cor, "          ", "BIC =", BIC, "\n", "McFaddens R =", McFcorr, "\n","\n")
  }
  return(reg.out)                   # returns reg.out object if function assigned to object
}

#' Censored two-group permutation test
#'
#' @description Performs a permutation test of differences in means between two groups of censored data.
#' @param y1 The column of data values plus detection limits
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y1` column, and 0 (or `FALSE`) indicates a detected value in `y1`.
#' @param grp Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param R The number of permutations used. Default is 9999
#' @param alternative indicates the alternative hypothesis and must be one of "`two.sided`", "`greater`" or "`less`". You can specify just the initial letter. Default is "`two.sided`".
#' @keywords permutation difference test
#' @export
#' @return Permutation test results with the number of permutations, range in group means and their difference, and range in `p-value`.
#' @details Because this is a permutation test it avoids the problem with MLE tests (`cen2means`) that assume a normal distribution.  No values are modeled as below zero and `p-values` are trustworthy. Ranges in means and p-values are due to interval-censoring of censored data means.
#'
#' @references
#' Good, P., 2000. Permutation Tests: A Practical Guide to Resampling Methods for Testing Hypotheses, 2nd ed, Springer Series in Statistics. Springer-Verlag, New York, NY. <https://doi.org/10.1007/978-1-4757-3235-1>
#'
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#' data(PbHeron)
#' cenperm2(PbHeron$Liver,PbHeron$LiverCen,PbHeron$DosageGroup,alternative="t")

cenperm2 <- function(y1, y2, grp, R = 9999, alternative = "two.sided") {
  yname <- deparse(substitute(y1))
  gname <- deparse(substitute(grp))
  xdat <- na.omit(data.frame(y1, y2, grp))
  x1 <- xdat[,1]
  x2 <- xdat[,2]
  Factor <- as.factor(xdat[,3])
  df <- length(levels(Factor))-1
  grpname <- levels(Factor)

  x1.lo <- x1*(1-x2)
  permdiff.lo <- vector(length = R)
  permdiff.hi <- vector(length = R)

  mu1hi <- mean(x1[Factor == grpname[1]])
  mu1lo <- mean(x1.lo[Factor == grpname[1]])
  mu2hi <- mean(x1[Factor == grpname[2]])
  mu2lo <- mean(x1.lo[Factor == grpname[2]])
  n1 <- length(x1[Factor == grpname[1]])
  n2 <- length(x1[Factor == grpname[2]])
  n <- n1+n2
  dbarhi <- mu1hi - mu2hi
  dbarlo <- mu1lo - mu2lo

  for (i in 1:R) {newFact <- sample(Factor)  # replace = FALSE
  permdiff.lo[i] <- mean(x1.lo[newFact == grpname[1]]) - mean(x1.lo[newFact == grpname[2]])
  permdiff.hi[i] <- mean(x1[newFact == grpname[1]]) - mean(x1[newFact == grpname[2]])
  }
  absdiff.hi = abs(dbarhi)
  absperm.hi = abs(permdiff.hi)
  if(alternative == 't'| alternative=="two.sided") pval.hi <- (1+sum(absdiff.hi<=absperm.hi))/(R+1)
  if(alternative == 'g'| alternative=="greater")  pval.hi <- (sum(permdiff.hi >= dbarhi)+1)/(R+1)
  if(alternative == 'l'| alternative=="less") pval.hi = (sum(permdiff.hi <= dbarhi)+1)/(R+1)

  absdiff.lo = abs(dbarlo)
  absperm.lo = abs(permdiff.lo)
  if(alternative == 't'| alternative=="two.sided") pval.lo <- (1+sum(absdiff.lo<=absperm.lo))/(R+1)
  if (alternative == 'g'| alternative=="greater")  pval.lo <- (sum(permdiff.lo >= dbarlo)+1)/(R+1)
  if (alternative == 'l'| alternative=="less") pval.lo = (sum(permdiff.lo <= dbarlo)+1)/(R+1)

  diff.lo <- min(dbarlo, dbarhi)
  diff.hi <- max(dbarlo, dbarhi)
  p.lo <- min(pval.lo, pval.hi)
  p.hi <- max(pval.lo, pval.hi)
  result <- data.frame(diff.lo, diff.hi, p.lo, p.hi)

  #added by PJ
  if(alternative%in%c('t','g','l')){alternative <- switch(alternative,"t"="two.sided","g"="greater","l"="less")}
  #  write results
  mean.txt <- paste("     Mean (", grpname[1], " - ", grpname[2], ") = ", sep="")
  cat(" Permutation test of mean CensData:", yname, "  by Factor:", gname, '\n', "   ", R, "Permutations", "    alternative =", alternative, "\n")
  cat ("  mean of", grpname[1], "=", signif(mu1lo, 4), "to", signif(mu1hi,4), "    mean of", grpname[2], "=", signif(mu2lo,4), "to", signif(mu2hi,4), "\n")
  cat( mean.txt, signif(dbarlo, 4), "to", signif(dbarhi, 4), "      p =", pval.lo, "to", pval.hi, '\n', "\n")
  return(invisible(result))
}

#' Censored data one-factor permutation test
#'
#' @description Performs a permutation test of differences in means between groups of censored data.
#' @param y1 The column of data values plus detection limits.
#' @param y2 The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y1` column, and 0 (or `FALSE`) indicates a detected value in `y1`.
#' @param grp Grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param R The number of permutations used. Default is 9999.
#'
#' @importFrom NADA ros cenros
#' @export
#' @return Permutation test results with the number of permutations, range in test statistics and `p-value` values through the various permutations. Group means are also listed.
#' @details Because this is a permutation test it avoids the problem with MLE tests (cenanova) that assume a normal distribution.  No values are modeled as below zero and group means and `p-values` are trustworthy.
#'
#' @references
#' Good, P., 2000. Permutation Tests: A Practical Guide to Resampling Methods for Testing Hypotheses, 2nd ed, Springer Series in Statistics. Springer-Verlag, New York, NY. <https://doi.org/10.1007/978-1-4757-3235-1>
#'
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @examples
#'
#' data(PbHeron)
#' cenpermanova(PbHeron$Liver,PbHeron$LiverCen,PbHeron$DosageGroup)

cenpermanova <- function(y1, y2, grp, R = 9999) {
  yname <- deparse(substitute(y1))
  gname <- deparse(substitute(grp))
  xdat <- na.omit(data.frame(y1, y2, grp))
  x1 <- xdat[,1]
  x2 <- xdat[,2]
  Factor <- as.factor(xdat[,3])
  df <- length(levels(Factor))-1
  grpname <- levels(Factor)
  x1.lo <- x1*(1-x2)
  permW.lo <- vector(length = R)
  permW.hi <- vector(length = R)
  group.means <- vector(length = nlevels(Factor))
  mean.names <- vector(length = nlevels(Factor))
  gpname <- as.character(levels(Factor))

  # Test statistic Wprime is the between group sum of squares minus invariant terms.
  # see "Permutation Tests, Second Edition" by Good (2000), page 44.
  Wprime.groups <- function (x, Factor) {
    Wprime <- 0
    gpname <- as.character(levels(Factor))
    for (i in 1:nlevels(Factor)) {
      W.grp <- (sum(x[Factor == gpname[i]]))^2 / length(x[Factor == gpname[i]])
      Wprime <- Wprime + W.grp
    }
    return(invisible(Wprime))
  }

  # Test statistics on original data
  teststat.lo <- Wprime.groups(x1.lo, Factor)
  teststat.hi <- Wprime.groups(x1, Factor)
  absW.hi = abs(teststat.hi)
  absW.lo = abs(teststat.lo)
  # permutations
  for (i in 1:R) {newFact <- sample(Factor)  # replace = FALSE
  permW.lo[i] <- Wprime.groups(x1.lo, newFact)
  permW.hi[i] <- Wprime.groups(x1, newFact)
  }
  # p values always two-sided
  absperm.lo = abs(permW.lo)
  pval.lo <- (1+sum(absW.lo<=absperm.lo))/(R+1)
  absperm.hi = abs(permW.hi)
  pval.hi <- (1+sum(absW.hi<=absperm.hi))/(R+1)

  # group means
  for (i in 1:nlevels(Factor)) {
    ros.out <- suppressWarnings(cenros(x1[Factor == gpname[i]], as.logical(x2[Factor == gpname[i]])))
    group.means[i] <- signif(mean(ros.out),4)
    mean.names[i] <- paste("mean(", gpname[i], ")", sep="")
  }
  names(group.means) <- mean.names

  #  write test results
  p.lo <- min(pval.lo, pval.hi)
  p.hi <- max(pval.lo, pval.hi)
  result <- data.frame(teststat.lo, teststat.hi, p.lo, p.hi, group.means)
  cat(" Permutation test of mean CensData:", yname, "  by Factor:", gname, '\n', "   ", R, "Permutations", "\n")
  cat( "Test Statistic =", signif(teststat.lo, 4), "to", signif(teststat.hi, 4), "      p =", pval.lo, "to", pval.hi, '\n', "\n")
  print(group.means, row.names = FALSE, print.gap = 3)
  return(invisible(result))
}

#' Prediction interval for censored data
#'
#' @description Computes prediction intervals for censored data assuming lognormal, gamma and normal distributions.
#' @param y.var The column of y (response variable) detected values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param pi.type Designation of either a `“two-sided”` interval (default) or a 1-sided `“upper”` or 1-sided `“lower”` interval.
#' @param conf Confidence coefficient of the interval, 0.95 (default).
#' @param newobs The number of new observations to be contained in the interval.
#' @param method Character string specifying the method of estimation. Default is `mle` (maximum likelihood). See details.
#' @keywords prediction interval
#' @export
#' @importFrom EnvStats elnormCensored predIntLnorm enormCensored predIntNorm
#' @return A table of prediction limits based on user provided confidence coefficient (`conf`) and prediction invterval type (`pi.type`)
#' @details Computes prediction intervals for three distributions.  This is a front-end to the individual functions from the EnvStats package.  By default all three are computed using maximum likelihood estimation (mle). The gamma distribution for censored data uses the Wilson-Hilferty approximation (normal distribution on cube roots of data). Other methods are available in EnvStats, but few methods are available for all three distributions. For info on other methods, see help for elnormCensored and enormCensored commands in EnvStats.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Millard, S.P., 2013. EnvStats: An R Package for Environmental Statistics. Springer-Verlag, New York.
#'
#' Krishnamoorthy, K., Mathew, T., Mukherjee, S., 2008. Normal-Based Methods for a Gamma Distribution, Technometrics, 50, 69-78.
#'
#' @seealso [EnvStats::enormCensored]
#'
#' @examples
#'
#' data(PbHeron)
#'
#' # Default
#' cenPredInt(PbHeron$Liver,PbHeron$LiverCen)
#'
#' # User defined confidence coefficient
#' cenPredInt(PbHeron$Liver,PbHeron$LiverCen, conf=0.5)
#'
#' # User defined confidence coefficient outside of acceptable range
#' # the procedure will stop and give an error.
#' \dontrun{cenPredInt(PbHeron$Liver,PbHeron$LiverCen, conf=1.1)}
#'
#' # User defined prediction interval type
#' cenPredInt(PbHeron$Liver,PbHeron$LiverCen,pi.type="lower")
#' cenPredInt(PbHeron$Liver,PbHeron$LiverCen,pi.type="upper")


cenPredInt <- function(y.var, cen.var, pi.type = "two-sided", conf = 0.95, newobs = 1, method = "mle")  {
  if(conf>1|conf<0){stop("Please select a confidence coefficient between 0 and 1.")}

  obj.lnorm <- elnormCensored (y.var, cen.var, method = method)
  obj.lnorm2 <- predIntLnorm (obj.lnorm, k=newobs, pi.type = pi.type, conf.level = conf)
  obj.norm <- enormCensored (y.var, cen.var, method = method)
  obj.norm2 <- predIntNorm (obj.norm, k=newobs, pi.type = pi.type, conf.level = conf)
  dat.gamma <- y.var^(1/3)
  obj.gamma <- enormCensored (dat.gamma, cen.var, method = method)
  obj.gamma2 <- predIntNorm (obj.gamma, k=newobs, pi.type = pi.type, conf.level = conf)
  pi1.text <- paste(conf*100, "% LPL", sep="")
  pi2.text <- paste (conf*100, "% UPL", sep="")

  if (pi.type == "two-sided") {
    title.txt <- paste (conf*100, "% Prediction Limits", sep ="", "\n")

    pi.lnorm <- obj.lnorm2$interval$limits
    pi.norm <- obj.norm2$interval$limits
    pi.gamma <- (obj.gamma2$interval$limits)^3
    pi.gamma[1] <- max(0, pi.gamma[1])
  }
  else if (pi.type == "lower") {  title.txt <- paste (conf*100, "% Lower Prediction Limit", sep ="", "\n")
  pi.lnorm <- obj.lnorm2$interval$limits[1]
  pi.norm <- obj.norm2$interval$limits[1]
  pi.gamma <- max(0, (obj.gamma2$interval$limits[1])^3 )
  }
  else {   title.txt <- paste (conf*100, "% Upper Prediction Limit", sep ="", "\n")
  pi.lnorm <- obj.lnorm2$interval$limits[2]
  pi.norm <- obj.norm2$interval$limits[2]
  pi.gamma <- (obj.gamma2$interval$limits[2])^3
  }

  cat (title.txt)

  Distribution <- c("Lognormal", "Gamma", "Normal")
  lpl <- c(obj.lnorm2$interval$limits[1], max(0, (obj.gamma2$interval$limits)[1]^3), obj.norm2$interval$limits[1])
  if (pi.type == "upper") {lpl[1:3] <- "NA"}
  upl <- c(obj.lnorm2$interval$limits[2], obj.gamma2$interval$limits[2]^3, obj.norm2$interval$limits[2])
  if (pi.type == "lower") {upl[1:3] <- "NA"}
  results <- data.frame(Distribution, lpl, upl)
  names(results) <- c("Distribution", pi1.text, pi2.text)
  return(results)
}


#' Q-Q Plot censored data
#'
#' @description Plots quantile-quantile (Q-Q) plot of censored data fitted a data distribution
#' @param y.var The column of `y` (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param dist One of three distributional shapes to fit to your data:  lognormal (`lnorm`), normal (`norm`) or gamma (`gamma`).
#' @param Yname Optional – input text in quotes to be used as the variable name on the Q-Q plot.  The default is the name of the `y.var` input variable.
#' @export
#'
#' @return A single Q-Q plot of data fitted to normal, lognormal or gamma distributions with Shapiro-Francia W value printed on plot.
#'
#' @importFrom EnvStats gofTestCensored qqPlotCensored distChooseCensored
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#'\dontrun{
#' data(Brumbaugh)
#' cenQQ(Brumbaugh$Hg,Brumbaugh$HgCen)
#'
#' # User defined distribution
#' cenQQ(Brumbaugh$Hg,Brumbaugh$HgCen,dist="gamma")
#'}



cenQQ <- function(y.var, cen.var, dist = "lnorm", Yname = yname)  {
  #added to stop if dist is not from the list
  #if(!(dist%in%c("norm","lnorm","gamma"))){stop(paste0(dist," distribution is not supported with this function, try again."))}

  yname <- deparse(substitute(y.var))
  cen.logical <- as.logical(cen.var)

  if (sum(as.integer(cen.var)) > 0)    # not all data are detects
  {var.choose <- distChooseCensored(y.var, cen.logical)
  norm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$norm$statistic, 3) )
  lnorm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$lnorm$statistic, 3) )
  gamma.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$gamma$statistic, 3) )

  if (dist == "norm") {
    qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
    mtext(norm.text)
    #  legend("bottomright", legend = norm.text)
  }

  if (dist == "lnorm")  {
    ylabel <- paste ("ln (", Yname, ")", sep = "")
    qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = ylabel, main = "Lognormal Q-Q Plot")
    mtext(lnorm.text)
    #  legend("bottomright", legend = lnorm.text)
  }

  if (dist == "gamma")  {
    qqPlotCensored(y.var, cen.logical, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
    mtext(gamma.text)
    #   legend("bottomright", legend = gamma.text)
  }
  }

  else    # all data are detects
  {var.choose <- distChoose(y.var, method = "sf", alpha = 0.05, choices = c("norm", "gamma", "lnorm"))
  norm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$norm$statistic, 3) )
  lnorm.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$lnorm$statistic, 3) )
  gamma.text <- paste("Shapiro-Francia W =", signif(var.choose$test.results$gamma$statistic, 3) )

  if (dist == "norm") {
    EnvStats::qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", xlab = "Normal Quantiles", ylab = Yname, main = "Normal Q-Q Plot")
    mtext(norm.text)
  }

  if (dist == "lnorm")  {
    ylabel <- paste ("ln (", Yname, ")", sep = "")
    EnvStats::qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", distribution = "lnorm", xlab = "Normal Quantiles", ylab = ylabel, main = "Lognormal Q-Q Plot")
    mtext(lnorm.text)
  }

  if (dist == "gamma")  {
    EnvStats::qqPlot(y.var, pch = 19, add.line = TRUE, line.col = "red", distribution = "gamma", estimate.params = TRUE, ylab = Yname, main = "Gamma Q-Q Plot")
    mtext(gamma.text)
  } }
}

#' Q-Q plot of censored regression residuals
#'
#' @description Plots a quantile-quantile (Q-Q) plot of censored regression residuals for simple or multiple regression.
#' @param y.var The column of `y` (response variable) values plus detection limits. Alternative, with interval-censord data, the column of the lower end of the interval.
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.  Alternatively, with interval-censored data the column of the high end of the interval.
#' @param x.vars One or more uncensored explanatory variable(s). See Details
#' @param LOG Indicator of whether to compute the regression in the original y units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @param intcens a logical value indicating the input data is interval-censored instead of a column of values plus a column of indicators.
#' @param main overal ltitle for the plot
#' @export
#' @return Q-Q Plot of model residuals and Shapiro-Francia test results.
#'
#' @importFrom survival survreg Surv
#' @importFrom EnvStats qqPlotCensored gofTestCensored
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Millard, S.P., 2013. EnvStats: An R Package for Environmental Statistics. Springer-Verlag, New York.
#'
#' Shapiro, S.S., Francia, R.S., 1972. An approximate analysis of variance test for normality. Journal of the American Statistical Association 67, 215–216.
#'
#' @examples
#' data(Brumbaugh)
#'
#' # One variable
#' cenregQQ(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh$PctWetland)
#'
#' # More than one variable for demostration purposes
#'cenregQQ(Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh[,c("PctWetland","SedLOI","Weight")])

cenregQQ <- function(y.var, cen.var, x.vars, LOG = TRUE, intcens = FALSE, main = NULL) {
  yname <- deparse(substitute(y.var))
  if (LOG == TRUE)  {lnvar <- log(y.var)
  flip.log <- max(lnvar) +1 - lnvar
  surv.log <- Surv(flip.log, as.logical(1-cen.var), type="right" )
  if (is.data.frame(x.vars))  {
    reg.out <- survreg(surv.log ~ ., data = x.vars, dist = "gaussian")

  }
  else {reg.out <- survreg(surv.log~x.vars, dist = "gaussian") }    # 1 x variable
  newcoeffs <- reg.out$coefficients * (-1)
  newcoeffs[1] <- max(lnvar)+1 + newcoeffs[1]
  ylog.pred <- max(lnvar) +1 - reg.out$linear.predictors
  ylog.resi <- lnvar - ylog.pred
  vtext<- paste("Quantiles of", yname, "residuals (log units)")
  testnorm <- gofTestCensored(ylog.resi,cen.var)
  if (is.null(main)) main <- "Lognormal Q-Q Plot of Residuals"
  ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic,5), "  p =", round(testnorm$p.value,6))
  qqPlotCensored(ylog.resi, cen.var, add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = main)
  mtext(ptext)
  }
  else
  {
    if (intcens == TRUE) {y.low <- y.var*(1-cen.var)
    surv.norm <- Surv(y.low, y.var, type="interval2")}
    else {flip.norm <- max(y.var) +1 - y.var
    surv.norm <- Surv(flip.norm, as.logical(1-cen.var) )}

    if (is.data.frame(x.vars))  {
      reg.out <- survreg(surv.norm ~ ., data = x.vars, dist = "gaussian")
    }

    else {reg.out <- survreg(surv.norm~x.vars, dist = "gaussian") }
    if (intcens == TRUE) {ynorm.pred <- reg.out$linear.predictors
    ynorm.resi <- y.low - ynorm.pred}
    else {ynorm.pred <- max(y.var) +1 - reg.out$linear.predictors
    ynorm.resi <- y.var - ynorm.pred}
    vtext<- paste("Quantiles of", yname, "residuals")
    testnorm <- gofTestCensored(ynorm.resi,cen.var)
    if (is.null(main)) main <- "Normal Q-Q Plot of Residuals"
    ptext <- paste("Shapiro-Francia W =", round(testnorm$statistic, 5), "  p =", round(testnorm$p.value, 5))
    qqPlotCensored(ynorm.resi, cen.var, add.line = T, prob.method = "modified kaplan-meier", ylab = vtext, main = main)
    mtext(ptext)
  }
}

#' Seasonal Kendall permutation test on censored data
#'
#' @param time Column of the time variable, either a sequence of days or decimal times, etc.  Will be the scale used for time in the trend analysis.
#' @param y The column of y (response variable) values plus detection limits
#' @param y.cen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param group Column of the season classifications. A factor in R, so usually though not necessarily a text variable.  If numeric, define as a factor before running the script.
#' @param LOG Indicator of whether to compute the regression in the original y units, or on their logarithms.  The default is to use the logarithms (LOG = `TRUE`).  To compute in original units, specify the option LOG = `FALSE` (or LOG = 0).
#' @param R The number of repetitions in the permutation process.  R is often between 999 and 9999 (+1 for the observed test statistic produces 1000 to 10000 repetitions). By default R=4999. Increasing R simply results in lower variation in the pvalues produced between runs.
#' @param nmin The minimum number of observations needed for the entire time period to be tested, per season.  For example, with 1 sample per year per season over an 8-year period, you have 8 observations for each season.  You can increase this number if you want a higher minimum.  Don’t decrease it below 4.  If there are fewer than nmin values that season is skipped and not included in the overall test & a note will be printed.
#' @param seaplots In addition to the plot of the overall Seasonal Kendall trend line, plots of the trend in individual seasons can also be drawn.
#'
#' @return Prints the Kendall trend test results for each season individually. The overall Seasonal Kendall test and Theil-Sen line results are both printed and returned.
#'
#' If `seaplots=TRUE` each season's trend line will be plotted along with the overall Seasonal Kendall (Akritas-Theil-Sen) line.
#' If `seaplots=FALSE` only the overall Seasonal Kendall (Akritas-Theil-Sen) line will be plotted on a data scatterplot.
#' @export
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Hirsch, R.M., Slack, J.R., Smith, R.A., 1982. Techniques of Trend Analysis for Monthly Water Quality Data, Water Res. Reseach 18, 107-121.
#'
#' @seealso [NADA::cenken]
#'
#' @examples
#'
#' data(Brumbaugh)
#'
#' #Artifical time and season variables for demonstration purposes
#' Brumbaugh$time=1:nrow(Brumbaugh)
#' Brumbaugh$sea=as.factor(round(runif(nrow(Brumbaugh),1,4),0))
#'
#' with(Brumbaugh,censeaken(time,Hg,HgCen,sea,seaplots = TRUE))


censeaken <- function(time, y, y.cen, group, LOG = FALSE, R = 4999, nmin = 4, seaplots = FALSE)
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
    z <- data.frame(x.seas, y.seas)
    kenplot(yc[[i]], cc[[i]], xc[[i]], xcen = rep(0, times=ntest), xnam=xname, ynam=yname, atsline = TRUE)
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

  kenplot(yyy, ccc, xxx, xcen = rep(0, times=nall), xnam=xname, ynam=yname, Title = "Seasonal Kendall Test")
  abline(intall, medslope, lwd=2, col = "blue")
  mtext ("Overall Trend Line", col = "blue")

  return (invisible(RESULTS))
}

#' Upper Tolerance interval for censored data
#'
#' @description Computes a one-sided upper tolerance interval for censored data assuming lognormal, gamma and normal distributions.
#' @param y.var The column of y (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detectionlimit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param conf Confidence coefficient of the interval, 0.95 (default).
#' @param cover Coverage, the percentile probability above which the tolerance interval is computed.  The default is 90, so a tolerance interval will be computed above the 90th percentile of the data.
#' @param method.fit The method used to compute the parameters of the distribution.  The default is maximum likelihood (`“mle”`). The alternative is robust ROS (`“rROS”`).  See Details.
#' @importFrom EnvStats elnormCensored predIntLnorm enormCensored predIntNorm eqlnormCensored eqnormCensored
#' @importFrom fitdistrplus fitdistcens
#' @export
#'
#' @return  Prints and returns the percentile (`cover`), upper tolerance limit (`conf`) and BIC of fit for lognormal, normal and approximated gamma distributions. Plots empirical and theoretical CDFs with BIC values as legend.
#' @details Computes upper one-sided tolerance intervals for three distributions.  This is a front-end to the individual functions from the EnvStats package.  By default all three are computed using maximum likelihood estimation (mle); robust ROS is available as an alternate method for all three distributions. The gamma distribution for censored data uses the Wilson-Hilferty approximation (normal distribution on cube roots of data). For more info on the relative merits of robust ROS versus mle, see Helsel (2011) and Millard (2013).
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Millard, S.P., 2013. EnvStats: An R Package for Environmental Statistics. Springer-Verlag, New York.
#'
#' Krishnamoorthy, K., Mathew, T., Mukherjee, S., 2008. Normal-Based Methods for a Gamma Distribution, Technometrics, 50, 69-78.
#'
#' @examples
#'
#' data(PbHeron)
#'
#' # Default
#' cenTolInt(PbHeron$Liver,PbHeron$LiverCen)
#'
#' # User defined conficence interval
#' cenTolInt(PbHeron$Liver,PbHeron$LiverCen,conf=0.75)
#'
#' # User defined percentile
#' cenTolInt(PbHeron$Liver,PbHeron$LiverCen,cover=0.5)
#'
#' # inputs outside acceptable ranges
#' \dontrun{
#' cenTolInt(PbHeron$Liver,PbHeron$LiverCen,cover=1.25)
#' cenTolInt(PbHeron$Liver,PbHeron$LiverCen,conf=1.1)
#' cenTolInt(PbHeron$Liver,PbHeron$LiverCen,method.fit="ROS")
#' }

cenTolInt <- function(y.var, cen.var, conf = 0.95, cover = 0.9, method.fit = "mle")
{
  if(conf>1|conf<0){stop("Confidence coefficient should be between 0 and 1, please try again.")}
  if(cover>1|cover<0){stop("Percentile should be between 0 and 1, please try again.")}
  if(!(method.fit%in%c("mle","rROS"))){stop("Please select 'mle' or 'rROS' for method fit.")}

  nameofy <- deparse(substitute(y.var))
  obj.lnorm <-  eqlnormCensored (y.var, cen.var, p=cover, method = method.fit, ci=TRUE, ci.type = "upper", conf.level = conf)
  obj.norm <- eqnormCensored (y.var, cen.var, p=cover, method = method.fit, ci=TRUE, ci.type = "upper", conf.level = conf)
  dat.gamma <- y.var^(1/3)
  obj.gamma <- eqnormCensored (dat.gamma, cen.var, p=cover, method = method.fit, ci=TRUE, ci.type = "upper", conf.level = conf)

lnorm.txt <- paste ("Lognormal ", cover*100,"th Pctl", "      ", conf*100, "% Upper Tolerance Limit", sep ="", "\n")
norm.txt <-  paste ("Normal ", cover*100,"th Pctl", "      ", conf*100, "% Upper Tolerance Limit", sep ="", "\n")
gamma.txt <-  paste ("~Gamma ", cover*100,"th Pctl", "     ", conf*100, "% Upper Tolerance Limit", sep ="", "\n")
pct.text <- paste(cover*100, "th Pctl", sep="")
ti.text <- paste (conf*100, "% UTL", sep="")

pct.lnorm <- obj.lnorm$quantiles
ti.lnorm <- obj.lnorm$interval$limits[2]
pct.norm <- obj.norm$quantiles
ti.norm <- obj.norm$interval$limits[2]
pct.gamma <- obj.gamma$quantiles^3
ti.gamma <- (obj.gamma$interval$limits[2])^3

left <- y.var*(1-as.integer(cen.var))
right <- y.var
var.frame <- data.frame(left, right)
y.lnorm <- fitdistcens(var.frame, "lnorm")
y.gamma <- fitdistcens(var.frame, "gamma")
y.norm <- fitdistcens(var.frame, "norm")
bic.3 = c(y.lnorm$bic, y.gamma$bic, y.norm$bic)
ti.best = c(ti.lnorm, ti.gamma, ti.norm)
pct.best = c(pct.lnorm, pct.gamma, pct.norm)

Distribution <- c("Lognormal", "Gamma", "Normal")
results <- data.frame(Distribution, pct.best, ti.best, bic.3, method.fit)
names(results) <- c("Distribution", pct.text, ti.text, "BIC", "Method")

bic.best = min(bic.3)
dist.best = floor(bic.best/bic.3)
pct.best = max(pct.best*dist.best)
ti.best = max(ti.best*dist.best)

cenCompareCdfs (y.var, cen.var, Yname = nameofy, dist3 = "norm")
abline (h=cover, lty = "dotted", col = "black")
return (results)
}

#' Trend analysis of censored data with a covariate
#'
#' @description Trend analysis after adjustment of censored data for a covariate.
#' @param y.var The column of y (response variable) values plus detection limits
#' @param y.cens The column of indicators, where 1 (or `TRUE`) indicates a detectionlimit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param x.var Column of a covariate (not time).  `y.var` will be smoothed versus `x.var` and residuals taken to subtract out the relationship between `y` and `x`.
#' @param time.var Column of the time variable, either a sequence of days or decimal times, etc.  Will be the scale used for time in ATS trend analysis.
#' @param link Default = `“identity”` which means it uses data in the original units. See details.
#' @param Smooth Type of smoother used in the GAM. Default is `“cs”`, shrinkage cubic regression splines. See details for other options.
#' @keywords trend analysis GAM spline
#' @export
#'
#' @importFrom mgcv gam
#' @importFrom NADA cenxyplot
#' @importFrom cenGAM tobit1
#' @return
#'
#' Prints three plots: Y data vs time with GAM Smooth, Residuals from GAM Smooth vs time, and ATS trend line of residuals vs time.
#'
#' Returns GAM residuals and ATS results on trend test of residuals (intercept, slope, Kendall's tau, p-value for trend)
#'
#' @details
#'
#' Default `link` = identity. Other options are available see `cenGAM::tobit1` for more options.
#'
#' Default `Smooth` is `"cs"` for shrinkage cubic regression splines. See `mgcv::smooth.terms` for other types of smoothing algorithms.  '"ts"' is a thin-plate regression spline and is also commonly used.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#'@seealso [mgcv::gam]
#'
#' @examples
#'
#' data(Brumbaugh)
#'
#' Brumbaugh$time=1:nrow(Brumbaugh)
#'
#' with(Brumbaugh,centrend(Hg,HgCen,SedTotHg,time.var=time))


centrend <- function(y.var, y.cens, x.var, time.var, link = "identity", Smooth = "cs") {
  yname <- deparse(substitute(y.var))
  xname <- deparse(substitute(x.var))
  tname <- deparse(substitute(time.var))
  y.txt <- paste(yname, "residuals")
  dl = y.cens * y.var * 1.001
  #  y.cens <- as.factor(y.cens)
  dat.all <- data.frame(y.var, x.var, time.var, dl, y.cens)
  dat.nonas <- na.omit(dat.all)
  colnames(dat.nonas) <- c(yname, xname, tname, "DL", "Cens = 1")

  # cs, ts, tp are probably best choices.
  gam.y <- gam(dat.nonas[,1] ~ s(dat.nonas[,2], bs = Smooth), family = tobit1(link = link, left.threshold = dat.nonas[,4]))
  # sorting fitted values
  o <- order(dat.nonas[,2] , dat.nonas [,1])
  #plot(dat.nonas[,1] ~ dat.nonas[,2], main = "1. Data and GAM Smooth", ylab = yname, xlab = xname)

  oldpar<- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(3,1))

   # plot 1.  Y data vs covariate, with smooth
  x.cens = rep(0, times=length(dat.nonas[,3]) )
  cenxyplot(dat.nonas[,2], as.logical(x.cens), dat.nonas[,1], as.logical(dat.nonas[,5]), main = "1. Data and GAM Smooth", ylab = yname, xlab = xname, pch = 19, cex = 0.7)
  # draw the smooth
  lines (dat.nonas[o,2], gam.y$fitted.values[o], col = 'red')

  # plot 2. Y residuals from smooth vs. covariate
  plot(gam.y$residuals ~ dat.nonas[,2], main = "2. Residuals from GAM Smooth", xlab = xname, ylab = y.txt)
  abline (h=0, col = "blue")
  cat("Trend analysis of", yname, "adjusted for", xname, "\n")

  # plot 3.  Y residuals vs time with ATS line
  ats.out <- ATS(gam.y$residuals, dat.nonas[,5], dat.nonas[,3], x.cens, LOG = FALSE, xlabel = tname, ylabel = y.txt)
  dat.out <- data.frame(gam.y$residuals, dat.nonas[,2])
  colnames(dat.out) <- c("GAMresidual", xname)
  dat.out <- cbind(dat.out, ats.out)
  return (invisible(dat.out))
}

#' Compute an ECDF and Distribution Parameters for Censored Data
#'
#' @description Computes the empirical cumulative distribution function (ECDF) for censored data. Estimates parameters of the distribution, including the mean and quantiles.
#' @param y1 Either the lowest possible concentrations (interval-censored format) or conccentrations plus detection limits for indicator formated data.
#' @param y2 Either the highest possible concentrations (interval-censored format) or censoring indicators (logical. 1 or `TRUE` = censored, 0 or FALSE = detected) for indicator formatted data.
#' @param conf The confidence coefficient for confidence intervals around the Kaplan-Meier mean and median. Default = 0.95.
#' @param qtls Probabilities for the quantiles to be estimated.  Defaults are (0.10, 0.25, 0.50, 0.75, 0.90).  You may add and/or substitute probabilities -- all must be between and not including 0 to 1..
#' @param Cdf Logical `TRUE`/`FALSE` indicator of whether to plot the empirical cumulative distribution function (cdf) using Kaplan-Meier quantiles.
#' @param printstats Logical `TRUE`/`FALSE` option of whether to print the resulting statisics in the console window, or not.  Default is TRUE.
#' @param Ylab Optional input text in quotes to be used as the variable name on the ecdf plot.  The default is the name of the `y1` input variable.
#'
#' @importFrom survival Surv survfit
#' @importFrom stats quantile
#' @return
#' If `printstats=TRUE`: Based on the provided `conf` value, Kaplan-Meier summary statistics (`mean`,`sd`,`median`), lower and upper confidence intervals around the mean and median value, sample size and percent of censored samples are returned. The specified quantile values are also printed and returned.
#'
#' If `Cdf=TRUE`: The ecdf of censored data is plotted.
#' @details Quantiles and parameters are estimated using the survfit function. The mean computed is the "restricted mean" (see help for the survfit function of the survival package).  Internally uses interval-censoring to avoid a small bias in the mean produced by the NADA package's cenfit function, which uses the reverse Kaplan-Meier procedure, converting left-censored to right-censored data prior to computing the ecdf and mean. See Gillespie et al. for more discussion.
#'
#' @export
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Gillespie, B.W., et al., 2010.  Estimating Population Distributions When Some Data Are Below a Limit of Detection by Using a Reverse Kaplan-Meier Estimator. Epidemiology 21, 564-570.
#'
#' @seealso [survival::survfit] [NADA::cenfit]
#' @examples
#'
#' data(Brumbaugh)
#'
#' cfit(Brumbaugh$Hg,Brumbaugh$HgCen)
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

#' Kendall's S-statistic for permutations of censored data
#'
#' @description Computes a Kendall rank correlation S-statistic for permutations of censored data. Collectively these represent the variation in S expected when the null hypothesis is true.  Called by censeaken. computeS is not expected to be of much use to users on its own.
#' @param x Column of the time variable, either a sequence of days or decimal times, etc.  Time data for one season.
#' @param y The column of y (response variable) values plus detection limits for one season.
#' @param ycen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y` column, and 0 (or `FALSE`) indicates a detected value in `y`.
#' @param seas Name of a single season classification. Usually though not necessarily a text variable.
#' @param R The number of repetitions in the permutation process.  R is often between 999 and 9999 (+ the 1 observed test statistic produces 1000 to 10000 realizations).
#'
#' @return An Rx1 matrix containing an S-value for each of the R data permutations.
#' @seealso [Kendall::Kendall]
#' @references
#' Helsel, D.R., Hirsch, R.M., Ryberg, K.R., Archfield, S.A., Gilroy, E.J., 2020. Statistical Methods in Water Resources. U.S. Geological Survey Techniques and Methods, book 4, chapter A3, 458p., https://doi.org/10.3133/tm4a3.

#' @examples
#' data(Brumbaugh)
#'
#' #Artifical time and season variables for demonstration purposes
#' Brumbaugh$time=1:nrow(Brumbaugh)
#' Brumbaugh$sea=as.factor(round(runif(nrow(Brumbaugh),1,4),0))
#'
#'
#' with(Brumbaugh,computeS(time,Hg,HgCen,sea,R=100))
#'
computeS<- function(x, y, ycen, seas = NULL, R=R) {
  #  called by censeaken.  Returns permutations of S for a single season
  nx = length(x); nj <- nx-1
  Sout <- matrix(0, nrow =R, ncol=1)
  delmax <- 0
  delmin <-0
  delsign<-0
  for (i in 1:R) {xperm <- sample(x)
  yperm <- y[order(xperm)]
  ycenperm <- ycen[order(xperm)]
  Ymin <- yperm*(1-ycenperm)   # 0 for NDs
  Ymax <- Ymin+(ycenperm*0.99*yperm)  # <1 becomes 0.99
  Sindiv <- rep(0, nx*(nj)/2)
  ns <- 1
  for (j in 1:nj) {
    for (k in (j+1):nx) {
      delmax <- sign(Ymax[k] - Ymax[j])
      delmin <- sign(Ymin[k] - Ymin[j])
      delsign <- abs(sign(delmax - delmin)) # if 0, delmax is correct.  if 1 Sindiv =0
      Sindiv[ns] <- (1-delsign)*delmax
      ns <- ns+1
    }
  }
  Sout[i,1] <- sum(Sindiv)
  }
  # print(Sout)
  return (Sout)
}

#' Censored data sample size
#'
#' @description Computes the equivalent sample size of censored data.  Observations at lower detection limits have a greater percent of the equivalent information of a detected value than observations at higher detection limits.
#' @param y.var The column of data values plus detection limits.
#' @param y.cen The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @keywords Sample Size censored
#' @export
#' @importFrom NADA censummary
#' @details
#' Based on "Method 2" of Dr. Brenda Gillespie's talk at ASA National Meeting 2019.  This method differs from hers in how the percentile probabilities for the detection limits are computed.  Probabilities here are computed using Regression on Order Statistics (ROS).
#'
#' Computes the equivalent n, the number of observations including censored values, as a measure of information content for data with nondetects.
#'
#' @return Prints summary statistics including
#' \itemize{
#' \item `n` sample size
#' \item `n.cen` number of censored data
#' \item `pct.cen` percent of data censored
#' \item `min` minimum reported value
#' \item `max` maximum reported value
#' }
#'
#' Summary of censored data including
#' \itemize{
#' \item `limit` detection limit
#' \item `n` number of censored values per limit
#' \item `uncen` number of detected values at or above the limit
#' \item `pexceed` proportion of data that exceeds the limit
#' }
#'
#' Summary of the equivalent sample size for detected and censored values.
#' \item `n.equiv` the equivalent number of observations
#' \item `n.cen.equiv` equivalent number of detected obs in the censored data
#' \item `n.detected` number of uncensored values
#'}
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Gillespie, B.W., Dominguez, A., Li, Y., 2019. Quantifying the information in values below the detection limit (left-censored data).  Presented at the 2019 Joint Statistical Meetings of the Amer. Stat. Assoc., Denver, CO., July 31, 2019.
#'
#' @seealso [NADA::censummary]
#'
#' @examples
#' data(Brumbaugh)
#'
#' equivalent_n(Brumbaugh$Hg,Brumbaugh$HgCen)

equivalent_n <- function(y.var, y.cen){

  ycen <- as.logical(y.cen)
  yname <- deparse(substitute(y.var))

aa <- censummary(y.var, ycen)
n.equiv <- sum(aa$limits[,2]*aa$limits[,4] + aa$limits[,3])
n_total <- sum(aa$limits[,2] + aa$limits[,3])
n_detected <- sum(aa$limits[,3])
n_cens <- sum(aa$limits[,2])
n.cen.equiv <- sum(aa$limits[,2]*aa$limits[,4])
n.detected <- n.equiv - n.cen.equiv
message(yname); print(aa)
equiv.out <- data.frame(n.equiv,  n.cen.equiv, n.detected)
cat("equivalent sample size:", "\n")
print(equiv.out, row.names = FALSE, print.gap = 3)
aa[["equivalent"]] <- equiv.out
return (invisible (aa))
}

#' Plot robust median ATS line for censored data
#'
#' @description
#' Function used by other functions to plot Akritas-Theil-Sen (ATS) line for censored data.  For one x variable regression. Both Y and X variables may be censored.
#' @param y1 The column of y (response variable) values plus detection limits
#' @param ycen The y-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and `0` (or `FALSE`) indicates a detected value in `y.var`.
#' @param x1 The column of x (explanatory variable) values plus detection limits
#' @param xcen The x-variable indicators, where 1 (or `TRUE`) indicates a detection limit in the `x.var` column, and `0` (or `FALSE`) indicates a detected value in `x.var`.
#' @param atsline Indicator of whether to draw the ATS line or not. Default is FALSE.
#' @param xnam Custom label for the x axis of plots.  Default is x variable column name.
#' @param ynam Custom label for the y axis of plots.  Default is y variable column name.
#' @param Title Custom title for plots.  Default is "Akritas - Theil - Sen line".
#' @importFrom NADA cenken
#' @return
#' Scatterplot of data plus ATS line.  Censored values are drawn for both X and Y variables as dashed lines up to the detection limits.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @seealso [NADA::cenken]
#'
#' @examples
#' \dontrun{
#' # Both y and x are censored
#' data(PbHeron)
#' with(PbHeron, kenplot(Blood, BloodCen, Kidney, KidneyCen))
#'
#' # x is not censored
#' data(Brumbaugh)
#' with(Brumbaugh, kenplot(Hg, HgCen, PctWetland,rep(0, times=length(PctWetland))))
#' }

kenplot <- function(y1, ycen, x1, xcen, atsline = FALSE, xnam=NULL, ynam=NULL, Title="Akritas - Theil - Sen line")  {
  alldat <- data.frame(y1, ycen, x1, xcen)

  if (!is.null(xnam)) {xnam = xnam}else(xnam=deparse(substitute(x1)))
  if (!is.null(ynam)) {ynam = ynam}else(ynam=deparse(substitute(y1)))

  xmin <- min(alldat[,3])
  xmax <- max(alldat[,3])
  ymin <- min(alldat[,1])
  ymax <- max(alldat[,1])

  if (ymin <= 0) {
    y1a <- y1+abs(min(y1))+1
    cka<- cenken(y1a, as.logical(ycen), x1, as.logical(xcen))
    int <- cka$intercept - (abs(min(y1))+1)
    slp <- cka$slope
    tau <- cka$tau
    pval <- cka$p
  }
  else { ck<- cenken(y1, as.logical(ycen), x1, as.logical(xcen))
  int <- ck$intercept;  slp <- ck$slope
  tau <-ck$tau;  pval<- ck$p
  }

  bothdetect <- as.integer(ycen)+as.integer(xcen)
  detected <- alldat[bothdetect == 0,]
  ynd <- alldat[as.integer(ycen) ==1 ,]  # set of censored y values
  nyc <- length(as.integer(ynd$ycen))    # number of censored y values
  xnd <- alldat[as.integer(xcen) ==1 ,]  # set of censored x values
  nxc <- length(as.integer(xnd$xcen))    # number of censored x values

#  oldpar <- par(no.readonly = TRUE)
#  on.exit(par(oldpar))

  plot(detected$x1, detected$y1, ylim = c(ymin, ymax), xlim = c(xmin, xmax), ylab = ynam, xlab = xnam, pch=19, cex=0.7, main=Title, xaxs="r", yaxs="r")
  if (atsline == TRUE) {
    abline(int, slp, col = "purple")}

  # vertial dashed lines for y censored
  if (nyc != 0) {
      for (i in 1:nyc ){
        dashy <- min(0-0.5, ymin)
        dashx <- ynd[i,3]
        dash <- data.frame (dashx, dashy)
        dash[2,1] <- ynd[i,3]
        dash[2,2] <- ynd[i,1]
        lines(dash, lty="dashed", col = "red")
      }
    }
    # horizontal dashed lines for x censored
    if (nxc != 0) {
      for (i in 1:nxc ){
        dashy <- xnd[i,1]
        dashx <- min(0, xmin-0.5)
        dash <- data.frame (dashx, dashy)
        dash[2,1] <- xnd[i,3]
        dash[2,2] <- xnd[i,1]
        lines(dash, lty="dashed", col = "red")
      }
    }
  }


#' Computes ranks of data with one or multiple detection limits
#'
#' @description Computes the within-column ranks of data having one or more detection limits. If multiple limits are present in a column, data are first re-censored at the highest detection detection limit.
#' @param dat.frame A data frame. Default format is paired = `TRUE`, where for 3 chemical parameters the input format is C1 I1 C2 I2 C3 I3, a concentration column followed by its censoring indicator column.
#' @param paired An option to specify paired = `FALSE`, where the input format would be C1 C2 C3 I1 I2 I3 where the C columns contain concentrations or a detection limit, and the I columns are their associated indicators, in the same order as the concentration columns.
#' @export
#' @return Returns columns of ranks of censored data in the same order as the paired columns of input data.  For 3 chemical parameters, the data frame returned will be R1 R2 R3 where R represents the ranks of the C1 C2 C3 input data accounting for the censoring indicated by columns I1 I2 I3.
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @examples
#' library(NADA) #For example data
#' data(PbHeron)
#'
#' ordranks(PbHeron[,4:15])
#'

ordranks <- function(dat.frame, paired = TRUE) {
  cols <- ncol(dat.frame)
  half <- cols/2
  dat.orig <- as.matrix(dat.frame)
  j=0;  ranks=0; DLmax=0; data.1 <- dat.orig; data.0 <- dat.orig
  if (paired) { for (i in seq(1, to=(cols-1), by=2)) {j = j+1
  data.0[,j] <- dat.orig[,i]
  data.0[,j+half] <- dat.orig[,i+1]
  nvec <- colnames(data.0)
  nvec <- nvec[seq(1, (cols-1), 2)]
  nvec <- paste("rnk.", nvec, sep="")
  } }
  else {data.0 <- dat.orig
  nvec <- colnames(data.0)
  nvec <- nvec[1:half]
  nvec <- paste("rnk.", nvec, sep="")
  }

  for (i in 1:half) {DLmax <- max(data.0[,i]*data.0[,i+half])
  data.1[,i] <- data.0[,i]*(1-data.0[,i+half])  # all <ND = 0 in data.1
  data.1[,i][data.1[,i] < DLmax] = 0   # all detects below DLmax = 0
  }
  data.2 <- data.1[,1:half]
  r.out <- data.2
  for (i in 1:half) {ranks <- rank(data.2[,i])
  if (i==1) {r.out <- ranks}
  else {r.out <- cbind(r.out, ranks)}
  }
  colnames(r.out) <- nvec
  return(r.out)
}

#' Partial plots for censored MLE regression
#'
#' @description Draws a partial plot for each X variable in regression of a censored Y variable against multiple X variables.
#' @param y.var The column of y (response variable) values plus detection limits.
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value is in `y.var`.
#' @param x.vars Multiple uncensored explanatory variable(s). See Details
#' @param LOG Indicator of whether to compute the censored regression in the original y units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original Y units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @param smooth.method Method for drawing a smooth on the partial plot.  Options are c("gam", "none"). "gam" is a censored generalized additive model using the cenGAM and mgcv packages.
#' @param gam.method	Method for computing the gam smooth.  See the mgcv package for options.  Default is a thinplate ("tp") spline.  "cs" is another good option.
#' @param multiplot  If TRUE, plots are drawn 6 per page.  If FALSE, all plots are drawn on a separate page.
#' @export
#' @return
#' When `x.vars` is one variable, a message is printed that partial plots are not possible with only one X variable and execution stops.
#' When `x.vars` is a data frame of more than one variable, partial plots are drawn for each X variable and text is printed comparing AICs for regression using the untransformed X variable with log and cube-root transforms of the X variable, as a supplement to evaluating linearity on the partial plots alone.
#'
#'
#' @importFrom survival survreg Surv
#' @importFrom cenGAM mgcv
#'
#' @details
#' Partial plots for uncensored data often are drawn with superimposed smooths. At times looking only at the data values without a smooth can better enable the human eye to determine whether the overall pattern is linear or not.  If this is the best method for you, use the smooth.method = "none" option to not draw a smooth.  The most common smooth used for uncensored data is loess, which does not recognize censored data and so uses the detection limit (DL) value itself.  This results in biased-high smooths that incorrectly treat values at the DLs equal to uncensored (detected) data. The partplots function in NADA2 was written to provide a better alternative, smoothing the partial residual pattern with a censored generalized additive model (gam).  The censored gam recognizes the nondetects as left-censored data with a maximum at the DL when computing the smooth. DLs may vary with each observation -- multiple DLs in a dataset are not a problem in routines of the NADA2 package.
#'
#' 'y.var': The default is that the Y variable will be log transformed.
#'
#' `x.vars`: Enter the name of a data frame of columns of the x variables. No extra columns unused in the regression allowed. Create this by `x.frame <- data.frame (Temp, Flow, Time)` for 3 variables (temperature, flow and time).
#'
#' Gray open circles represent censored data and are the residual between the detection limit and the predicted value from the censored regression.  As predictions recognize that the detection limit is an upper limit, predicted values on the regression line are most often below the detection limit, leading to positive residuals.  Note that the true residual could be anywhere below that value.  That fact is recognized by the censored regression but is difficult to represent on a plot.
#'
#' AIC for regression models with untransformed X, log and cube-root transforms of X are printed to evaluate which of the three transformationss results in the ‘best’ model.
#'

#' @seealso [survival::survreg]
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Cook, R.D., 1993. Exploring Partial Residual Plots, Technometrics 35, 351-362.
#'
#' @examples
#'
#' data(Brumbaugh)
#'
#' # For demostration purposes
#' partplots (Brumbaugh$Hg,Brumbaugh$HgCen,Brumbaugh[,c("SedMeHg","PctWetland")])

partplots <- function (y.var, cen.var, x.vars, LOG = TRUE, smooth.method = "gam", gam.method = "tp", multiplot = TRUE)
{  yname <- deparse(substitute(y.var))
nonas <- na.omit(cbind(y.var, cen.var, x.vars))
xnona <- data.frame(nonas[,-(1:2)])
pch.symb <- (nonas[, 2] + 16) - nonas[, 2] *16
col.symb <- (nonas[, 2] + 1) + nonas[, 2] *6
oldpar<- par(no.readonly=T)
on.exit(par(oldpar))
if (multiplot == TRUE) {par(mfrow = c(3,2))}

# function to put data into format used by the gam function tobit1 model
dat4cengam <- function(y.orig, cen.orig) {
  y.gam <- y.orig*(1-cen.orig)     # low end of possible values. <DL data = 0.
  y.gam [y.gam == 0] <- -Inf       # <DL data = -Inf
  dl.gam <- y.orig*cen.orig         # 0 for detects.  DL for nondetects.
  dat.gam <- invisible(data.frame(y.gam, dl.gam))
  return(dat.gam)
}

if (LOG == TRUE)  {
  if (min(nonas[,1]) <= 0) stop('Cannot take logs of zero or negative values.  Use LOG=FALSE')
  lnvar <- log(nonas[,1])    # Y in log units (default)
  yname.log = paste("ln(", yname, ")", sep = "")
  flip.log <- max(lnvar) - lnvar
  surv.log <- Surv(flip.log, as.logical(1-nonas[,2]) )

  if (is.data.frame(x.vars))  {             # multiple x variables
    yname.part <- paste(yname.log, ".partial", sep="")
    for (i in 1:ncol(xnona))  {temp.x <- data.frame(xnona[, -i])
    xname.part <- paste(colnames(xnona[i]), ".partial", sep="")
    y.part <- survreg(surv.log ~ ., data = temp.x, dist = "gaussian")
    x.part <- lm(xnona[, i] ~ ., data = temp.x)
    flip.resi <- flip.log - y.part$linear.predictors
    # y in original, not flipped, units
    ylog.resi <- -1*flip.resi
    ylog.pred <- max(lnvar) - y.part$linear.predictors
    # y coefficients in original units
    y.coefficients <- y.part$coefficients * (-1)
    y.coefficients[1] <- max(lnvar) + y.part$coefficients[1]

    #  plot y resids vs x resids with smooth

    # option 1:  no smooth
    plot(x.part$residuals, ylog.resi, xlab = xname.part, ylab = yname.part, pch = pch.symb, col = col.symb)

    # option 2:  censored GAM
    if (smooth.method == "gam")
    { data.gam <- dat4cengam(ylog.resi, nonas[,2])
    ylog.gam <- gam(data.gam[,1] ~s(x.part$residuals), bs = gam.method, family = tobit1(link = "identity", left.threshold = data.gam[,2]))
    o <- order(x.part$residuals , ylog.resi)
    lines (x.part$residuals[o], ylog.gam$fitted.values[o], col = 'red')}

    if (multiplot == FALSE) {title(paste("Open circles are", yname, "nondetects"))}

    # compute and print r2 and AIC for original and possible x transforms
    cat(colnames(xnona[i]), "\n", "untransformed", "\n")
    notrans <- cencorreg(nonas[,1], nonas[,2], xnona, verbose = FALSE)
    AIC.none <- -2*notrans$loglik[2] + (2*notrans$df +1)

    cat("cube root", "\n")
    x.sign <- sign(xnona[i])
    x.cube <- abs(xnona[i])**(1/3) * x.sign
    cube.xnona <- cbind(x.cube, temp.x)
    cube.trans <- cencorreg(nonas[,1], nonas[,2], cube.xnona, verbose = FALSE)
    AIC.cube <- -2*cube.trans$loglik[2] + (2*cube.trans$df +1)

    cat("log transform", "\n")
    if (min(xnona[i]) > 0) {
      x.log <- log(xnona[i])
      log.xnona <- cbind(x.log, temp.x)
      log.trans <- cencorreg(nonas[,1], nonas[,2], log.xnona, verbose = FALSE)
      AIC.log <- -2*log.trans$loglik[2] + (2*log.trans$df +1) }
    else {AIC.log = AIC.cube
    cat("Cannot take logs of zero or negative values.", "\n")}
    AIC.diff <- AIC.none - min(AIC.log, AIC.cube)
    cat("Decrease in AIC from transformation of", colnames(xnona[i]), "=", max(0,AIC.diff), "\n", "\n")
    }  # end of cycle thru x variables

  }  # end of multiple x vars
  else (stop("For only one x variable partial plots not needed."))
}    # end of logs.

#  In original Y units
else {
  flip.y <- max(nonas[,1]) - nonas[,1]
  surv.y <- Surv(flip.y, as.logical(1-nonas[,2]) )

  if (is.data.frame(x.vars))  {
    yname.part <- paste(yname, ".partial", sep="")
    # multiple x variables
    for (i in 1:ncol(xnona))  {temp.x <- data.frame(xnona[, -i])
    xname.part <- paste(colnames(xnona[i]), ".partial", sep="")
    y.part <- survreg(surv.y ~ ., data = temp.x, dist = "gaussian")
    x.part <- lm(xnona[, i] ~ ., data = temp.x)
    flip.resi <- flip.y - y.part$linear.predictors
    # y in original, not flipped, units
    y.resi <- -1*flip.resi
    y.pred <- max(nonas[,1]) - y.part$linear.predictors
    # y coefficients in original units
    y.coefficients <- y.part$coefficients * (-1)
    y.coefficients[1] <- max(nonas[,1]) + y.part$coefficients[1]

    # option 1: no smooth
    plot(x.part$residuals, y.resi, xlab = xname.part, ylab = yname.part, pch = pch.symb, col = col.symb)

    # option 2:  censored GAM
    if (smooth.method == "gam")
    { data.gam <- dat4cengam(y.resi, nonas[,2])
    y.gam <- gam(data.gam[,1] ~s(x.part$residuals), bs = gam.method, family = tobit1(link = "identity", left.threshold = data.gam[,2]))

    o <- order(x.part$residuals , y.resi)
    lines (x.part$residuals[o], y.gam$fitted.values[o], col = 'red')}

    if (multiplot == FALSE) {title(paste("Open circles are", yname, "nondetects"))}

    # compute and print r2 and AIC for original and possible x transforms
    cat(colnames(xnona[i]), "\n", "untransformed", "\n")
    notrans <- cencorreg(nonas[,1], nonas[,2], xnona, verbose = FALSE)
    AIC.none <- -2*notrans$loglik[2] + (2*notrans$df +1)

    cat("cube root", "\n")
    x.sign <- sign(xnona[i])
    x.cube <- abs(xnona[i])**(1/3) * x.sign
    cube.xnona <- cbind(x.cube, temp.x)
    cube.trans <- cencorreg(nonas[,1], nonas[,2], cube.xnona, verbose = FALSE)
    AIC.cube <- -2*cube.trans$loglik[2] + (2*cube.trans$df +1)

    cat("log transform", "\n")
    if (min(xnona[i]) > 0) {
      x.log <- log(xnona[i])
      log.xnona <- cbind(x.log, temp.x)
      log.trans <- cencorreg(nonas[,1], nonas[,2], log.xnona, verbose = FALSE)
      AIC.log <- -2*log.trans$loglik[2] + (2*log.trans$df +1) }
    else {AIC.log = AIC.cube
    cat("Cannot take logs of zero or negative values.", "\n")}
    AIC.diff <- AIC.none - min(AIC.log, AIC.cube)
    cat("Decrease in AIC from transformation of", colnames(xnona[i]), "=", max(0,AIC.diff), "\n", "\n")

    }  # end of cycle thru multiple x variables
  }  # end of multiple variables
  else {stop("For only one x variable partial plots not needed.")}
} # end of original units
par(oldpar)
}

#' PbHeron
#'
#' @description
#' From the `NADA` R-Package.
#'
#' Lead concentrations in the blood and several organs of herons in Virginia.
#'
#' Objective is to determine the relationships between lead concentrations in the blood and various organs. Do concentrations reflect environmental lead concentrations, as represented by dosing groups? There is one detection limit, at 0.02 ug/g. Used in Chapters 10 and 11 of the NADA book.
#'
#' @usage data(PbHeron)
#'
#' @docType data
#' @keywords dataset
#' @name PbHeron
#' @source Golden et al., 2003, Environmental Toxicology and Chemistry 22, pp. 1517-1524.
#' @references
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.

"PbHeron"

#' Censored Empirical Cumulative Distribution Function
#'
#' @description Plots cdfs of one or more groups of censored data.  Illustrates the differences between groups for group tests such as those done using cen1way or cenanova.
#' @param y.var The column of data values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param group Optional - grouping or factor variable. Can be either a text or numeric value indicating the group assignment.
#' @param logscale Indicator of whether to plot the data in original y units, or on a log scale.  The default is to use the original units, which usually makes more sense for a cdf plot.
#' @param Ylab Optional – input text in quotes to be used as the variable name on the cdf plot.  The default is the name of the `y.var` input variable.
#' @keywords CDF
#' @return Plots an empirical cumulative distribution function. If `group=NULL` then a ECDF with 95% confidence interval is produced. If `group` is identified then ECDFs are produced for each group.
#' @export
#' @importFrom NADA cenfit flip Cen
#' @importFrom graphics axis box
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#'@examples
#' data(PbHeron)
#'
#' # with groups
#' with(PbHeron,plotcdf(Liver,LiverCen,DosageGroup))
#'
#' # all data
#' with(PbHeron,plotcdf(Liver,LiverCen))

plotcdf <- function(y.var, cen.var, group=NULL, logscale=FALSE, Ylab=varname) {
  varname <- deparse(substitute(y.var))
  log = ""
  oldpar<- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (logscale == TRUE) {log="x"}
  if ( is.null(group) ) {

    cenind <- as.logical(cen.var)
    cenfit.out <- cenfit(y.var, cenind)

    plot(cenfit.out@survfit, log=log, xlab=Ylab, lwd = 2, axes=F, main = "Cumulative Distribution Function")
    abline(h=1, lwd=3, col = "white")
    box(bty = "o")
    axis(side = 1)
    axis(side = 2)
  }
  else {
    Factor <- as.factor(group)
    cenind <- as.logical(cen.var)
    cenfit.out <- cenfit(y.var, cenind, Factor)
    grpname <- deparse(substitute(group))
    i <- length(levels(Factor))
    clrs <- c (1:i)

    if (logscale == TRUE)  {
        plot(cenfit.out@survfit, log=log, xlab=Ylab, lwd = 2, col=clrs, axes=F, main = "Cumulative Distribution Function")
        abline(h=1, lwd=3, col = "white")
        box(bty = "o")
        axis(side = 1)
        axis(side = 2)
    }
    else {    # logscale is false
        plot(cenfit.out@survfit, log=log, xlab=Ylab, lwd = 2, col=clrs, axes=F, main = "Cumulative Distribution Function")
        abline(h=1, lwd=3, col = "white")
        box(bty = "o")
        axis(side = 1)
        axis(side = 2)
    }
    legend("bottomright",levels(Factor), lty=1:i, lwd=2, text.col=clrs, col = clrs, title = grpname)
  }
  }

#' Test for difference in left-censored samples
#' @description Performs a nonparametric Paired Prentice-Wilcoxon test of whether the median difference between two columns of paired censored data equals 0 (O'Brien and Fleming, 1987)
#' @param xd The first column of data values plus detection limits
#' @param xc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the xd column, and 0 (or `FALSE`) indicates a detected value in `xd`.
#' @param yd The second column of data values plus detection limits
#' @param yc The column of censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the yd column, and 0 (or `FALSE`) indicates a detected value in `yd`
#' @param alternative The usual notation for the alternate hypothesis.  Default is `“two.sided”`.  Options are `“greater”` or `“less”`.
#'
#' @importFrom survival survfit Surv
#' @importFrom stats na.exclude pnorm
#' @return Paired Prentice-Wilcoxon test results including Z-statistic, n (sample size), p-value and median difference
#' @export
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' O’Brien, P.C., Fleming, T.R., 1987. A Paired Prentice-Wilcoxon Test for Censored Paired Data. Biometrics 43, 169–180. https://doi.org/10.2307/2531957
#'
#' @seealso [survival::survfit] [survival::Surv]
#'
#' @examples
#' data(PbHeron)
#' ppw.test(PbHeron$Liver,PbHeron$LiverCen,PbHeron$Bone,PbHeron$BoneCen)
#'

ppw.test <- function(xd, xc, yd, yc, alternative="two.sided")
{
  xname <- deparse(substitute(xd))
  yname <- deparse(substitute(yd))
  OBrienFleming=TRUE
  ## Error checks
  if(length(xd) != length(yd))
    stop("Lengths of x and y must be the same for paired data.")
  keep <- !(is.na(xd) | is.na(yd) | is.na(xc) | is.na(yc))
  x <- xd[keep]
  y <- yd[keep]
  if(any(c(x, y) < 0))
    stop("Negative values in x or y")
  N <- length(x)
  if(OBrienFleming) {
    xd <- x
    xc <- xc[keep]
    yd <- y
    yc <- yc[keep]

    for(i in seq(N)) {
      if(xc[i] && yc[i]) { # Both censored, make same
        xyMax <- max(xd[i], yd[i])
        xd[i] <- yd[i] <- xyMax
      }
      else if(xc[i] && xd[i] > yd[i]) { # x censored > y
        yd[i] <- xd[i]
        yc[i] <- TRUE
      }
      else if(yc[i] && yd[i] > xd[i]) { # y censored > x
        xd[i] <- yd[i]
        xc[i] <- TRUE
      }
    } # Done, no need to check uncensored observations
  }
  ## Test requires stacked data
  group = factor(c(rep(xname, length(x)), rep(yname, length(y))),
                 levels=c(xname, yname))
  values <- c(xd, yd)
  cenflag <- as.logical(c(xc, yc))
  alternative <- pmatch(alternative, c("two.sided", "greater", "less"))
  if(is.na(alternative))
    stop('Invalid choice for alternative, must match "two.sided", "greater", or "less"')
  ## Define the test:
  PPW.test <- function(values, cenflag, group) {
    ## data must be in exactly 2 groups of equal size and stacked
    ## x first, then y (for alternative not two.sided)
    ## required adjustment to values
    adjust <- min(diff(sort(unique(values))))/10.
    values <- ifelse(cenflag, values-adjust, values)
    adjust <- adjust / length(values)
    dupes <- unique(values[!cenflag][duplicated(values[!cenflag])])
    in.values <- values
    if(length(dupes) > 0) { #there are dupes
      for(i in seq(along=dupes)) {
        sel <- which(values == dupes[i])
        mult <- 0L:(length(sel) - 1L)
        in.values[sel] <- dupes[i] + adjust * mult
      }
    }
    ## create data frame and add observed value at 0 to compute
    ## correct probs
    df <- data.frame(values=c(0, -in.values), cenflag=c(T, !cenflag))
    kmout <- survfit(Surv(values, cenflag) ~ 1, data=df, na.action=na.exclude)
    kmout$time <- - kmout$time # convert back to actual values
    ## compute mean survival for tied values
    St <- kmout$surv
    if(length(dupes) > 0L) { #there are dupes
      for(i in seq(along=dupes)) {
        ## take advantage of the fact that the output are sorted
        sel <- which(values == dupes[i])
        mult <- seq(0L, 1L - length(sel))
        sel <- which(kmout$time == dupes[i])
        St[sel] <- mean(St[sel + mult])
      }
    }

    ## Define the link between the observed data and the kaplan meier table
    ## and compute needed stats.
    link <- match(values, kmout$time)
    St <- St[link]
    Uncen <- values[!cenflag]
    UncenSt <- St[!cenflag]
    UncenSt <- UncenSt[order(Uncen)]
    Uncen <- sort(Uncen)
    Score <- 1. - 2*St
    ## This is not fast, but it works
    for(i in which(cenflag)) # fix each each censored score
      Score[i] <- 1. - UncenSt[which(Uncen > values[i])][1L]
    ## Compute d
    ## Reverse sense of d so that alternatives are in same direction--
    ## gives different sense of d from the original
    Score <- matrix(Score, ncol=2L)
    d <-  Score[, 2L] - Score[, 1L]
    return(list(Z=sum(d)/sqrt(sum(d^2)), Scores=Score, Diffs=d))
  } # end of PPW.test

  ret1 <- PPW.test(values, cenflag, group)
  stat <- ret1$Z
  names(stat) <- "Paired Prentice Z"
  meth <- "Paired Prentice-Wilcoxon test"
  param <- length(group)/2L
  names(param) <- "n"
  ## add finishing touches
  if(alternative == 1L)   # alternative is two-sided
  {  pvalue <- (1. - pnorm(abs(stat))) * 2.
     altern <- paste(xname, "not equal to", yname, sep = " ") }
  else if(alternative == 2L)  # alternative is greater than
  {  pvalue <- 1. - pnorm(stat)
     altern <- paste(xname, ">", yname, sep = " ") }
  else # alternaitve is less than
  {  pvalue <- pnorm(stat)
     altern <- paste(xname, "<", yname, sep = " ") }
  names(pvalue) <- "p value"
  mu <- 0
  names(mu) <- "difference"
  Scoremat <- cbind(ret1$Scores, ret1$Diffs)
  colnames(Scoremat) <- c("xScore", "yScore", "d")
  ## For diagnostic plot, create min and max differences
  d1 <- xd - ifelse(yc, 0, yd)
  d2 <- ifelse(xc, 0, xd) - yd
  mind <- pmin(d1, d2)
  maxd <- pmax(d1, d2)
  Scoremat <- cbind(Scoremat, minDiff=mind, maxDiff=maxd)

  sv.out <- survfit(Surv(Scoremat[,4], Scoremat[,5], type="interval2")~1)
  med.diff <- min(sv.out$time [sv.out$surv <= 0.50])

  # getting number of signif digits
  y.count <- vector(length=length(yd))
  for (i in 7:1) {y.count[yd == signif(yd, i)] <- i}
  # computing median difference
  median.diff <- paste ("Median difference equals", signif(med.diff, max(y.count)))

  retval <- list(statistic = stat, parameters = param,
                 p.value = pvalue, null.value = mu,
                 alternative = c("two.sided", "greater", "less")[alternative],
                 method = meth, data.name = paste(xname, "and", yname, sep = " "),
                 PPWmat=Scoremat)
  #oldClass(retval) <- c("htest", "ppw")

#  print(retval)
  txt <- paste("Paired Prentice Wilcoxon test for (x:", xname, " - ", "y:", yname, ") equals 0", "\n", "    alternative: ", altern, "\n", sep = "")
  txt2 <- paste("n =", param, "  Z =", signif(stat, 4), "  p-value =", signif(pvalue, 4))

  cat(txt, "\n", txt2, "\n")
  cat(" ", median.diff, "\n")
  return(invisible(retval))
}

#' Computes confidence intervals on regression on order statistics (ROS) mean
#'
#' @description Uses ROS model output from the NADA package and computes the Zhou and Gao 1997 modified Cox’s method two-sided confidence interval around the mean for a lognormal distribution.  Computes a t-interval for a gaussian ROS model output.
#' @param cenros.out an ROS model output object (see details)
#' @param conf Confidence coefficient of the interval (Default is 0.95)
#' @return Prints a lower (LCL) and upper (UCL) confidence interval based on the `conf` provided (Default is 95%)
#' @details
#' This function uses an ROS model output based on the `ros` funtion in the `NADA` package.  The lognormal distribution is the default for the NADA package but a gaussian distribution is optional.
#' For more detail on ROS modeling see the `ros` help file (`?NADA::ros`).
#'
#' For implementation of `ROSci(...)` see the examples below.
#' @importFrom stats qt sd
#' @export
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Lee, L., Helsel, D., 2005. Statistical analysis of water-quality data containing multiple detection limits: S-language software for regression on order statistics. Computers & Geosciences 31, 1241–1248. <https://doi.org/10.1016/j.cageo.2005.03.012>
#'
#' Zhou, X.-H., Gao, S., 1997. Confidence Intervals for the Log-Normal Mean. Statistics in Medicine 16, 783–790. <https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7<783::AID-SIM488>3.0.CO;2-2>
#'
#' @seealso [NADA::ros]
#'
#' @examples
#' data(Brumbaugh)
#' myros <- NADA::ros(Brumbaugh$Hg,Brumbaugh$HgCen)
#'
#' summary(myros)
#'
#' # ROS Mean
#' mean(myros)
#'
#' # 95% CI around the ROS mean
#' ROSci(myros)


ROSci <- function(cenros.out, conf=0.95) {
  p <- 1-((1-conf)/2)
  n <- length(cenros.out$obs)
  c <- as.character(100*conf)
  ciname1 <- paste ("   LCL", c, sep = "")
  ciname2 <- paste ("    UCL", c, sep = "")
  cinames <- paste (ciname1, ciname2)
  if (cenros.out$forwardT == "log") {
    # For the default lognormal distribution, Cox’s method is used.
    scale <- as.vector(cenros.out[1]$coefficients[2])
    #  coefficient #2 is the scale, the sd of the logs
    bhat <- log(mean(cenros.out))
    #  Mean log plus one-half scale^2 ;  the mean transformed to log units
    gamz <- qt(p,(n-1)) * sqrt((scale^2/n) + (((0.5)*scale^4)/(n-1)))
    #  Zhou and Gao 1997 modified Cox’s method for the CI of a lognormal distribution
    cilog <- c(exp(bhat - gamz), exp(bhat + gamz))
    cat(cinames, "    assuming a lognormal distribution", "\n")
    cat(cilog, sep = "  ")}

  else { # for a gaussian distribution;
    tstat <- qt(p,(n-1))
    halfw <- tstat*(sd(cenros.out)/sqrt(n))
    ciNorm <- c(mean(cenros.out) - halfw, mean(cenros.out) + halfw)
    cat(cinames, "    assuming a gaussian distribution", "\n")
    cat(ciNorm, sep = "  ")
  }

}


#' Plot U-score Nonmetric Multidimensional Scaling of censored data
#'
#' @description Plots an NMDS of uscores output from the `uscores` or `uscoresi` functions.
#' @param uscor A data frame of uscores or ranks of uscores produced by either the `uscores(...)` or `uscoresi(...)` functions
#' @param group Optional grouping variable. Sites will be represented by different colored symbols for each group.
#' @param title Optional title for the NMDS graph.
#' @param legend.pos For when group is specified, the location of the legend on the graph showing the colors representing each group’s data.  Default is “bottomleft”.  Alternatives are “topright” and “centerleft”, etc.
#' @importFrom vegan metaMDS
#' @importFrom stats dist
#' @return Prints an NMDS plot of censored data groupings based on U-scores
#' #' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @export
#'
#' @examples
#' data(PbHeron)
#'
#' PbHeron.u <- uscores(PbHeron[,4:15])
#' uMDS(PbHeron.u)
#'
#' # With group specific
#' uMDS(PbHeron.u,group=PbHeron$DosageGroup)

uMDS <- function(uscor, group = NULL, title=NULL, legend.pos = "bottomleft") {
  uname <- deparse(substitute(uscor))
  if (is.null(title)) {title <- paste("NMDS of", uname)}
  rownumb <- nrow(uscor)
  row.labels <- as.character(1:rownumb)
  uclid <- dist(uscor)  # euclidean distance matrix of uscores.
  u.mds <- metaMDS(uclid, zerodist = "add")
  if (is.null(group)) {plot(u.mds, type = "p", main = title)
    text(u.mds$points, pos=4, cex = 0.9, offset = 0.25)
  }
  else {
    gp.factor <- as.factor(group)
    gp.nlevels <- c(1:nlevels(gp.factor))
    gp.col <- as.integer(gp.factor)
    gp.symb <- 18+gp.col
    plot.col=ordiplot(u.mds, type="none", display="sites", pch = 19, main=title)
    points(plot.col, "sites", pch=19, col=gp.col)
    text(plot.col, "sites", pos=4, cex = 0.8)
    legend(legend.pos, legend=levels(gp.factor), bty="n",col = gp.nlevels, pch = 19)
  }
}

#' U-scores for (non-interval, sinle-column) Censored Data
#'
#' @description Computes the column of uscores from 2 columns of data in the indicator value format. Multiple detection limits allowed.  Called by the uscores function, uscore (this function) is not expected to be of much use to users on its own.
#' @param y The column of data values plus detection limits
#' @param ind The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y` column, and 0 (or `FALSE`) indicates a detected value in `y`.
#' @param rnk  A `TRUE`/`FALSE` variable on whether to compute the multivariate pattern on the uscores, or the ranks of the uscores.  Default is rnk=`TRUE`, use the ranks. rnk = `FALSE` returns the uscores.
#' @export
#' @return Returns a single column of uscores or the ranks of uscores for a single pair of (concentration, indicator) censored data columns.
#'
#' @examples
#' data(Brumbaugh)
#' uscore(Brumbaugh$Hg,Brumbaugh$HgCen)

Usc <- function(y, ind, rnk=TRUE){
  x <- na.omit(data.frame (y, ind))
  n=length(x$y)
  ylo=(1-as.integer(x$ind))*x$y
  # yadj=y
  yadj=x$y-(sign(x$y-ylo)*0.001*x$y)
  overlap=x$y
  Score=overlap
  for (j in 1:n) {
    for (i in 1:n ){
      overlap[i]=sign(sign(yadj[i]-ylo[j])+sign(ylo[i]-yadj[j]))
    }
    Score[j] = -1*sum(overlap) # -1 so that low values = low scores
  }

  if (rnk) {uscore=rank(Score)} else {uscore = Score}
  return(uscore)
}

#' Interval-censored U-Score
#'
#' @description Interval-censored computation of uscores and their ranks for 1 parameter.  Called by uscoresi. Usci is not expected to be of much use to users on its own.
#' @param ylo The lower end of the concentration interval
#' @param yhi The upper end of the concentration interval
#' @param rnk A `TRUE`/`FALSE` variable on whether to compute the multivariate pattern on the uscores, or the ranks of the uscores.  Default is rnk=`TRUE`, use the ranks. rnk = `FALSE` returns the uscores.
#' @return Returns a single column of uscores or the ranks of uscores for a single pair of (low, high) interval-censored data columns.
#'
#' @export


Usci <- function(ylo, yhi, rnk=T){
  x <- na.omit(data.frame (ylo, yhi))
  n = length(x$ylo)
  yadj=x$yhi-(sign(x$yhi-x$ylo)*0.001*x$yhi) #sets a <1 to be <1
  overlap=x$yhi  # sets correct dimensions
  Score=overlap  # sets correct dimensions
  for (j in 1:n) {
    for (i in 1:n ){
      overlap[i]=sign(sign(yadj[i]-x$ylo[j])+sign(x$ylo[i]-yadj[j]))
    }
    Score[j] = -1*sum(overlap) # -1 so that low values = low scores
  }

  if (rnk) {Uscore=rank(Score)} else {Uscore = Score}
  # print(Score)
  # print(Uscore)
  return(Uscore)
}

#' Uscores for multiple columns of censored data
#'
#' @description Computes uscores of censored data in the indicator format. Multiple DLs allowed.
#' @param dat.frame A data frame. Default format is paired = `TRUE`, where for 3 chemical parameters the input format is C1 I1 C2 I2 C3 I3, a concentration column followed by its corresponding indicator column.
#' @param paired When paired = `FALSE`, the input format is C1 C2 C3 I1 I2 I3 where the C columns contain concentrations or a detection limit, and the I columns are their associated indicators, in the same order as the concentration columns.
#' @param rnk A logical `TRUE`/`FALSE` variable on whether to compute the multivariate pattern on the uscores, or the ranks of the uscores.  Default is rnk=`TRUE`, use the ranks. rnk = `FALSE` returns the uscores.
#'
#' @return A matrix of uscores, one column for each chemical parameter.  If there is only one chemical parameter a vector of uscores is returned.
#' @export
#'
#' @examples
#' data(PbHeron)
#'
#' uscores(PbHeron[,4:15])

uscores <- function(dat.frame, paired = TRUE, rnk=TRUE) {
  NAtext <- "Rows were deleted due to NAs. Will need to delete same rows in grouping variable before running NMDS or anosim."
  cols <- ncol(dat.frame)
  if (cols >2) { half <- cols/2     # multiple pairs of columns
  dat.nonas <- as.matrix(na.omit (dat.frame))
  j=0;  uscr=0; DLmax=0; data.1 <- dat.nonas; data.0 <- dat.nonas
  if (paired) { for (i in seq(1, to=(cols-1), by=2)) {j = j+1
  data.0[,j] <- dat.nonas[,i]
  data.0[,j+half] <- dat.nonas[,i+1]   # reorganizing paired columns
  nvec <- colnames(data.0)
  nvec <- nvec[seq(1, (cols-1), 2)]
  nvec <- paste("usc.", nvec, sep="")
  } }
  else {data.0 <- dat.nonas   # no reorg needed.  columns in blocks.
  nvec <- colnames(data.0)
  nvec <- nvec[1:half]
  nvec <- paste("usc.", nvec, sep="")
  }

  u.out <- data.0[,1:half]
  for (i in 1:half) {u <- Usc(data.0[,i], data.0[,i+half], rnk=rnk)
  if (i==1) {u.out <- u}
  else {u.out <- cbind(u.out, u)}
  }
  # check for deletion of NAs
  if (length(dat.nonas[,1]) < length(dat.frame[,1])) {warning(NAtext)}

  colnames(u.out) <- nvec
  return(u.out)
  }

  else {    # only 1 pair of columns
    x <- na.omit(dat.frame)
    names(x) <- c("y", "ind")
    n=length(x$y)
    ylo=(1-as.integer(x$ind))*x$y
    yadj=x$y-(sign(x$y-ylo)*0.001*x$y)
    overlap=x$y  # sets correct dimensions
    Score=overlap  # sets correct dimensions
    for (j in 1:n) {
      for (i in 1:n ){
        overlap[i]=sign(sign(yadj[i]-ylo[j])+sign(ylo[i]-yadj[j]))
      }
      Score[j] = -1*sum(overlap)     # -1 so that low values = low scores
    }

    if (rnk) {uscore=rank(Score)} else {uscore = Score}
    # check for deletion of NAs
    if (length(x[,1]) < length(dat.frame[,1])) {warning(NAtext)}

    return(uscore)
  }
}


#' U-scores for interval-censored data (multiple columns)
#'
#' @description Computes uscores within columns of interval-censored data (the "i").  Data may have one or more detection limits.
#' @param dat.frame A data frame. Default format is: paired = `TRUE`, the default input format is ylo1 yhi1 ylo2 yhi2 ylo3 yhi3, where ylo is the low end of the interval of possible values, and yhi is the high end. There is a pair of columns for each censored parameter.
#' @param paired An option to specify paired = `FALSE`, where the format would be ylo1 ylo2 ylo3 yhi1 yhi2 yhi3, low values for each parameter followed by the high values in the same order.
#' @param rnk rnk=`TRUE` returns the ranks of uscores.  rnk = `FALSE` returns the uscores themselves. Default is rnk = `TRUE` to return the ranks.
#' @param Cnames Cnames =`1` uses the "lo" column names to name the uscores columns (the default).  Cnames = `2` uses the "hi" column names.
#' @details Input is a data.frame of paired low and high possible range of values, in an interval-censored format. ylo = the lower end of the interval is the first (left) column in the pair.  yhi is the upper end of the interval, at the second (right) column in the pair. For a detected value, ylo=yhi.  For a ND,  ylo != yhi. The uscore is the number of observations known to be lower minus the number of observations known to be higher. Ties, including those such as <1 vs <3 or 4 vs 4 or <3 vs 2, are given a 0 value in the uscore computation.
#'
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @export
#'


uscoresi <- function(dat.frame, paired = TRUE, rnk=TRUE, Cnames = 1) {
  NAtext <- "Rows were deleted due to NAs. Will need to delete same rows in grouping variable before running NMDS or anosim."
  cols <- ncol(dat.frame)
  # for multiple pairs of columns
  if (cols > 2) {half <- cols/2
  dat.nonas <- as.matrix(na.omit (dat.frame))
  j=0; data.0 <- dat.nonas
  if (paired) { for (i in seq(1, to=(cols-1), by=2)) {j = j+1
  data.0[,j] <- dat.nonas[,i]
  data.0[,j+half] <- dat.nonas[,i+1]   # reorganizing paired columns.
  nvec <- colnames(data.0)
  if (Cnames == 1) {nvec <- nvec[seq(1, (cols-1), 2)] }
  else {nvec <- nvec[seq(2, cols, 2)] }
  nvec <- paste("usc.", nvec, sep="")
  } }
  # no reorg needed,  Not paired, columns in blocks.
  else { data.0 <- dat.nonas
  nvec <- colnames(data.0)
  nvec <- nvec[1:half]
  nvec <- paste("usc.", nvec, sep="")
  }

  u.out <- data.0[,1:half]
  for (i in 1:half) {u <- Usci(data.0[,i], data.0[,i+half], rnk=rnk)
  if (i==1) {u.out <- u}
  else {u.out <- cbind(u.out, u)}
  }

  # check for deletion of NAs
  if (length(dat.nonas[,1]) < length(dat.frame[,1])) {warning(NAtext)}

  colnames(u.out) <- nvec
  return(u.out)
  }
  else {                   # only 1 pair of columns
    x <- na.omit(data.frame (dat.frame))
    names(x) <- c("ylo", "yhi")
    n = length(x$ylo)
    yadj=x$yhi-(sign(x$yhi-x$ylo)*0.001*x$yhi)   #  sets a <1 to be <1
    overlap=x$yhi  # sets correct dimensions
    Score=overlap  # sets correct dimensions
    for (j in 1:n) {
      for (i in 1:n ){
        overlap[i]=sign(sign(yadj[i]-x$ylo[j])+sign(x$ylo[i]-yadj[j]))
      }
      Score[j] = -1*sum(overlap) # -1 so that low values = low scores
    }

    if (rnk) {Uscore=rank(Score)} else {Uscore = Score}
    # check for deletion of NAs
    if (length(x[,1]) < length(dat.frame[,1])) {warning(NAtext)}

    return(Uscore)
  }
}


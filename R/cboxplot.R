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

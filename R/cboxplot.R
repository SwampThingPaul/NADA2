#' Draws censored boxplots
#'
#' @description
#' Draws boxplots for left-censored data with one ore more detection limit(s). Portions below the maximum detection limit(s) are not shown by default, as their percentiles are not known.
#' @param x1 The column of x (response variable) values plus detection limits.
#' @param x2 The x-variable censoring indicators, where 1 (or `TRUE`) indicates a detection limit in the y1 column, and 0 (or `FALSE`) indicates a detected value in y1.
#' @param xgroup An optional column of a grouping variable.  Draws side-by-side boxplots if this variable is present.
#' @param LOG `TRUE`/`FALSE` indicator of whether to plot the Y axis data on the original scale (`FALSE`) or log scale (`TRUE`).
#' @param show `TRUE`\/`FALSE` indicator of whether to show estimated values computed using ROS for the portion of the box below the maximum DL (`TRUE`), or just leave the lower portion blank (`FALSE`).
#' @param minmax `TRUE`/`FALSE` indicator of whether to draw outliers individually. Default is to show outliers. Setting `minmax = TRUE` will draw the whiskers out to the max and min of the dataset.
#' @param ordr A vector indicating the order of boxes to be drawn on the boxplot, if not in alphabetical order (the default).  Example: for 4 boxplots for groups A, B, C, D, to change the order to the reverse type ordr = c(4, 3, 2, 1).  Example 2: To change the order to A, C, D, B, type ordr = c(1, 3, 4, 2)
#' @param Ylab Y axis label text, if something is wanted other than the Y variable name in the dataset.
#' @param Xlab X axis label text, if something is wanted other than the group variable name in the dataset.
#' @param Title Text to show as the graph title.  Default is blank.
#' @param dl.loc Location indicator of where to plot the "MaxDL=" text on some versions of the plot.  Possible entries are “topleft”, “topright”, “topcenter”, and the corresponding “bottom” text.
#' @param dl.col Color of the max detection limit line(s), and the legend text stating the max DL.  Default is “red”, but all standard R colors may be used.
#' @param dl.lty line type of max detection limit line(s).
#' @param dl.lwd line wide of max detection limit line(s).
#' @param bxcol Color for interior of boxplots. Specify just one color if all boxes are to be the same color.  If a different color is desired for each of three boxplots, as one example, use bxcol = c(“red”, “white”, “blue”) etc.
#' @param Ymax Maximum Y value to be shown on the plot.  Used to cut off high outliers on plot and better show the bulk of the boxplots.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @param Hlines Data to add horizontal reference lines to the boxplot. Required input is a data frame of 4 columns.  See Details.
#' @param ... Additional arguments passed to the default boxplot function

#' @details If maximum detection limits vary among groups, separate maxDL lines will be drawn for each group's boxplot. If one group has fewer than 3 detected observations its boxplot will not be drawn.  Its detection limits will not count when computing the maximum limit.  However, if only one boxplot is drawn for the entire dataset by not specifying a group variable, the detection limits from the portion that is the mostly ND group will be used when computing the maximum limit.
#' @details The reuired input to draw additional horizontal lines (Hlines option) is a data frame with 4 columns of input, one row per horizontal line.  More than one line may be drawn.  Column one is the Y axis value for the line.  Column 2 is the line color, column 3 is the line type (lty) and column 4 is the text to be added just above the line.  To add one line at a value of 40, for example, use Hlines = yline, after defining yline = data.frame(c(40, "purple", "dotted", "New Health Std")).  To draw two lines, define yline as yline = data.frame(matrix(c(40, "purple", "dotted", "New Health Std", 70, "blue", "longdash", "Old Health Std"), ncol = 4, byrow=TRUE))) . If no text is wanted use " " for the column 4 entry for that line. See ?par under lty for standard line types.
#' @export
#' @importFrom graphics boxplot lines plot polygon
#' @importFrom grDevices adjustcolor
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


cboxplot <- function(x1, x2, xgroup=NULL, LOG =FALSE, show=FALSE, ordr = NULL,
                     Ylab=yname, Xlab = gname, Title = NULL,
                     dl.loc = "topright", dl.col = "red",dl.lwd = 2,
                     dl.lty = "longdash",bxcol = "white",
                     Ymax = NULL, minmax = FALSE, printstat = FALSE, Hlines = NULL,...) {

  ## handle inputs
  if (is.null(xgroup)) {
    ydat <- na.omit(data.frame(x1, x2))
    y1 <- ydat[,1]; y2 <- ydat[,2]; group <- NULL
  } else {
    ydat <- na.omit(data.frame(x1, x2, xgroup))
    y1 <- ydat[,1]; y2 <- ydat[,2]; group <- ydat[,3]
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  yname <- deparse(substitute(x1))
  if (is.null(Ylab))  Ylab <- yname
  y1 <- as.numeric(y1)
  y2 <- as.logical(y2)   # TRUE = censored

  ## ----------------------------------
  ## no censored values at all - default to boxplot
  ## ----------------------------------
  if (!any(y2)) {
    warning("No censored data detected. Passing arguments to standard boxplot().")
    Ylab <- deparse(substitute(x1))
    Xlab <- NA
    if (is.null(xgroup)) {
      boxplot(y1, na.action = na.omit, ylab = Ylab, xlab = Xlab,
                        main = Title, col = bxcol, log = ifelse(LOG, "y", ""),...)
      rslt <- data.frame(MaximumDL="none", Group="all")
    } else {
      gname <- deparse(substitute(xgroup))
      if (!is.null(ordr)) group <- factor(group, levels=ordr)
      glabs <- levels(factor(group))
      boxplot(y1 ~ group, na.action = na.omit, ylab = Ylab, xlab = gname,
                        names = glabs, main = Title, col = bxcol, log = ifelse(LOG, "y", ""),...)
      rslt <- data.frame(MaximumDL="none", Group=glabs)
    }
    return(invisible(rslt))
  }

  box.fill <- adjustcolor( "white", alpha.f = 0.6)
  if (show==TRUE) {bdl.col <- box.fill}
  else {bdl.col <- "white"}
  yname <- deparse(substitute(x1))
  if (is.null(Ylab))  Ylab <- yname
  y1 <- as.numeric(y1)
  y2 <- as.integer(y2)
  xmin = 0
  Ylim <- NULL
  grp.all = "1"
  gname = NULL
  rslt<-data.frame(MaximumDL=NA,Group=NA)

  if (sum(y2) > 0)    # not all data are detects
  { dlmax <- max(y1[y2 == 1])
  dltxt <- paste("Max DL=", signif(dlmax, 5), sep="")
  nonas <- na.omit(data.frame(y1, y2))         # omits NAs for the sake of the ros function
  y.nona <- nonas[,1]
  nd <- (nonas[,2] == 1)   # nd is vector of logical T/F
  if(sum(1-nonas[,2]) <= 2 )  {
    stop("Note: Data had fewer than 3 detects and so cannot be plotted")}
  try(if(sum(as.integer(nonas[,2])) == 0) stop("All data are detected.  Use the boxplot function."))

  if (LOG == TRUE)  {
    if (is.null(group) == TRUE) {    # log scale, no group
      xx = c(0.5, 1.5, 1.5, 0.5)
      y.ros <- suppressWarnings(ros(y.nona, nd))
      y.min <- log(min(y.ros$modeled))
      yy = c( y.min, y.min, log(dlmax), log(dlmax))
      Ylab = c(Ylab, "(natural logs)" )
      if (is.null(Ymax) == FALSE) { Ylim = c(y.min, log(Ymax))}
      if (minmax != TRUE) {
        # boxplot with separate outliers beyond 1.5*IQR
        graphics::boxplot(log(y.ros$modeled), ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplecol = "white", outcex = 0.8, main =Title, col=bxcol, ylim = Ylim,...)}
      else {
        #  boxplot drawn to max and min.  minmax = TRUE
        graphics::boxplot(log(y.ros$modeled), ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplewex = 0.2, range=0,  main =Title, col=bxcol, ylim = Ylim,...)}

      polygon (xx, yy, col = bdl.col, border=bdl.col)
      abline (h=y.min, col=bdl.col, lwd=8)
      abline(h = log(dlmax), col = dl.col,
             lty = dl.lty, lwd = dl.lwd)
      text(1.35, log(dlmax), labels = dltxt, pos=1, col=dl.col, cex=0.8)
      rslt<-data.frame(MaximumDL=dlmax, Group= "all")
    }

    else {              # log scale, with groups
      nonas <- na.omit(data.frame(y1, y2, group))  # omits NAs for the sake of the ros function
      y.nona <- nonas[,1]
      nd <- (nonas[,2] == 1)   # nd is vector of logical T/F
      y.grp <- nonas[,3]
      grp <- as.factor(y.grp)
      gname <- deparse(substitute(xgroup))
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
          warning("Note: One group had fewer than 3 detects and so cannot be plotted.")

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
      if (minmax != TRUE) {
        # boxplot with separate outliers beyond 1.5*IQR
        boxplot(log(y.all)~grp.all, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplecol = "white", outcex = 0.8, main =Title, col=bxcol, ylim = Ylim,...)}
      else {
        #  boxplot drawn to max and min. minmax = TRUE
        boxplot(log(y.all)~grp.all, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplewex = 0.2, range=0,  main =Title, col=bxcol, ylim = Ylim,...)}

      multiDL <- ifelse (sd(maxDL) != 0, TRUE, FALSE)
      if (multiDL == FALSE) {     # plot all with same max DL
        # Everything below the max DL is gray, not black
        polygon(xx, yy, col = bdl.col, border=bdl.col)
        abline(h = log(dlmax), col = dl.col,
               lty = dl.lty, lwd = dl.lwd)
        legend(dl.loc, legend = dltxt, bty="n", text.col=dl.col, cex=0.8)
        rslt[,1] <- dl.loc
        rslt[,2] <- t(glevels)
      }
      else {
        #        rslt<-data.frame(MaximumDL=NA,Group=NA)
        for (i in 1:gpnum) {xgrp <- c(i-0.5, i+0.5, i+0.5, i-0.5)   # different max DL per group
        ygrp = c(log(ymin.all), log(ymin.all), log(maxDL[i]), log(maxDL[i]))
        polygon(xgrp, ygrp, col = bdl.col, border = bdl.col)
        lx = c(i-0.45, i+0.45)
        ly = c(log(maxDL[i]), log(maxDL[i]))
        lines(lx, ly, col = dl.col, lty = dl.lty, lwd = dl.lwd)
        rslt[i,1]=maxDL[i]
        rslt[i,2]=levels(grp.all)[i]
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
      if (minmax != TRUE) {
        # boxplot with separate outliers beyond 1.5*IQR
        boxplot(y.ros$modeled, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplecol = "white", outcex = 0.8, main =Title, col=bxcol, ylim = Ylim,...)}
      else {
        #  boxplot drawn to max and min.  minmax = TRUE
        boxplot(y.ros$modeled, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplewex = 0.2, range=0,  main =Title, col=bxcol, ylim = Ylim,...)}

      # Everything below the max DL is gray, not black
      polygon(xx, yy, col = bdl.col, border=bdl.col)
      abline (h=y.min, col=bdl.col, lwd=8)
      abline(h=dlmax, col = dl.col, lty="longdash", lwd=2)
      text(1.35, dlmax, labels = dltxt, pos=3, col=dl.col, cex=0.8)
      rslt<-data.frame(MaximumDL=dlmax, Group= "all")
    }

    else {            #  LOG = FALSE, with groups
      nonas <- na.omit(data.frame(y1, y2, group))  # omits NAs for the sake of the ros function
      y.nona <- nonas[,1]
      nd <- (nonas[,2] == 1)   # nd is vector of logical T/F
      y.grp <- nonas[,3]

      grp <- as.factor(y.grp)
      gname <- deparse(substitute(xgroup))
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
          warning("Note: One group had fewer than 3 detects and so cannot be plotted.")

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
      if (minmax != TRUE) {
        # boxplot with separate outliers beyond 1.5*IQR
        boxplot(y.all~grp.all, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplecol = "white", outcex = 0.8, main =Title, col=bxcol, ylim = Ylim,...)}
      else {
        #  boxplot drawn to max and min. minmax = TRUE
        boxplot(y.all~grp.all, ylab = Ylab, xlab = Xlab, log="", whisklty = "solid", staplewex = 0.2, range=0,  main =Title, col=bxcol, ylim = Ylim,...)}

      multiDL <- ifelse (sd(maxDL) != 0, TRUE, FALSE)
      if (multiDL == FALSE) {      # all groups use same max DL
        polygon(xx, yy, col = bdl.col, border=bdl.col)
        abline (h=0, col=bdl.col, lwd=8)
        abline(h=dlmax, col = dl.col, lty="longdash", lwd=2)
        legend(dl.loc, legend = dltxt, bty="n", text.col=dl.col, cex=0.8)
      }
      else {
        #      rslt<-data.frame(MaximumDL=NA,Group=NA)
        for (i in 1:gpnum) {xgrp <- c(i-0.5, i+0.5, i+0.5, i-0.5)   # different max DL per group
        ygrp = c(ymin.all, ymin.all, maxDL[i], maxDL[i])
        polygon(xgrp, ygrp, col = bdl.col, border = bdl.col)
        lx = c(i-0.45, i+0.45)
        ly = c(maxDL[i], maxDL[i])
        lines(lx, ly, col = dl.col, lty="longdash", lwd=1.5)
        rslt[i,1]=maxDL[i]
        rslt[i,2]=levels(grp.all)[i]
        }
      }
    }
  }
  }   # end of when there are nondetects

  else   # when there are no nondetects
  { LOG <- ifelse (LOG, "y", "")
  if (is.null(group) == TRUE)    # no group
  {  boxplot(y1, na.action = na.omit, ylab = Ylab, col = bxcol, main = Title, log=LOG,...)
    rslt<-data.frame(MaximumDL="none", Group= "all")
  }
  else        # with groups
  {gname <- deparse(substitute(group))
  if (is.null(ordr) == FALSE)  {
    group = factor(group, levels = ordr)
    }
  glabs <- levels(group)
  boxplot(y1~group, na.action = na.omit, ylab = Ylab, xlab = gname, names = glabs, col = bxcol, main = Title, log=LOG,...)
  rslt[,1] <- "none"
  rslt[,2]= glabs
  }
  }

  if (!is.null(Hlines))
  { if (!is.null(group)) {xpos <- 0.6 + (0.05 * length(levels(grp)))}
    else {xpos = 0.6}
    for (j in 1:nrow(Hlines))
    { if (LOG == FALSE) {abline(h=Hlines[j,1], col=Hlines[j,2], lty = Hlines[j,3])
      text(x=as.vector(xpos), y=as.vector(as.numeric(Hlines[j,1])), labels = Hlines[j,4], pos=3, col=Hlines[j,2], cex=0.8, offset = 0.25)}
      else {abline(h=log(as.numeric(Hlines[j,1])), col=Hlines[j,2], lty = Hlines[j,3])
        text(x=as.vector(xpos), y=log(as.vector(as.numeric(Hlines[j,1]))), labels = Hlines[j,4], pos=3, col=Hlines[j,2], cex=0.8, offset = 0.25)}
    }
  }
  if(printstat==TRUE&sum(y2)>0){return(rslt)}else if(printstat==FALSE&sum(y2)>0){invisible(rslt)}
  else{warning("Dataset does not have any censored data.")}
}

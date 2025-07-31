## From NADA package (./summary.R), code revised for S3 standards

#' Summary Statistics for Censored Data
#'
#' Produces basic, and hopefully useful, summary statistics on censored data.
#'
#' This is a generic function with methods for numeric vectors and grouped data.
#'
#' @param obs A numeric vector of observations.
#' @param censored A logical vector indicating `TRUE` where an observation in `obs` is censored (a less-than value), and `FALSE` otherwise.
#' @param groups Optional grouping variable (factor or numeric) used to compute summaries by subset.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class `"censummary"` containing summary statistics,
#' or a list of such objects if grouped.
#'
#' @details The generic function dispatches to methods based on the presence
#' and type of `groups`:
#' \describe{
#'   \item{`censummary.default`}{Basic summary for numeric vectors.}
#'   \item{`censummary.factor`}{Summary grouped by a factor variable.}
#'   \item{`censummary.numeric`}{Summary grouped by a numeric grouping variable (converted to factor).}
#' }
#'
#' @references
#' Helsel, Dennis R. (2005). *Nondetects and Data Analysis: Statistics for Censored Environmental Data*. John Wiley and Sons, USA, NJ.
#'
#' @export
censummary <- function(obs, censored, groups = NULL, ...) {
  if (missing(groups)) {
    return(censummary.default(obs, censored, ...))
  } else if (is.factor(groups)) {
    return(censummary.factor(obs, censored, groups, ...))
  } else if (is.numeric(groups)) {
    return(censummary.numeric(obs, censored, groups, ...))
  } else {
    stop("Unsupported 'groups' type: must be missing, factor, or numeric.")
  }
}

#' @rdname censummary
#' @export
censummary.default <- function(obs, censored, ...) {
  ret <- cohn(obs, censored)
  ret$n <- length(obs)
  ret$n.cen <- sum(censored)

  props <- c(ret$n, ret$n.cen, pctCen(obs, censored), min(obs), max(obs))
  names(props) <- c("n", "n.cen", "pct.cen", "min", "max")

  limits <- data.frame(ret$limit, ret$C, ret$A, ret$P)
  colnames(limits) <- c("limit", "n", "uncen", "pexceed")

  structure(
    list(all = props, limits = limits),
    class = "censummary"
  )
}

#' @rdname censummary
#' @export
censummary.factor <- function(obs, censored, groups, ...) {
  ret <- lapply(levels(groups), function(i) {
    x <- obs[groups == i]
    xc <- censored[groups == i]
    censummary(x, xc, ...)
  })
  names(ret) <- levels(groups)
  structure(ret, class = "NADAList")
}

#' @rdname censummary
#' @export
censummary.numeric <- function(obs, censored, groups, ...) {
  censummary.factor(obs, censored, as.factor(groups), ...)
}


#' @export
#' @method print censummary
print.censummary <- function(x, ...) {
  for (i in names(x)) {
    cat(i, ":\n", sep = "")
    print(x[[i]])
    cat("\n")
  }
}


#' Summary Statistics for Censored Data using ROS, MLE, and K-M
#'
#' A convenience function that produces a comparative table of
#' summary statistics obtained using the `ros()`, `cenmle()`,
#' and `cenfit()` routines. These methods are Regression on
#' Order Statistics (ROS), Maximum Likelihood Estimation (MLE), and
#' Kaplan-Meier (K-M).

#' @param ... to pass argument
#'
#' @return A data frame with the summary statistics.
#'
#' @details If the data do not fulfill the criteria for the application of
#' any method, no summary statistics will be produced.
#'
#' @export
censtats <- function(...){
  stop("function in development")
}
# @param obs A numeric vector of observations.
# @param censored A logical vector indicating `TRUE` where an observation in `obs` is censored (a less-than value) and `FALSE` otherwise.
# censtats <- function(obs, censored) {
#   skm <- cenfit(obs, censored)
#   sros <- ros(obs, censored)
#   smle <- cenmle(obs, censored)
#
#   med <- c(median(skm), median(sros), median(smle))
#   sd <- c(sd(skm), sd(sros), sd(smle)[1])
#   mean <- c(mean(skm)[1], mean(sros), mean(smle)[1])
#
#   len <- c(length(obs), sum(censored), pctCen(obs, censored))
#   names(len) <- c("n", "n.cen", "pct.cen")
#   print(len)
#
#   data.frame(median = med, mean = mean, sd = sd,
#              row.names = c("K-M", "ROS", "MLE"))
# }


#' Produces a censored boxplot
#'
#' Draws a boxplot with the highest censoring threshold shown as a horizontal line.
#' Any statistics below this line are invalid and must be estimated using methods for censored data.
#'
#' @param obs A numeric vector of observations.
#' @param cen A logical vector indicating TRUE where an observation in `obs` is
#'   censored (a less-than value) and FALSE otherwise.
#' @param group A factor vector used for grouping `obs` into subsets (each group
#'   will be a separate box).
#' @param log A logical indicating if the y-axis should be in log units.
#'   Default is `TRUE`.
#' @param range A numeric value determining how far the plot whiskers extend
#'   from the box. If positive, the whiskers extend to the most extreme data
#'   point which is no more than `range` times the interquartile range from
#'   the box. The default is 0, which extends whiskers to the min and max.
#' @param ... Additional arguments passed to [graphics::boxplot()].
#'
#' @importFrom graphics boxplot
#'
#' @return The output of the default [graphics::boxplot()] method.
#'
#' @references
#' Helsel, Dennis R. (2005). *Nondetects and Data Analysis; Statistics for
#' censored environmental data*. John Wiley and Sons, Hoboken, NJ.
#'
#' @export
cenboxplot <- function(obs, cen, group, log = TRUE, range = 0, ...) {
  log <- if (log) "y" else ""

  if (missing(group)) {
    graphics::boxplot(ros(obs, cen), log = log, range = range, ...)
  } else {
    modeled <- groups <- character()
    for (i in levels(as.factor(group))) {
      mod <- suppressWarnings(
        ros(obs[group == i], cen[group == i])$modeled
      )
      grp <- rep(i, length(mod))
      modeled <- c(modeled, mod)
      groups <- c(groups, grp)
    }
    graphics::boxplot(as.numeric(modeled) ~ as.factor(groups),
            log = log, range = range, ...)
    ret <- data.frame(ros.model = as.numeric(modeled), group = groups)
  }

  abline(h = max(obs[cen]), col = "red", lty = 2)
  invisible(ret)
}



#' Produces a censored x-y scatter plot
#'
#' Draws an x-y scatter plot with censored values represented by
#' dashed lines spanning from the censored threshold to zero.
#'
#' @param x A numeric vector of observations.
#' @param xcen A logical vector indicating TRUE where an observation in `x` is
#'   censored (a less-than value) and FALSE otherwise.
#' @param y A numeric vector of observations.
#' @param ycen A logical vector indicating TRUE where an observation in `y` is
#'   censored (a less-than value) and FALSE otherwise.
#' @param log A character string specifying which axes are logarithmic:
#'   `"x"` for x-axis, `"y"` for y-axis, `"xy"` or `"yx"` for both axes,
#'   or `""` for none (default).
#' @param lty The line type for the lines representing censored-data ranges.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @importFrom graphics segments
#'
#' @references
#' Helsel, Dennis R. (2005). *Nondetects and Data Analysis; Statistics for
#' censored environmental data.* John Wiley and Sons, Hoboken, NJ.
#'
#' @export
cenxyplot <- function(x, xcen, y, ycen, log = "", lty = "dashed", ...) {
  plot(x, y, log = log, type = "n", ...)
  points(x[!ycen], y[!ycen], ...)

  x0 <- x1 <- x[ycen]
  y0 <- rep(cenpar.usr(log)[3], length(x[ycen]))
  y1 <- y[ycen]

  segments(x0, y0, x1, y1, lty = lty, ...)
}

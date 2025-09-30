#' @title Regression on Order Statistics (ROS)
#' @name ros
#'
#' @description Perform regression on order statistics for left-censored data.
#'
#' @param obs Numeric vector of observations or formula of the form `response ~ censor`, where `response` is numeric values and `censor` is a logical or binary indicator (TRUE if censored).
#' @param censored Logical vector of left-censored indicators.
#' @param data A `data.frame` containing the variables in the formula.
#' @param forwardT Name of transformation function (e.g., "log", "trueT").
#' @param reverseT Name of back-transformation function (e.g., "exp", "trueT").
#' @param na.action Function to handle missing values.
#' @param ... Additional arguments passed to the default method.
#'
#' @details Code for this function is originally from the NADA package developed by R. Lopaka Lee and Dennis Helsel.
#' By default, ros performs a log transformation prior to, and after operations over the data.
#' This can be changed by specifying a forward and reverse transformation function using the forwardT
#' and reverseT parameters. No transformation will be performed if either forwardT or reverseT are set to NULL.
#'
#' The procedure first computes the Weibull-type plotting positions of the combined uncensored
#' and censored observations using a formula designed for multiply-censored data (see hc_ppoints).
#' A linear regression is formed using the plotting positions of the uncensored observations and
#' their normal quantiles. This model is then used to estimate the concentration of the censored
#' observations as a function of their normal quantiles. Finally, the observed uncensored values
#' are combined with modeled censored values to corporately estimate summary statistics of the entire
#' population. By combining the uncensored values with modeled censored values, this method is more
#'  resistant of any non-normality of errors, and reduces any transformation errors that may be incurred.
#'
#' @return A list with:
#' \describe{
#'   \item{`modeled`}{Numeric vector of uncensored + imputed censored values.}
#'   \item{`modeled.censored`}{Imputed values only for censored observations.}
#'   \item{`uncensored`}{Original uncensored values.}
#'   \item{`censored`}{Original censored values.}
#'   \item{`censored.ranks`}{Censored ranks used in estimation.}
#'   \item{`uncensored.ranks`}{Uncensored ranks used in estimation.}
#'   \item{`model`}{Fitted linear model object.}
#' }
#' @export
#' @importFrom stats model.frame qnorm predict.lm
#'
#' @examples
#' df <- data.frame(
#'   conc = c(0.2, 0.5, 1.0, 0.4, 2.0, 0.3),
#'   censored = c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
#' )
#' ros(conc ~ censored, data = df)
#'

ros <- function(obs, censored = NULL,data=NULL, forwardT = "log", reverseT = "exp",
                 na.action = getOption("na.action")) {

  if (inherits(obs, "formula")) {
    mf <- model.frame(obs, data = data)

    if (ncol(mf) != 2) {
      stop("Formula must be of the form response ~ censor")
    }

    obs <- mf[[1]]
    censored <- mf[[2]]
  }
  if (is.null(obs) || is.null(censored)) stop("Both x and y must be provided")
  if (is.null(forwardT) || is.null(reverseT)) {
    forwardT <- reverseT <- "trueT"
  } else {
    if (!exists(forwardT)) stop("Cannot find forwardT function: ", forwardT)
    if (!exists(reverseT)) stop("Cannot find reverseT function: ", reverseT)
  }

  if (anyNA(obs) || anyNA(censored)) {
    if (is.null(na.action)) na.action <- "na.exclude"
    obs <- do.call(na.action, list(obs))
    censored <- do.call(na.action, list(censored))
  }

  if (mean(censored) > 0.8) {
    warning("More than 80% of data is censored.")
  }

  ix <- which(obs > max(obs[!censored]))
  if (length(ix)) {
    warning("Dropped censored values that exceed max of uncensored values.")
    obs <- obs[-ix]
    censored <- censored[-ix]
  }

  ix <- order(obs)
  obs <- obs[ix]
  censored <- censored[ix]

  pp <- hc_ppoints(obs, censored, na.action)
  pp.nq <- qnorm(pp[!censored])
  obs.transformed <- get(forwardT)(obs[!censored])
  lm_model <- lm(obs.transformed ~ pp.nq, na.action = na.action)

  modeled <- obs
  new_dat <- data.frame(pp.nq = qnorm(pp[censored]))
  modeled[censored] <- get(reverseT)(predict.lm(object = lm_model, newdata = new_dat))

  result <- list(
    obs = obs,
    censored = censored,
    modeled = modeled,
    pp = pp,
    model = lm_model,
    forwardT = forwardT,
    reverseT = reverseT
  )
  class(result) <- "ros"
  return(result)
}


#' @describeIn ros A wrapper for `ros()` with simplified argument handling.
#' @export
cenros <- function(x, ...) {
  ros(x, ...)
}

#' Class `"ros"`
#'
#' An S3 class returned by the `ros()` function, representing the result of regression on order statistics for left-censored data.
#'
#' @name ros-class
#' @docType class
#' @format An object of class `"ros"` with elements:
#' \describe{
#'   \item{`obs`}{The original observations.}
#'   \item{`censored`}{Logical vector indicating which observations were censored.}
#'   \item{`modeled`}{Vector with uncensored and imputed values.}
#'   \item{`pp`}{Plotting positions.}
#'   \item{`model`}{The fitted linear model object.}
#'   \item{`forwardT`}{Transformation function used.}
#'   \item{`reverseT`}{Back-transformation function used.}
#' }
#' @keywords internal
NULL


#' @importFrom stats median
#' @export
print.ros <- function(x, ...) {
  n <- length(x$modeled)
  n.cen <- sum(x$censored)
  median.cen <- median(x$modeled)
  mean.cen <- mean(x$modeled)
  sd.cen <- sd(x$modeled)

  ret  <-  c(n, n.cen, median.cen, mean.cen, sd.cen)
  names(ret)  <- c("n", "n.cen", "median", "mean", "sd")

  print(ret,...)
  invisible(ret)
}


#' @export
summary.ros <- function(object, plot = FALSE, ...) {
  ret <- summary(object$model)
  if (plot) {
    plot(object$model, ...)
  }
  return(ret)
}

#' @rdname ros
#' @method plot ros
#' @export
#' @param x ros2 model object
#' @param plot.censored default = FALSE, if set to true it will also plot censored data
#' @param lm.line will plot linear model line
#' @param grid will add grid
#' @param ylab default is "Value" but custom text can be added
#' @param pch default set to 16, codes consistent with `points` and `plot` functions
#' @param ... arguments passed to plot function
#' @importFrom graphics axTicks axis
#' @importFrom stats qnorm

plot.ros <- function(x,
                      plot.censored = FALSE,
                      lm.line = TRUE,
                      grid = TRUE,
                      ylab = "Value",
                      pch = 16,
                      ...) {
  # Extract modeled values and censored indicators
  modeled <- x$modeled
  censored <- x$censored
  pp <- x$pp

  # Normal quantiles for plotting positions
  pp.nq.uncen <- qnorm(pp[!censored])
  pp.nq.cen <- qnorm(pp[censored])

  y.uncen <- modeled[!censored]
  y.cen <- modeled[censored]

  # Plot limits
  ylim <- range(c(y.uncen, y.cen), na.rm = TRUE)
  xlim <- range(c(pp.nq.uncen, pp.nq.cen), na.rm = TRUE)

  # Use log scale on y axis if forwardT is "log"
  logy <- ifelse(!is.null(x$forwardT) && x$forwardT == "log", "y", "")

  # Base plot with uncensored points
  plot(pp.nq.uncen, y.uncen,
       xlim = xlim, ylim = ylim,
       xlab = "Normal Quantiles",
       ylab = ylab,
       pch = pch,
       log = logy,
       ...)

  # Add censored points if requested
  if (plot.censored && length(y.cen) > 0) {
    points(pp.nq.cen, y.cen, pch = 1, col = "red", ...)
  }

  # Add regression line if requested
  if (lm.line) {
    lines.ros(x, col = "blue", lwd = 2)
  }

  # Add grid lines if requested
  if (grid) {
    abline(h = axTicks(2), lty = "dotted", col = "gray")
    atv <- qnorm(c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
    abline(v = atv, lty = "dotted", col = "gray")

    # Add top axis for Percent Chance of Exceedance
    axis(3, at = atv, labels = rev(c("5", "10", "25", "50", "75", "90", "95")), las = 2)
    mtext("Percent Chance of Exceedance", side = 3, line = 3)
  }
}

#' @export
#' @importFrom stats qnorm
lines.ros <- function(x, col = "blue", lwd = 2, ...) {
  pp <- x$pp
  # Normal quantiles covering the plotting positions range
  minmax.nq <- qnorm(range(pp))

  # Predict modeled values for those quantiles
  preds <- predict(x$model, newdata = data.frame(pp.nq = minmax.nq))

  # Apply reverse transform if present
  revT <- x$reverseT
  if (!is.null(revT) && exists(revT)) {
    preds <- get(revT)(preds)
  }

  lines(minmax.nq, preds, col = col, lwd = lwd, ...)
}

#' @export
as.data.frame.ros <- function(x, ...) {
  data.frame(
    obs = x$obs,
    censored = x$censored,
    pp = x$pp,
    modeled = x$modeled
  )
}

#' Compute plotting positions for censored and uncensored data
#'
#' @param obs A numeric vector of observations.
#' @param censored A logical vector indicating which values are censored.
#' @param na.action A function to handle NA values (default from options).
#' @return A numeric vector of plotting positions.
#' @importFrom stats ppoints
#' @keywords internal

hc_ppoints <- function(obs, censored, na.action = getOption("na.action")) {
  # from hc.ppoints
  if (!is.logical(censored)) stop("censored indicator must be logical vector!\n")

  if (anyNA(obs)||anyNA(censored)) {
    if (is.null(na.action)) na.action <- "na.omit"
    obs  <- do.call(na.action, list(obs))
    censored  <- do.call(na.action, list(censored))
  }

  pp  <- numeric(length(obs))

  if (!any(censored)) pp <- ppoints(obs)
  else
  {
    cn  <- cohn(obs, censored)
    pp[!censored]  <- hc_ppoints_uncen(obs, censored, cn, na.action)
    pp[censored]   <- hc_ppoints_cen(obs, censored, cn, na.action)
  }

  return(pp)
}

#' @title Plotting Positions for Uncensored Observations (Cohn Method)
#' @description
#' Computes plotting positions for uncensored observations based on the
#' Cohn method grouping. This is used in hydrologic statistics and
#' censoring methods.
#'
#' @param obs A numeric vector of observed values.
#' @param censored A logical vector indicating which observations are censored (`TRUE` for censored, `FALSE` otherwise).
#' @param cn An optional list containing Cohn grouping information (usually from [cohn()]); if missing, it will be computed internally.
#' @param na.action A function to handle missing values (default is `getOption("na.action")`).
#'
#' @return A numeric vector of plotting positions corresponding to the uncensored observations.
#'
#' @keywords internal
hc_ppoints_uncen <- function(obs, censored, cn = NULL, na.action = getOption("na.action")) {
  ## From NADA hc.ppoints.uncen
  if (!is.logical(censored)) stop("censored indicator must be logical vector!\n")
  if (anyNA(obs)||anyNA(censored)) {
    if (is.null(na.action)) na.action  <- "na.omit"
    obs  <-  do.call(na.action, list(obs))
    censored <- do.call(na.action, list(censored))
  }

  if (missing(cn)) { cn <- cohn(obs, censored) }
  #cn = cohn(obs, censored)

  nonzero <- (cn$A != 0)
  A <- cn$A[nonzero]
  B <- cn$B[nonzero]
  P <- cn$P[nonzero]
  limit <- cn$limit[nonzero]

  pp  <- numeric()
  n <- length(limit)
  for (i in 1:n){
    R <- 1:A[i]

    k <- P[i+1]
    if (is.na(k)) k <- 0

    for (r in seq_along(R)){
      pp  <- c(pp, (1 - P[i]) + ((P[i] - k) * R[r])/(A[i] + 1))
    }
  }
  return(pp)
}


#' @title Plotting Positions for Censored Observations (Cohn Method)
#' @description
#' Computes plotting positions for censored observations using the
#' Cohn method grouping. This is primarily for hydrologic data with left-censoring.
#'
#' @param obs A numeric vector of observed values.
#' @param censored A logical vector indicating which observations are censored (`TRUE` for censored, `FALSE` otherwise).
#' @param cn An optional list containing Cohn grouping information (usually from [cohn()]); if missing, it will be computed internally.
#' @param na.action A function to handle missing values (default is `getOption("na.action")`).
#'
#' @return A numeric vector of plotting positions corresponding to the censored observations.
#'
#' @keywords internal
hc_ppoints_cen <- function(obs, censored, cn = NULL, na.action = getOption("na.action")) {
  ## From NADA hc.ppoints.cen

  if (!is.logical(censored)) stop("censored indicator must be logical vector!\n")
  if (anyNA(obs)||anyNA(censored)) {
    if (is.null(na.action)) na.action  <- "na.omit"
    obs <- do.call(na.action, list(obs))
    censored  <-  do.call(na.action, list(censored))
  }

  if (missing(cn)) { cn <- cohn(obs, censored) }

  C <- cn$C
  P <- cn$P
  limit <- cn$limit

  if (P[1] == 1)
  {
    C <- C[-1]
    P <- P[-1]
    limit <- limit[-1]
  }

  pp  <- numeric()
  for (i in 1:length(limit))
  {
    c.i <- C[i]
    for (r in 1:c.i)
    {
      pp <- c(pp, (1 - P[i]) * r/(c.i + 1))
    }
  }
  return(pp)
}



#' Calculate Cohn Numbers
#'
#' Computes the A, B, C, P values used in probability plotting for censored data.
#'
#' @param obs A numeric vector of observations.
#' @param censored A logical vector indicating which observations are censored.
#' @param na.action Function to handle missing values (default is getOption("na.action")).
#' @return An object of class `"cohn"` containing the Cohn statistics.
#' @keywords internal
#' @examples
#'
#' obs <- c(1, 2, 3, 4, 5, 6, 7)
#' cens <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
#' NADA2:::cohn(obs, cens)
#'
cohn <- function(obs, censored, na.action = getOption("na.action")) {
  # improved speed from vectorization
  if (!is.logical(censored)) stop("censored must be a logical vector.")

  # Apply na.action
  if (anyNA(obs) || anyNA(censored)) {
    if (is.null(na.action)) na.action <- na.omit
    obs <- do.call(na.action, list(obs))
    censored <- do.call(na.action, list(censored))
  }

  uncen <- obs[!censored]
  cen   <- obs[censored]

  limit <- sort(unique(cen))
  if (any(uncen < limit[1])) limit <- c(0, limit)

  n <- length(limit)

  # Compute A, B, C vectorized
  A <- vapply(1:n, function(j) {
    if (j == n) sum(uncen >= limit[j])
    else sum(uncen >= limit[j] & uncen < limit[j + 1])
  }, numeric(1))

  B <- vapply(1:n, function(j) sum(obs <= limit[j]) - sum(uncen == limit[j]), numeric(1))
  C <- vapply(1:n, function(j) sum(cen == limit[j]), numeric(1))

  # === Numerically stable P ===
  frac <- A / (A + B)
  P <- numeric(n)
  P[n] <- frac[n]
  if (n > 1) {
    P[(n-1):1] <- Reduce(
      function(nextP, j) frac[j] * (1 - nextP) + nextP,
      (n-1):1,
      init = P[n],
      accumulate = TRUE
    )[-1]
  }

  structure(list(A = A, B = B, C = C, P = P, limit = limit),
            class = "cohn")
}
#' @export
print.cohn <- function(x, ...) {
  cat("Cohn Probabilistic Plotting Statistics\n")
  print(data.frame(limit = x$limit, A = x$A, B = x$B, C = x$C, P = x$P), row.names = FALSE)
  invisible(x)
}
#' @export
summary.cohn <- function(object, ...) {
  cat("Summary of Cohn Probabilistic Plotting Statistics:\n")
  print(data.frame(
    Limit = object$limit,
    A = object$A,
    B = object$B,
    C = object$C,
    P = round(object$P, 4)
  ), row.names = FALSE)
  invisible(object)
}


#' Compute Kendall's Tau Correlation and ATS Line for Censored Data
#'
#' Computes Kendall's tau for singly (y only) or doubly (x and y) censored data.
#' Also computes the Akritas-Theil-Sen (ATS) nonparametric regression line,
#' with the Turnbull estimate of the intercept.
#'
#' @param y A numeric vector of observations or a formula.
#' @param ycen A logical vector indicating `TRUE` where an observation in `y` is
#'   censored (a less-than value) and `FALSE` otherwise. Can be omitted for the case
#'   where `y` is not censored.
#' @param x A numeric vector of observations.
#' @param xcen default is NULL for trend analysis purposes. If included, a
#'   logical vector indicating `TRUE` where an observation in `x` is
#'   censored (a less-than value) and `FALSE` otherwise.
#' @details
#' If using the formula interface, the `ycen`, `x`, and `xcen`
#' arguments are not specified. All information is instead provided through a
#' formula via the `y` argument. The formula must use a `Cen` object
#' as the response (on the left-hand side of `~`), and predictors (optional)
#' on the right-hand side separated by `+`. See examples below.
#'
#' Kendall's tau is a nonparametric correlation coefficient that measures
#' monotonic association between `y` and `x`. For left-censored data,
#' concordant and discordant pairs are evaluated wherever possible. For example,
#' with increasing `x` values, a change in `y` from <1 to 10 is considered
#' an increase (concordant), while a change from <1 to 0.5 or <5 is considered a tie.
#'
#' The ATS line is defined as the slope resulting in Kendall's tau of 0
#' between residuals `(y - slope * x)` and `x`. The routine uses
#' an iterative bisection search to find this slope. The intercept is the
#' median residual, computed using the Turnbull estimate for interval-censored
#' data as implemented in the \pkg{Icens} package.
#'
#' @return A list containing:
#' \item{tau}{Kendall's tau correlation coefficient.}
#' \item{slope}{The estimated slope from the ATS line.}
#' \item{p.value}{P-value for testing the null hypothesis of no association.}
#'
#' @references
#' Helsel, D. R. (2005). *Nondetects and Data Analysis: Statistics for Censored Environmental Data*.
#' John Wiley & Sons.
#'
#' Akritas, M. G., Murphy, S. A., & LaValley, M. P. (1995).
#' The Theil-Sen Estimator with Doubly Censored Data and Applications to Astronomy.
#' *Journal of the American Statistical Association*, 90, 170–177.
#'
#' @examples
#' # Both y and x are censored
#' data(PbHeron)
#' with(PbHeron,cenken(Blood,BloodCen,Kidney,KidneyCen))
#'
#' # x is not censored
#' data(TCEReg)
#' with(TCEReg, cenken(log(TCEConc), TCECen, PopDensity))
#'
#' # Synthetic time-series with trend analysis
#' set.seed(123)
#'
#' ## Parameters
#' n <- 15  # 15 years of data
#' time <- 1:n
#'
#' ## Components
#' trend <- 0.235 * time
#' noise <- rnorm(n, mean = 5, sd = 1.5)
#'
#' syn_dat <- data.frame(Yr = 1989 + time, value = trend + noise)
#' syn_dat$censored <- syn_dat$value < quantile(syn_dat$value, 0.2)
#'
#' with(syn_dat,cenken(value,censored,Yr))
#' \dontrun{
#' plot(value~Yr,syn_dat,pch=21,bg=ifelse(syn_dat$censored==TRUE,"red","blue",cex=1.5))
#' abline(h=quantile(syn_dat$value, 0.2),lty=2,col="red")
#' }


#'
#' @keywords regression
#' @export

cenken <- function(y, ycen = NULL, x = NULL, xcen = NULL) {
  UseMethod("cenken")
}

#' @export
cenken.default <- function(y, ycen, x,xcen = NULL) {
  if (is.null(xcen)) {
    xcen <- rep(FALSE, length(x))
  } else if (length(xcen) != length(x)) {
    stop("Length of 'xcen' must match length of 'x'")
  }

  kendallATS(y, ycen, x, xcen)

}

#' @export
cenken.formula <- function(y, ycen = NULL, x = NULL,xcen = NULL) {
  stop("Only default interface is currently supported in this S3 implementation.")
}

#' @export
print.cenken <- function(x, ...) {
  cat("Censored Kendall's Tau Regression:\n")
  cat("Slope     :", x$slope, "\n")
  cat("Intercept :", x$intercept, "\n")
  cat("Tau       :", x$tau, "\n")
  cat("p-value   :", x$p, "\n")
  invisible(x)
}

#' @export
lines.cenken <- function(x, ...) {
  abline(a = x$intercept, b = x$slope, ...)
}

#' @keywords internal
kendallATS <- function(y, ycen, x, xcen, tol = 1e-7, iter = 1e3) {
  # Apply small offset to y for tied values when censored
  y <- y - (min(y) / 1000) * ycen

  # Root-finding algorithm to maximize Kendall's tau
  iter_s <- function(lb, ub, x, xcen, y, ycen, tol, iter) {
    step_s <- function(lb, ub, x, xcen, y, ycen) {
      b <- c(lb, (lb + ub) / 2, ub)
      res <- sapply(b, function(slope) y - slope * x)
      s_vals <- apply(res, 2, ktau_s, x = x, xcen = xcen, ycen = ycen)
      list(b = b, s = s_vals)
    }

    bs <- step_s(lb, ub, x, xcen, y, ycen)

    for (i in seq_len(iter)) {
      b <- bs$b
      s <- bs$s

      if (s[1] * s[2] <= 0) {
        s[3] <- s[2]
        b[3] <- b[2]
      } else {
        s[1] <- s[2]
        b[1] <- b[2]
      }

      if (s[2] == 0 || abs(b[3] - b[1]) <= tol) break
      bs <- step_s(b[1], b[3], x, xcen, y, ycen)
    }

    bs
  }

  # Get initial slope bounds
  k_bounds <- ktau_b(x, xcen, y, ycen)

  # Find slope using iterative refinement
  bs <- iter_s(k_bounds[1], k_bounds[2], x, xcen, y, ycen, tol, iter)
  lbs <- iter_s(bs$b[1], bs$b[2], x, xcen, y, ycen, tol, iter)
  ubs <- iter_s(bs$b[2], bs$b[3], x, xcen, y, ycen, tol, iter)

  # Midpoint of lower and upper bounds as final slope
  slope <- 0.5 * (ubs$b[2] + lbs$b[2])

  # Compute intercept using Turnbull estimator
  intercept <- turnbull(y, ycen, x, xcen, slope)

  # Final tau and p-value
  tau_result <- ktau_p(x, xcen, y, ycen)

  # Return result as class "cenken"
  structure(list(
    slope = slope,
    intercept = intercept,
    tau = tau_result$tau,
    p = tau_result$p
  ), class = "cenken")
}


#' @keywords internal
turnbull <- function(y, ycen, x, xcen, slope, tol = .Machine$double.eps) {
  resid <- y[!xcen] - slope * x[!xcen]
  A <- cbind(L = resid - (y[!xcen] * ycen[!xcen]), R = resid)
  em <- EM(A, tol = tol)
  surv <- rev(cumsum(rev(em$pf)))[-1]
  int <- em$intmap[2,][min(which(surv <= 0.5))]
  return(int)
}

#' @keywords internal
ktau_b <- function(x, xcen, y, ycen) {
  cx <- !xcen
  cy <- !ycen
  slopes <- unlist(lapply(seq_along(x), function(i) (y[i] - y[1:i]) / (x[i] - x[1:i])))
  return(range(slopes[is.finite(slopes)]))
}

#' @keywords internal
ktau_s <- function(x, xcen, y, ycen) {
  xy1 <- sign(outer(y, y, "-"))
  xy2 <- xy1 * sign(outer(x, x, "-"))
  itot <- 0.5 * sum((outer(ycen, -ycen, "-") == 0) * xy2)
  xy3 <- outer(ycen, ycen, "-")
  xy4 <- outer(y, y, "<=")
  itot <- itot + sum(((xy3 * xy4) == 1) * xy2)
  return(itot)
}

#' @keywords internal
ktau_p <- function(x, xcen, y, ycen) {
  xx <- x
  cx <- xcen
  yy <- y
  cy <- ycen
  n  <- length(xx)

  delx <- min(diff(sort(unique(xx))))/1000
  dely <- min(diff(sort(unique(yy))))/1000
  dupx <- xx - delx * cx
  diffx <- outer(dupx, dupx, "-")
  diffcx <- outer(cx, cx, "-")
  xplus <- outer(cx, -cx, "-")
  dupy <- yy - dely * cy
  diffy <- outer(dupy, dupy, "-")
  diffcy <- outer(cy, cy, "-")
  yplus <- outer(cy, -cy, "-")
  signyx <- sign(diffy * diffx)
  tt <- (sum(1 - abs(sign(diffx))) - n)/2
  uu <- (sum(1 - abs(sign(diffy))) - n)/2
  cix <- sign(diffcx) * sign(diffx)
  cix <- ifelse(cix <= 0, 0, 1)
  tt <- tt + sum(cix)/2
  signyx <- signyx * (1 - cix)
  ciy <- sign(diffcy) * sign(diffy)
  ciy <- ifelse(ciy <= 0, 0, 1)
  uu <- uu + sum(ciy)/2
  signyx <- signyx * (1 - ciy)
  xplus <- ifelse(xplus <= 1, 0, 1)
  yplus <- ifelse(yplus <= 1, 0, 1)
  diffx <- abs(sign(diffx))
  diffy <- abs(sign(diffy))
  tplus <- xplus * diffx
  uplus <- yplus * diffy
  tt <- tt + sum(tplus)/2
  uu <- uu + sum(uplus)/2

  itot <- sum(signyx * (1 - xplus) * (1 - yplus))
  kenS <- itot/2
  tau <- (itot)/(n * (n - 1))

  J <- n * (n - 1)/2
  taub <- kenS/(sqrt(J - tt) * sqrt(J - uu))

  varS <- n * (n - 1) * (2 * n + 5)/18

  intg <- 1:n
  dupx <- xx - delx * cx
  dupy <- yy - dely * cy
  dorder <- order(dupx)
  dxx <- dupx[dorder]
  dcx <- cx[dorder]
  dorder <- order(dupy)
  dyy <- dupy[dorder]
  dcy <- cy[dorder]
  tmpx <- dxx - intg * (1 - dcx) * delx
  tmpy <- dyy - intg * (1 - dcy) * dely
  rxlng <- rle(rank(tmpx))$lengths
  nrxlng <- table(rxlng)
  rxlng <- as.integer(names(nrxlng))
  x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
  x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
  x3 <- nrxlng * rxlng * (rxlng - 1)
  rylng <- rle(rank(tmpy))$lengths
  nrylng <- table(rylng)
  rylng <- as.integer(names(nrylng))
  y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
  y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
  y3 <- nrylng * rylng * (rylng - 1)
  delc <- (sum(x1) + sum(y1))/18 - sum(x2) * sum(y2)/(9 * n * (n - 1) * (n - 2)) - sum(x3) * sum(y3)/(2 * n * (n - 1))
  x4 <- nrxlng * (rxlng - 1)
  y4 <- nrylng * (rylng - 1)
  tmpx <- intg * dcx - 1
  tmpx <- ifelse(tmpx < 0, 0, tmpx)
  nrxlng <- sum(tmpx)
  rxlng <- 2
  x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
  x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
  x3 <- nrxlng * rxlng * (rxlng - 1)
  tmpy <- intg * dcy - 1
  tmpy <- ifelse(tmpy < 0, 0, tmpy)
  nrylng <- sum(tmpy)
  rylng <- 2
  y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
  y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
  y3 <- nrylng * rylng * (rylng - 1)
  deluc <- (sum(x1) + sum(y1))/18 - sum(x2) * sum(y2)/(9 * n * (n - 1) * (n - 2)) - sum(x3) * sum(y3)/(2 * n * (n - 1)) - (sum(x4) + sum(y4))
  dxx <- dxx - intg * dcx * delx
  dyy <- dyy - intg * dcy * dely
  rxlng <- rle(rank(dxx))$lengths
  nrxlng <- table(rxlng)
  rxlng <- as.integer(names(nrxlng))
  x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
  x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
  x3 <- nrxlng * rxlng * (rxlng - 1)
  rylng <- rle(rank(dyy))$lengths
  nrylng <- table(rylng)
  rylng <- as.integer(names(nrylng))
  y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
  y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
  y3 <- nrylng * rylng * (rylng - 1)
  delu <- (sum(x1) + sum(y1))/18 - sum(x2) * sum(y2)/(9 * n * (n - 1) * (n - 2)) - sum(x3) * sum(y3)/(2 * n * (n -  1))

  varS <- varS - delc - deluc - delu

  p.val <- 2 * (1 - pnorm((abs(kenS - sign(kenS)))/sqrt(varS)))

  return(list(tau=tau, p=p.val))
}


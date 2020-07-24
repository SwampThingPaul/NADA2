#' Computes confidence intervals on regression on order statistics (ROS) mean
#'
#' @description Regression on Ordered Statistics (ROS) is designed to evaluate censored analytical chemistry data. This function evaluates ROS results and computes the confidence interval around the mean.
#' @param cenros.out an ROS model output object (see details)
#' @param conf Confidence coefficient of the interval (Default is 95%)
#' @return Prints a lower (LCL) and upper (UCL) confidence interval based on the `conf` provided (Default is 95%)
#' @details
#' This function uses an ROS model output based on the `ros` funtion in the `NADA` package.
#'
#' For more detail on ROS modeling see the `ros` help file (`?NADA::ros`).
#'
#' For implementation of `ROSci(...)` see the examples below.
#'
#' @export
#'
#' @references
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
#'
#' Lee, L., Helsel, D., 2005. Statistical analysis of water-quality data containing multiple detection limits: S-language software for regression on order statistics. Computers & Geosciences 31, 1241–1248. <https://doi.org/10.1016/j.cageo.2005.03.012>
#'
#' Zhou, X.-H., Gao, S., 1997. Confidence Intervals for the Log-Normal Mean. Statistics in Medicine 16, 783–790. <https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7<783::AID-SIM488>3.0.CO;2-2>
#'
#' @seealso [NADA::ros]
#'
#' @examples
#' library(NADA) #For example data
#'
#' data(HgFish)
#' myros <- NADA::ros(HgFish$Hg,HgFish$HgCen)
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
    cat(cinames, "\n")
    cat(cilog, sep = "  ")}

  else { # for a gaussian distribution;
    tstat <- qt(p,(n-1))
    halfw <- tstat*(sd(cenros.out)/sqrt(n))
    ciNorm <- c(mean(cenros.out) - halfw, mean(cenros.out) + halfw)
    cat(cinames, "\n")
    cat(ciNorm, sep = "  ")
  }

}

#' Computes confidence intervals on regression on order statistics (ROS) mean
#'
#' @description Uses ROS model output from the `NADA` package and computes the Zhou and Gao 1997 modified Cox’s method two-sided confidence interval around the mean for a lognormal distribution.  Computes a t-interval for a gaussian ROS model output.
#' @param cenros.out an ROS model output object (see details)
#' @param conf Confidence coefficient of the interval (Default is 0.95)
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @return Prints a lower (LCL) and upper (UCL) confidence interval based on the `conf` provided (Default is 95%)
#' @details
#' This function uses an ROS model output based on the `ros` function in the `NADA` package.  The lognormal distribution is the default for the NADA package but a gaussian distribution is optional here.
#' For more detail on ROS modeling see the `ros` help file (`?NADA::ros`).
#'
#' For implementation of `ROSci(...)` see the examples below.
#' @importFrom stats qt sd
#' @export
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Lee, L., Helsel, D., 2005. Statistical analysis of water-quality data containing multiple detection limits: S-language software for regression on order statistics. Computers & Geosciences 31, 1241–1248. \doi{https://doi.org/10.1016/j.cageo.2005.03.012}
#'
#' Zhou, X.-H., Gao, S., 1997. Confidence Intervals for the Log-Normal Mean. Statistics in Medicine 16, 783–790. \doi{https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7<783::AID-SIM488>3.0.CO;2-2}
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
#' mean(myros$modeled)
#'
#' # 95% CI around the ROS mean
#' ROSci(myros)


ROSci <- function(cenros.out, conf=0.95,printstat=TRUE) {
  p <- 1-((1-conf)/2)
  n <- length(cenros.out$obs)
  c <- as.character(100*conf)

  if (cenros.out$forwardT == "log") {
    # For the default lognormal distribution, Cox’s method is used.
    scale <- as.vector(cenros.out[1]$coefficients[2])
    #  coefficient #2 is the scale, the sd of the logs
    bhat <- log(mean(cenros.out$modeled))
    #  Mean log plus one-half scale^2 ;  the mean transformed to log units
    gamz <- qt(p,(n-1)) * sqrt((scale^2/n) + (((0.5)*scale^4)/(n-1)))
    #  Zhou and Gao 1997 modified Cox’s method for the CI of a lognormal distribution
    cilog <- c(exp(bhat - gamz), exp(bhat + gamz))
    rslt<-data.frame(LCL=cilog[1],UCL=cilog[2])

    if(printstat==TRUE){
    cat("Assuming a lognormal distribution", "\n")
    print(rslt)
    }

    }

  else { # for a gaussian distribution;
    tstat <- qt(p,(n-1))
    halfw <- tstat*(sd(cenros.out)/sqrt(n))
    ciNorm <- c(mean(cenros.out$modeled) - halfw, mean(cenros.out$modeled) + halfw)

    rslt<-data.frame(LCL=ciNorm[1],UCL=ciNorm[2])
    if(printstat==TRUE){
    cat("Assuming a gaussian distribution", "\n")
    print(rslt)
    }

  }
return(invisible(rslt))
}

#' Prediction interval for censored data
#'
#' @description Computes prediction intervals for censored data assuming lognormal, gamma and normal distributions.
#' @param y.var The column of y (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detectionlimit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param pi.type Designation of either a `“two-sided”` interval (default) or a 1-sided `“upper”` or 1-sided `“lower”` interval.
#' @param conf Confidence coefficient of the interval, 0.95 (default).
#' @param newobs The number of new observations to be contained in the interval.
#' @param method character string specifying the method of estimation. Default is `mle` (maximum likelihood).
#' @keywords prediction interval
#' @export
#' @importFrom EnvStats elnormCensored predIntLnorm enormCensored predIntNorm
#' @return A table of prediction limits based on user provided confidence coefficient (`conf`) and prediction invterval type (`pi.type`)
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
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


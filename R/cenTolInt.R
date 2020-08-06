#' Upper Tolerance interval for censored data
#'
#' @description Computes a one-sided upper tolerance interval for censored data assuming lognormal, gamma and normal distributions.
#' @param y.var The column of y (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detectionlimit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param conf Confidence coefficient of the interval, 0.95 (default).
#' @param cover Coverage, the percentile probability above which the tolerance interval is computed.  The default is 90, so a tolerance interval will be computed above the 90th percentile of the data.
#' @param method.fit The method used to compute the parameters of the distribution.  The default is maximum likelihood (`“mle”`). The alternative is robust ROS (`“rROS”`).
#' @importFrom EnvStats elnormCensored predIntLnorm enormCensored predIntNorm eqlnormCensored eqnormCensored
#' @importFrom fitdistrplus fitdistcens
#' @export
#'
#' @return
#' Plot of Empirical and theoretical CDFs with BIC values provided
#'
#' Prints the percentile (`cover`) and upper tolerance limit (`conf`) for lognormal, normal and approximated gamma distributions.
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Helsel, D.R., 2005. Nondetects and Data Analysis: Statistics for Censored Environmental Data, 1st ed. John Wiley and Sons, USA, N.J.
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
results <- data.frame(Distribution, pct.best, ti.best, bic.3)
names(results) <- c("Distribution", pct.text, ti.text, "BIC")

bic.best = min(bic.3)
dist.best = floor(bic.best/bic.3)
pct.best = max(pct.best*dist.best)
ti.best = max(ti.best*dist.best)

cenCompareCdfs (y.var, cen.var, Yname = nameofy, dist3 = "norm")
abline (h=cover, lty = "dotted", col = "black")
return (results)
}

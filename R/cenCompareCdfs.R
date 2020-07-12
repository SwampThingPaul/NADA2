#' Comparison of empirical cdf
#'
#' @description Plots the empirical cdf and cdfs of three theoretical distributions, fit by maximum likelihood estimation (MLE).
#' @param y.var The column of y (response variable) values plus detection limits
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param dist3 Name of the third distribution to be plotted, default is `normal` distrubtion. Other distributions include `lnorm`(for log-normal), `gamma` and `weibull`(for Weibull).
#' @param Yname Optional – input text in quotes to be used as the variable name.  The default is the name of the `y.var` input variable.
#' @keywords ANOVA
#' @export
#' @return prints a plot of the empirial CDFs with BIC value for each distribution.
#' @importFrom fitdistrplus cdfcompcens fitdistcens
#'
#' @examples
#'
#' library(NADA) #For example data
#'
#' data(TCEReg)
#' cenCompareCdfs(TCEReg$TCEConc, TCEReg$TCECen)
#'
#' # With log-normal distribution
#' cenCompareCdfs(TCEReg$TCEConc, TCEReg$TCECen,dist3="lnorm")
#'
#' # Using an distribution not supported by this function (yet)
#' cenCompareCdfs(TCEReg$TCEConc, TCEReg$TCECen,dist3="beta")
#'
#' # With Yname specified
#' cenCompareCdfs(TCEReg$TCEConc, TCEReg$TCECen,Yname="TCE Conc (ug/L)\nLong Island, NY USA")

cenCompareCdfs <- function(y.var, cen.var, dist3="normal", Yname = yname)  {
  #added to stop if dist3 is not from the list
  if(!(dist3%in%c("normal","lnorm","gamma","weibull"))){stop(paste0(dist3," distribution is not supported with this function, try again."))}

  dist.vals <- c("normal","lnorm","gamma","weibull")
  dist.vals.text <- c("Normal","Lognormal","Gamma","Weibull")

  yname <- deparse(substitute(y.var))
  left <- y.var*(1-as.integer(cen.var))
  right <- y.var
  var.frame <- data.frame(left, right)

  y.dist1 <- fitdistcens(var.frame, "lnorm")
  y.dist2 <- fitdistcens(var.frame, "gamma")
  y.dist3 <- fitdistcens(var.frame, dist3)

  bic.dist1 <- paste("Lognormal BIC =", signif(y.dist1$bic, 3) )
  bic.dist2 <- paste("Gamma BIC =", signif(y.dist2$bic, 3) )
  bic.dist3 <- paste(dist.vals.text[match(dist3,dist.vals)],"BIC =", signif(y.dist3$bic, 3) )

  cdfcompcens(list(y.dist1, y.dist2, y.dist3), legendtext=c(bic.dist1, bic.dist2, bic.dist3), xlab = Yname, fitlty = c(1, 5, 3), lwd = c(1, 1, 2))

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
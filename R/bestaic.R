#' Find the lowest AIC multiple regression model
#'
#' @description Computes (2^k-1) censored regression models and their AIC statistics.  Prints out the lowest AIC models and the terms used.
#' @param y.var The column of y (response variable) values plus detection limits.
#' @param cen.var The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value is in `y.var`.
#' @param x.vars One or more uncensored explanatory variable(s). See Details
#' @param LOG Indicator of whether to compute the regression in the original y units, or on their logarithms.  The default is to use the logarithms (`LOG = TRUE`).  To compute in original units, specify the option `LOG = FALSE` (or `LOG = 0`).
#' @param n.models The number of models with their AIC values to be printed in the console window.  All (2^k-1) models are computed internally. This sets how many "best" (lowest AIC) models have output printed to the console.
#' @export
#' @return
#' Prints number of `x.vars`, lists `x.vars` and AIC values.
#'
#' @importFrom survival survreg Surv
#'
#' @details
#'
#' `x.vars`: If 1 x variable only, enter its name.  If multiple x variables, enter the name of a data frame of columns of the x variables. No extra columns unused in the regression allowed. Create this by `x.frame <- data.frame (Temp, Flow, Time)` for 3 variables (temperature, flow and time).
#'
#' AIC of each model is printed from lowest to highest AIC to help evaluate the ‘best’ regression model. n.models determines how many lines of model info is printed.
#'
#' LOG: The default is that the Y variable will be log transformed (LOG = TRUE).
#'
#' @seealso [survival::survreg]
#'
#' @references
#' Helsel, D.R., 2011. Statistics for censored environmental data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' @examples
#'\dontrun{
#' data(Brumbaugh)
#'
#' # Multiple regression
#' bestaic(Brumbaugh$Hg, Brumbaugh$HgCen, Brumbaugh[, c("SedMeHg","PctWetland", "SedAVS")])
#' }

bestaic <- function (y.var, cen.var, x.vars, LOG = TRUE, n.models = 10)
  {
  yname <- deparse(substitute(y.var))
  nonas <- na.omit(cbind(y.var, cen.var, x.vars))
  xnona <- data.frame(nonas[,-(1:2)])
  xnames <- names(xnona)
  nxvars <- ncol(xnona)      # number of x variables
  num.model = 0
  nmodels = (2^nxvars)-1

for (i in 1:nxvars)  {
  combs <- combn(nxvars, i)  # i is number of xvars in regression model
  ncolcombs <- ncol(combs)   # number of combinations for that model size
  for (j in 1:ncolcombs)
  { var.set <- t(combs)[j,]  # i elements in var.set
  frame <- xnona[, var.set[1:i]]   # for only 1 x variable
  if (i != 1) {frame <- data.frame(xnona[, var.set[1:i]])
  colnames(frame) <- xnames[var.set[1:i]] }
  reg.model <- cencorreg(nonas[,1], nonas[,2], frame, LOG = LOG, verbose = 0)
  names.xvars <- xnames[var.set[1:i]]
  aic <- -2*reg.model$loglik[2] + (2*reg.model$df +1)
  n.xvars <- length(names.xvars)
  n.xvars <- paste(n.xvars, "   ")

  for(m in 1:i) { x.name <- names.xvars[m]
  if (m == 1) {model.xvars <- x.name}
  #  else {name.temp <- paste(", ", x.name, sep = "")
  else {name.temp <- paste(" ", x.name, sep = "")
  model.xvars <- paste(model.xvars, name.temp, sep="")    }
  }  # end of getting model xvariables named

  num.model = num.model+1
  if (num.model == 1) {model.out <- data.frame (n.xvars, model.xvars, aic)}
  else {model.next <- data.frame (n.xvars, model.xvars, aic)
  model.out <- rbind(model.out, model.next) }
  }    # end of computing all models
}
o <- order(model.out$aic)
models.order <- model.out[o,]
models.best <- models.order[(1:n.models),]

  cat("Evaluating", nmodels, "models and printing the", n.models, "lowest AIC models", "\n")
  print(models.best, row.names = FALSE)

}

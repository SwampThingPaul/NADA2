## From NADA package (./All.R), code revised for S3 standards

#' Default Probability Values for Quantile Estimation
#'
#' A numeric vector of default probabilities used in quantile calculations,
#' commonly used in censored data methods.
#'
#' @format Numeric vector
#' @keywords data
#' @export
#' @examples
#' \dontrun{
#' NADAprobs
#' quantile(1:10, probs = NADAprobs)
#' }
NADAprobs <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)


#' @title NADA List Class
#'
#' @description A simple S3 class used to group related NADA objects as a list.
#'
#' @param ... Objects to include in the list.
#'
#' @return An object of class `NADAList`, which is a list of the input components.
#'
#' @export
NADAList <- function(...) {
  structure(list(...), class = "NADAList")
}


#' @title Show Method for NADAList Objects
#' @description Prints named elements of a `NADAList` object.
#' @param x A `NADAList` object.
#' @param ...	further arguments passed to or from other methods. (not used)
#' @export
print.NADAList <- function(x, ...) {
  tag <- names(x)
  for (i in seq_along(x)) {
    cat(tag[i], "\n")
    print(x[[i]])
    cat("\n")
  }
  invisible(x)
}


#' Internal: Get Plotting Limits in Log or Linear Scale
#'
#' @param log Character specifying log scale: "x", "y", or "xy".
#' @return Numeric vector of plot limits adjusted for log scale.
#' @keywords internal
cenpar.usr <- function(log) {
  usr <- par("usr")
  switch(log,
         xy = {usr <- 10^usr},
         x  = {usr[1:2] <- 10^usr[1:2]},
         y  = {usr[3:4] <- 10^usr[3:4]}
  )
  usr
}


#' Internal: Split Qualified Numeric Vector Into Observations and Censoring
#'
#' @param v Character vector with qualifying symbols like "<0.5".
#' @param qual.symbol Character symbol used to indicate censoring (default "<").
#' @return A list with elements `obs` (numeric) and `cen` (logical).
#' @keywords internal
splitQual <- function(v, qual.symbol = "<") {
  v <- as.character(v)
  obs <- as.numeric(sub(qual.symbol, "", v))
  cen <- rep(FALSE, length(obs))
  cen[grep(qual.symbol, v)] <- TRUE
  list(obs = obs, cen = cen)
}


#' Internal: Calculate Percent Censored Observations
#'
#' @param obs Numeric observations.
#' @param censored Logical vector indicating censored values.
#' @param na.action Function or character string specifying NA handling (default getOption("na.action")).
#' @return Percent (0-100) of censored observations.
#' @keywords internal
pctCen <- function(obs, censored, na.action = getOption("na.action")) {
  if (!is.logical(censored)) stop("censored must be logical vector!\n")
  if (is.null(na.action)) na.action <- "na.omit"
  obs <- do.call(na.action, list(obs))
  censored <- do.call(na.action, list(censored))
  100 * (length(obs[censored]) / length(obs))
}


#' Censored data sample size
#'
#' @description Computes the equivalent sample size of censored data.  Observations at lower detection limits have a greater percent of the equivalent information of a detected value than observations at higher detection limits.
#' @param y.var The column of data values plus detection limits.
#' @param y.cen The column of indicators, where 1 (or `TRUE`) indicates a detection limit in the `y.var` column, and 0 (or `FALSE`) indicates a detected value in `y.var`.
#' @param printstat Logical `TRUE`/`FALSE` option of whether to print the resulting statistics in the console window, or not.  Default is `TRUE.`
#' @keywords Sample Size censored
#' @export
#' @details
#' Based on "Method 2" of Dr. Brenda Gillespie's talk at ASA National Meeting 2019.  This method differs from hers in how the percentile probabilities for the detection limits are computed.  Probabilities here are computed using Regression on Order Statistics (ROS).
#'
#' Computes the equivalent n, the number of observations including censored values, as a measure of information content for data with nondetects.
#'
#' @return Prints summary statistics including
#' \itemize{
#' \item `n` sample size
#' \item `n.cen` number of censored data
#' \item `pct.cen` percent of data censored
#' \item `min` minimum reported value
#' \item `max` maximum reported value
#' }
#'
#' Summary of censored data including
#' \itemize{
#' \item `limit` detection limit
#' \item `n` number of censored values per limit
#' \item `uncen` number of detected values at or above the limit
#' \item `pexceed` proportion of data that exceeds the limit
#' }
#'
#' Summary of the equivalent sample size for detected and censored values.
#' \itemize{
#' \item `n.equiv` the equivalent number of observations
#' \item `n.cen.equiv` equivalent number of detected obs in the censored data
#' \item `n.detected` number of uncensored values
#'}
#' @references
#' Helsel, D.R., 2011. Statistics for Censored Environmental Data using Minitab and R, 2nd ed. John Wiley & Sons, USA, N.J.
#'
#' Gillespie, B.W., Dominguez, A., Li, Y., 2019. Quantifying the information in values below the detection limit (left-censored data).  Presented at the 2019 Joint Statistical Meetings of the Amer. Stat. Assoc., Denver, CO., July 31, 2019.
#' @examples
#' data(Brumbaugh)
#'
#' equivalent_n(Brumbaugh$Hg,Brumbaugh$HgCen)

equivalent_n <- function(y.var, y.cen,printstat=TRUE){
  # --- STOP IF NO CENSORED DATA ---
  if (!any(y.cen)) {
    stop("equivalent_n() requires at least one censored observation. No censored data detected.")
  }

  ycen <- as.logical(y.cen)
  yname <- deparse(substitute(y.var))

  # Summarize censored data
  aa <- censummary(y.var, ycen)

  # Prepare summary tables
  all_df <-data.frame(aa[1])
  all_df$variable <-rownames(all_df)
  rownames(all_df)<- NULL
  colnames(all_df[c(2,1)])<-c("variable","value")
  limits_df<-data.frame(aa[2])
  colnames(limits_df) <- c("limit","n","uncen","pexceed")

  aa <- list(all=all_df,limits=limits_df)

  # Equivalent sample size calculations
  n.equiv <- sum(aa$limits[,2]*aa$limits[,4] + aa$limits[,3])
  n_total <- sum(aa$limits[,2] + aa$limits[,3])
  n_detected <- sum(aa$limits[,3])
  n_cens <- sum(aa$limits[,2])
  n.cen.equiv <- sum(aa$limits[,2]*aa$limits[,4])
  n.detected <- n.equiv - n.cen.equiv
  equiv.out <- data.frame(n.equiv,  n.cen.equiv, n.detected)

  # Print results if requested
  if(printstat==TRUE){
  cat("data:",yname,"\n")
  print(aa)
  cat("equivalent sample size:", "\n")
  print(equiv.out, row.names = FALSE, print.gap = 3)
  }

  aa[["equivalent"]] <- equiv.out
  return (invisible (aa))
}

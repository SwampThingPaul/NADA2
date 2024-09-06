#' Thiamethoxam concentrations in pollen
#'
#' @description
#' Thiamethoxam concentrations in pollen from the Ontario Pollen Monitoring Network study https://data.ontario.ca/en/dataset/pollen-monitoring-network-study
#'
#' Variables are:
#' \describe{
#' \item{Thiamethoxam: }{Thiamethoxam concentration in Concentrations in microgram per gram.}
#' \item{ThiaCens: }{Censoring indicator.  1 denotes that the value in column 1 is a reporting limit not a specific concentration.}
#' \item{SamplingEvent: }{A grouping variable from the sample design.  A concentration is from 1 of 4 events in time.}
#' \item{ThiaAbvBelow: }{A binary variable denoting whether the Thiamethoxam concentration is above or below 0.05 ug/g.}
#' }
#'
#' @usage data(Pollen_Thia)
#' @aliases Pollen_Thia
#'
#' @docType data
#' @keywords dataset
#' @name Pollen_Thia
#' @source Ontario Ministry of Agriculture, Food and Rural Affairs (Pollen Monitoring Network)

"Pollen_Thia"

#' WildfireRisk: Spatial Prediction of Fire Risk
#'
#' This package implements the wildfire risk mapping method of Price et al.
#' (2015). Given a set of locations for which to assess risk; a model relating
#' fire risk to landscape variables; and a set of raster layers for those
#' variables; the package provides functions to scan the neighbourhood of each
#' location and apply the given model to estimate the point risk.
#'
#' @section Things to know before using:
#' Presently the package can only work with models similar to that of Price et
#' al. (2015).
#'
#' @references
#' Price O, Borah R, Bradstock R, Penman T (2015) An empirical wildfire risk
#' analysis: the probability of a fire spreading to the urban interface in
#' Sydney, Australia. International Journal of Wildland Fire 24, 597.
#' doi:10.1071/WF14160.
#' \href{http://www.publish.csiro.au/?paper=WF14160}{Article at journal website}
#'
#' @docType package
#' @name WildfireRisk
#'
NULL

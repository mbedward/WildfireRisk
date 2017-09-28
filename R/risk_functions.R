#' Calculate risk of wildfire at given locations
#'
#' This function takes a data frame of scan lines created with \code{make_scan_lines}
#' and samples time since fire and forest cover values along each line. It then uses
#' these values to calculate the probability of wildfire travelling along each line.
#' The calculation is based on the regression equation of Price \emph{et al.} 2015.
#' The function then calculates a second probability, for comparison, with the time
#' since fire value for each line set to \code{tsf.max} (a high value representing
#' long unburnt).
#'
#' If the input scan lines have a coordinate reference system, and argument \code{as.sf}
#' is \code{TRUE} (default), then the same CRS will be set for the output. However,
#' note that this function does not check that the scan lines and raster layers are in
#' compatible projections, or that the value for sample.spacing is sensible for
#' the map units being used. This is up to you!
#'
#'
#' @param lines An \code{sf} object (spatial data frame) with scan lines
#'   for each location at which to predict fire risk.
#'
#' @param tsf A raster of time since fire values.
#'
#' @param forest A raster of forest cover where 0 indicates
#'   non-forest cells and all other values indicate forest cells.
#'
#' @param sample.spacing The spacing between sample points along scan lines.
#'
#' @param ffdi Constant FFDI value to use for calculations (default = 50).
#'
#' @param kbdi Constant KBDI value to use for calculations (default = 100).
#'
#' @param tsfmax High value of time since fire to use for calculation of risk
#'   if all scan lines are long unburnt (default = 50).
#'
#' @param as.sf If \code{TRUE} (default), the returned data frame will be an
#'   \code{sf} object where the geometry column has the input scan lines.
#'   If \code{FALSE}, the returned object is a simple data frame with no
#'   geometry column.
#'
#' @return A data frame giving, for each location and scan line, the
#'   variables used for the risk calculation and the predicted probabilities
#'   for both observed and long-unburnt time since fire (columns 'pobs' and 'pmax').
#'
#' @importFrom dplyr %>% do group_by mutate right_join summarize ungroup
#'
#' @export
#'
calculate_risk <- function(lines,
                           tsf, forest,
                           sample.spacing,
                           ffdi = 50, kbdi = 100, tsfmax = 50,
                           as.sf = TRUE) {

  forest <- raster::calc(forest, fun = function(x) x != 0)

  layers <- list(tsf = tsf, forest = forest)

  # If the scan lines have a map projection, remove it temporarily
  # so that we don't get pesky errors about incompatible units
  # in the calculations here.
  crs.in <- st_crs(lines)
  st_crs(lines) <- NA_crs_


  # Risk function
  riskfn <- function(tsf_mean, forest_p, distance, is_west) {
    plogis(58.658 +
             -9.118 * log(distance) +
             -40.2 * forest_p +
             0.041 * kbdi +
             0.117 * tsf_mean +
             0.0931 * ffdi +
             6.4 * log(distance) * forest_p +
             2.253 * is_west)
  }


  dat <- lines %>%

    # sample rasters and calculate summary statistics for the
    # scan lines at each location
    group_by(locationid) %>%

    do({
      samples <- sample_raster(layers, ., sample.spacing) %>%

        as.data.frame() %>%

        group_by(lineid) %>%

        summarize(tsf_mean = mean(tsf, na.rm = TRUE),
                  forest_p = mean(forest, na.rm = TRUE))
    }) %>%

    # turn back into an 'sf' object (perhaps there is some
    # way of avoiding having to do this?)
    ungroup() %>%
    right_join(lines) %>%
    sf::st_sf() %>%

    # calculate remaining line variables and probability
    # of fire travel for observed and maximum time since fire values.
    mutate(distance = sf::st_length(geometry),
           is_west = as.integer( is_line_west(geometry) ),
           pobs = riskfn(tsf_mean, forest_p, distance, is_west),
           pmax = riskfn(tsfmax, forest_p, distance, is_west) )


  if (!as.sf) {
    dat <- sf::st_sf(dat)
    dat$geometry <- NULL
  } else {
    st_crs(dat) <- crs.in
  }

  dat
}


#' Summarize scan line risk values
#'
#' This function takes an \code{sf} data frame of risk values for scan lines,
#' as returned by \code{calculate_risk}, and summarizes the probability values
#' over the set of lines for each location by: mean; median; lower quartile and
#' upper quartile.
#'
#' @param risk An \code{sf} data frame of risk values for scan lines as
#'   returned by \code{calculate_risk}.
#'
#' @return An \code{sf} data frame of summary risk statistics for each location.
#'
#' @seealso \code{\link{calculate_risk}}
#'
#' @importFrom dplyr do group_by
#'
#' @export
#'
summarize_risk <- function(risk) {
  # Helper function to retrieve central point from a
  # set of scan lines
  firstpoint <- function(lines) {
    m <- st_coordinates(lines)
    data.frame(x = m[1,1], y = m[1,2])
  }

  # Helper function to calculate quantile
  q <- function(x, prob) quantile(x, prob, na.rm = TRUE)

  # Get point locations
  loc <- risk %>%
    group_by(locationid) %>%
    do(firstpoint(.$geometry))

  # Summary statistics for each location
  pstats <- risk %>%
    # drop scan lines
    as.data.frame() %>%

    # calculate mean probabilities
    group_by(locationid) %>%

    summarize(pobs_mean = mean(pobs, na.rm = TRUE),
              pobs_25 = q(pobs, 0.25),
              pobs_50 = q(pobs, 0.5),
              pobs_75 = q(pobs, 0.75),

              pmax_mean = mean(pmax, na.rm = TRUE),
              pmax_25 = q(pmax, 0.25),
              pmax_50 = q(pmax, 0.5),
              pmax_75 = q(pmax, 0.75) ) %>%

    # join location data
    left_join(loc, by = "locationid") %>%

    # convert to a spatial (sf) object with point geometry
    st_as_sf(coords = c("x", "y"))


  # Set coordinate reference system via the EPSG code for Zone 55
  st_crs(pstats) <- st_crs(risk)

  pstats
}

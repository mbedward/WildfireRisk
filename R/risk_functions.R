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
#' If the input scan lines have a coordinate reference system defined, it will be set
#' for the output.
#'
#' @note This function does not check that the scan lines and raster layers are in
#' compatible projections, or that the value for sample.spacing is sensible for the
#' map units being used.
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
                           ffdi = 50, kbdi = 100, tsfmax = 50) {

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
      sample_raster(layers, ., sample.spacing) %>%

        as.data.frame() %>%

        group_by(lineid) %>%

        summarize(tsf_mean = mean(tsf, na.rm = TRUE),
                  forest_p = mean(forest, na.rm = TRUE))
    }) %>%


  dat <- dat %>%

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

  st_crs(dat) <- crs.in

  dat
}


#' Summarize scan line risk values
#'
#' This function takes an \code{sf} data frame of risk values for scan lines,
#' as returned by \code{calculate_risk}, and summarizes the probability values
#' over the set of lines for each location by mean and (optionally) quantiles.
#'
#' Scan line risk values are probabilities in the strict sense. In constrast,
#' the summary statistics for locations derived by this function, such as the
#' mean of scan line probabilities are not probabilities despite being
#' restricted to the unit inteval.
#'
#' @param risk An \code{sf} data frame of risk values for scan lines as
#'   returned by \code{calculate_risk}.
#'
#' @param quantiles A vector of quantiles (probabilities with values in [0,1])
#'   to calculate in addition to the mean. If \code{NULL} or an empty vector
#'   no quantiles are calculated.
#'
#' @return An \code{sf} data frame of summary risk statistics for each location.
#'
#' @seealso \code{\link{calculate_risk}}
#'
#' @importFrom dplyr %>% do group_by
#'
#' @export
#'
#' @examples
#' \dontrun{
#' risk.lines <- calculate_risk(scanlines, tsf, forest, 100)
#'
#' risk.locations <- summarize_risk(risk.lines, quantiles = c(0.1, 0.9))
#' }
#'
summarize_risk <- function(risk, quantiles = c(0.25, 0.75)) {
  # Helper function to retrieve central point from a
  # set of scan lines
  firstpoint <- function(lines) {
    m <- st_coordinates(lines)
    data.frame(x = m[1,1], y = m[1,2])
  }

  has.quantiles <- !is.null(quantiles) & length(quantiles) > 0

  if (has.quantiles) {
    qnames <- names(quantile(1, quantiles)) %>% stringr::str_replace("\\%", "")
  }

  # Get point locations
  loc <- risk %>%
    group_by(locationid) %>%
    do(firstpoint(.$geometry))


  # Helper function to calculate mean and quantiles and
  # return them as a data frame
  fn <- function(x, varname) {
    d <- data.frame(mu = mean(x, na.rm = TRUE))
    colnames(d) <- paste(varname, "mean", sep = "_")

    if (has.quantiles) {
      q <- quantile(x, probs = quantiles, na.rm = TRUE)
      q <- t(q)
      colnames(q) <- paste(varname, qnames, sep = "_")

      d <- cbind(d, q)
    }

    d
  }


  # Summary statistics for each location
  pstats <- risk %>%
    # drop scan lines
    as.data.frame() %>%

    # calculate mean probabilities
    group_by(locationid) %>%

    do({
      dobs <- fn(.$pobs, "pobs")
      dmax <- fn(.$pmax, "pmax")
      cbind(dobs, dmax)
    }) %>%

    ungroup() %>%

    # join location data
    left_join(loc, by = "locationid") %>%

    # convert to a spatial (sf) object with point geometry
    st_as_sf(coords = c("x", "y"))


  # Set coordinate reference system via the EPSG code for Zone 55
  st_crs(pstats) <- st_crs(risk)

  pstats
}


#' Profile polygons by scan line risk values
#'
#' @param blocks Polygons (or multi-polygons) to profile. These can be provided
#'   as an \code{sf} data frame or as the path to a file of vector data
#'   (e.g. ESRI shapefile).
#'
#' @param risk An \code{sf} data frame of risk values for scan lines as
#'   returned by \code{calculate_risk}.
#'
#' @param intersections An optional list of intersection data returned by a
#'   previous call to this function.
#'
#' @param strict.crs If TRUE, the blocks and risk objects must have
#'   exactly the same coordinate reference system defined. If FALSE,
#'   the function will make relaxed assumptions about the reference
#'   systems as described for function \code{lines_in_blocks}.
#'
#' @param quiet If TRUE, suppress messages and warnings about coordinate
#'   reference systems of objects.
#'
#' @return A list with two elements: blocks and intersections. The blocks
#'   element is an \code{sf} object with data from the input \code{risk} object
#'   plus columns giving the number of scan lines intersecting each block;
#'   number of locations to which the intersected scan lines belong; and
#'   the mean value of scan line probabilities. The intersections element is
#'   a list of vectors, one for each block, where vector values are indices
#'   of intersected scan lines. This can be used with the function on
#'   subsequent calls to avoid re-doing the intersection between blocks and
#'   lines.
#'
#' @importFrom dplyr %>% do group_by mutate row_number select ungroup
#'
#' @seealso \code{\link{lines_in_blocks}}
#'
#' @export
#'
profile_blocks <- function(blocks, risk,
                           intersections = NULL,
                           strict.crs = FALSE,
                           quiet = FALSE) {

  blocks <- .get_sf_object(blocks)

  if (is.null(intersections))
    intersections <- lines_in_blocks(blocks, risk, strict.crs, quiet)
  else if (!is.list(intersections) | length(intersections) != nrow(blocks))
    stop("Argument intersections should be a list of vectors with the\n",
         "  number of list elements equal to the number of blocks")

  res <- blocks %>%
    mutate(i__ = row_number()) %>%

    group_by(i__) %>%

    do({
      ii <- intersections[[.$i__]]
      x <- as.data.frame(risk[ii, ])

      if (length(ii) == 0)
        data.frame(nlines = 0, nlocations = 0, pobs_mean = NA)
      else
        data.frame(nlines = length(ii),
                   nlocations = n_distinct(x$locationid),
                   pobs_mean = mean(x$pobs, na.rm = TRUE) )

    }) %>%

    ungroup() %>%
    dplyr::select(-i__)

  list(blocks = st_sf(data.frame(blocks, res)),
       intersections = intersections)
}


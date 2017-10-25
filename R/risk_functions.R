#' Calculate a scan line risk value
#'
#' This function applies the equation of Price \emph{et al.} 2015 to
#' calculate risk (probability of fire travel along scan line) given
#' line values for mean time since fire, proportion of forest cover,
#' distance (line length) and whether the line has a westerly orientation.
#'
#' Generally, you won't call this function directly. It used by the
#' \code{calculate_risk} and \code{treat_blocks} functions.
#'
#' @param tsf_mean Mean value of time since fire for points sampled along
#'   the scan line.
#'
#' @param forest_p Proportion of points sampled along the scan line that
#'   were in forest.
#'
#' @param distance Line length in metres.
#'
#' @param is_west \code{TRUE} if the line has a westerly orientation;
#'   \code{FALSE} otherwise.
#'
#' @param ffdi Constant FFDI value to use for calculations (default = 50).
#'
#' @param kbdi Constant KBDI value to use for calculations (default = 100).
#'
#' @seealso \code{\link{calculate_risk}} \code{\link{treat_blocks}}
#'
#' @export
#
line_risk <- function(tsf_mean, forest_p, distance, is_west,
                      ffdi = 50, kbdi = 100) {
  plogis(58.658 +
           -9.118 * log(distance) +
           -40.2 * forest_p +
           0.041 * kbdi +
           0.117 * tsf_mean +
           0.0931 * ffdi +
           6.4 * log(distance) * forest_p +
           2.253 * is_west)
}




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
#' The coordinates of points along each scan line and corresponding values of
#' time since fire and forest cover are cached in a SQLite database for later
#' use. The path to this file is recorded as the \code{"dbname"} attribute of the
#' returned \code{risk} object. If you rename or move the database (e.g. when
#' transferring files to another system) you should manually update the attribute
#' with the new full path to the database. For example:
#' \code{attr(risk, "dbname") <- "c:/foo/bar/pointdata.db"}.
#'
#' @note This function assumes that scan lines and raster layers are in
#' compatible projections with metres as the distance unit. If the input scan lines
#' have a coordinate reference system defined, it will be set for the output.
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
#' @param point.db The path for the SQLite database file in which to cache
#'   sample point data for scan lines. If the specified file already exists
#'   it will be over-written and a warning message issued.
#'
#' @param ffdi Constant FFDI value to use for calculations (default = 50).
#'
#' @param kbdi Constant KBDI value to use for calculations (default = 100).
#'
#' @param tsfmax High value of time since fire to use for calculation of risk
#'   if all scan lines are long unburnt (default = 50).
#'
#'
#' @return A \code{risk} object (a type of \code{sf} spatial data frame) containing,
#'   for each location and scan line, the variables used for the risk calculation
#'   and the predicted probabilities for both observed and long-unburnt time since
#'   fire (columns 'pobs' and 'pmax'). The path/filename of the associated SQLite
#'   database with point data for the scan lines is recorded as the \code{"dbname"}
#'   attribute.
#'
#' @importFrom dplyr %>% do group_by mutate right_join summarize ungroup
#'
#' @export
#'
calculate_risk <- function(lines,
                           tsf, forest,
                           sample.spacing,
                           dbname = "pointdata.db",
                           ffdi = 50, kbdi = 100, tsfmax = 50) {

  forest <- raster::calc(forest, fun = function(x) !is.na(x) & x != 0)

  layers <- list(tsf = tsf, forest = forest)


  # If the scan lines have a map projection, remove it temporarily
  # so that we don't get pesky errors about incompatible units
  # in the calculations here.
  crs.in <- st_crs(lines)
  st_crs(lines) <- NA_crs_


  # Database connection for point data
  if (file.exists(dbname)) {
    warning("Overwriting existing file: ", dbname, "\n")
    unlink(dbname)
  }

  con <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = dbname)

  # Sample rasters at points along each scan line

  first.write <- TRUE

  dat <- lines %>%
    group_by(locationid) %>%

    do({
      pdat <- sample_raster(layers, ., sample.spacing)
      if (first.write) {
        RSQLite::dbWriteTable(con, "pointdata", pdat)
        RSQLite::dbGetQuery(con, "CREATE INDEX index_location ON pointdata (locationid)")
        first.write <- FALSE
      } else {
        RSQLite::dbWriteTable(con, "pointdata", pdat, append = TRUE)
      }

      pdat %>%
        group_by(locationid, lineid) %>%
        summarize(tsf_mean = mean(tsf, na.rm = TRUE),
                  forest_p = mean(forest, na.rm = TRUE))
    })

  RSQLite::dbDisconnect(con)

  dat <- dat %>%
    # turn back into an 'sf' object (perhaps there is some
    # way of avoiding having to do this?)
    ungroup() %>%
    right_join(lines, by = c("locationid", "lineid")) %>%
    sf::st_sf() %>%

    # Calculate remaining line variables and probability
    # of fire travel for observed and maximum time since fire values.
    #
    # The 'as.numeric' call for distance is to discard the units
    # returned by st_length and avoid having to define units for
    # other quantities in the helper function 'riskfn'
    #
    mutate(distance = as.numeric( sf::st_length(geometry) ),
           is_west = as.integer( is_line_west(geometry) ),
           pobs = riskfn(tsf_mean, forest_p, distance, is_west),
           pmax = riskfn(tsfmax, forest_p, distance, is_west) )

  st_crs(dat) <- crs.in

  class(dat) <- c("risk", class(dat))
  attr(dat, "dbname") <- normalizePath(dbname)

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


# Private helper to do block x line intersections if required
#
.do_block_line_intersect <- function(blocks, lines, intersections,
                                     strict.crs = FALSE, quiet = TRUE) {
  if (is.null(intersections))
    intersections <- lines_in_blocks(blocks, risk, strict.crs, quiet)

  else if (!is.list(intersections) | length(intersections) != nrow(blocks))
    stop("Argument intersections should be a list of vectors with the\n",
         "  number of list elements equal to the number of blocks")

  intersections
}



#' Simulate prescribed burning of landscape blocks and update risk values
#'
#' This function takes one or more blocks (polygons representing management units)
#' and simulates prescribed burning by setting time since fire within
#' each block to zero. Data for scan lines passing through the blocks is
#' then updated and the risk value (probablity of fire traversing the line)
#' recalculated for each line. The line geometries and their prior risk
#' values are taken from the \code{risk} object which is an \code{sf}
#' data frame as returned by \code{calculate_risk}.
#'
#' Sample point locations and prior data for each scan line are retrieved
#' from a database file specified by the \code{"dbname"} attribute of
#' the input \code{risk} object. See \code{calculate_risk} for more details.
#'
#' @param risk An \code{sf} data frame of scan lines with associated
#'   baseline risk values as returned by function \code{calculate_risk}.
#'
#' @param blocks An \code{sf} object of polygons for which to simulate
#'   prescribed burning.
#'
#' @return A data frame of updated scan line risk values with columns:
#'   locationid, lineid, ptreat.
#'
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
#'
#' @examples
#' \dontrun{
#' # Load a polygon shapefile of prescribed burning blocks
#' blocks <- st_read("burnblocks.shp")
#'
#' # Simulate treatment for three selected blocks
#' risk.updated <- treat_blocks(blocks[c(1, 5, 42), ], risk)
#' }
#'
#' @export
#'
treat_blocks <- function(blocks, risk) {

  con <- dbConnect(SQLite(), dbname = attr(risk, "dbname"))

  intersections <- lines_in_blocks(blocks, risk, by = "line")

  crossing.lines <- which(sapply(intersections, length) > 0)

  # For each line passing through one or more blocks, identify sample points
  # within blocks, set time since fire of those points to zero, and
  # re-calculate line risk

  sql <- "select * from pointdata where locationid = :loc and lineid = :line"

  ptreat <- risk$pobs

  for (iline in crossing.lines) {
    blk.indices <- intersections[[iline]]

    if (length(blk.indices) > 0) {
      # Get points for this scan line
      params <- list(loc = risk$locationid[iline],
                     line = risk$lineid[iline])

      pdat <- dbGetQuery(con, sql, params)

      # Set time since fire of points within blocks
      # to zero
      pdat <- st_as_sf(pdat, coords = c("x", "y"), crs = st_crs(blocks))
      ii <- st_intersects(pdat, blocks[blk.indices, ])
      ii <- sapply(ii, length) > 0

      pdat$tsf[ii] <- 0
      tsf.treated <- mean(pdat$tsf, na.rm = TRUE)

      p_treated[iline] <- line_risk(
        tsf_mean = tsf.treated,
        forest_p = risk$forest_p[iline],
        distance = risk$distance[iline],
        is_west = risk$is_west[iline])
    }
  }

  dbDisconnect(con)

  data.frame(risk[, c("locationid", "lineid")], ptreat)
}


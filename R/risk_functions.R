#' Calculate a scan line risk value
#'
#' Applies the equation of Price \emph{et al.} 2015 to calculate risk
#' (probability of fire travel along scan line) given line values for mean time
#' since fire, proportion of forest cover, distance (line length) and whether
#' the line has a westerly orientation. This function primarily serves as a
#' helper for the \code{calculate_risk} and \code{treat_blocks} functions.
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
calculate_line_risk <- function(tsf_mean, forest_p, distance, is_west,
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
#' The coordinates and data for points sampled along each scan line are written
#' to a SQLite database into a table named \code{'pointdata'}. Data for scan
#' lines, including end-points, location identifiers, values used for the risk
#' calculation (e.g. mean time since fire) and risk values are also written to
#' the database into a table named \code{'linedata'}. The path to this file is
#' recorded as the \code{"dbname"} attribute of the returned \code{risk} object.
#' If you rename or move the database (e.g. when transferring files to another
#' system) you should manually update the attribute with the new full path to
#' the database.
#' For example: \code{attr(risk, "dbname") <- "c:/foo/bar/pointdata.db"}.
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
#'   fire (columns 'pobs' and 'pmax'). These data are also written to the SQLite
#'   database file. The path/filename of the database with point data is recorded
#'   as the \code{"dbname"} attribute of the returned object.
#'
#' @importFrom dplyr %>% do group_by mutate right_join summarize ungroup
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbWriteTable
#' @importFrom RSQLite SQLite
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


  # Database connection for cached point and line summary data
  if (file.exists(dbname)) {
    warning("Overwriting existing database: ", dbname, "\n")
    unlink(dbname)
  }

  con <- dbConnect(RSQLite::SQLite(), dbname = dbname)

  # Sample rasters at points along each scan line

  first.write <- TRUE

  dat <- lines %>%
    group_by(locationid) %>%

    do({
      pdat <- sample_raster(layers, ., sample.spacing)
      if (first.write) {
        dbWriteTable(con, "pointdata", pdat)
        dbExecute(con, "CREATE INDEX index_pointloc ON pointdata (locationid)")

        first.write <- FALSE
      } else {
        dbWriteTable(con, "pointdata", pdat, append = TRUE)
      }

      pdat %>%
        group_by(locationid, lineid) %>%
        summarize(tsf_mean = mean(tsf, na.rm = TRUE),
                  forest_p = mean(forest, na.rm = TRUE))
    })

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
    # other quantities in the helper function 'calculate_line_risk'
    #
    mutate(distance = as.numeric( sf::st_length(geometry) ),
           is_west = as.integer( is_line_west(geometry) ),
           pobs = calculate_line_risk(tsf_mean, forest_p, distance, is_west),
           pmax = calculate_line_risk(tsfmax, forest_p, distance, is_west) )

  st_crs(dat) <- crs.in

  class(dat) <- c("risk", class(dat))
  attr(dat, "dbname") <- normalizePath(dbname)

  # Write the line level data to the database. We can't write
  # the risk/sf object directly because RSQLite doesn't like
  # the geometry column, so we create a plain data base with
  # columns for line end-points
  endpoints <- dat %>%
    st_coordinates() %>%
    as.data.frame() %>%
    group_by(L1) %>%
    summarize(x0 = first(X), y0 = first(Y), x1 = last(X), y1 = last(Y))

  tbl <- dat %>%
    as.data.frame() %>%
    cbind(endpoints) %>%
    select(-geometry, -L1)

  dbWriteTable(con, "linedata", tbl)
  dbExecute(con, "CREATE INDEX index_lineloc ON linedata (locationid)")

  dbDisconnect(con)

  dat
}


#' Summarize scan line risk values by location
#'
#' This function takes a set of risk values for scan lines,
#' as returned by \code{calculate_risk}, and summarizes the probability values
#' over the set of lines for each location by mean and (optionally) specified
#' quantiles.
#'
#' @param line.risk Either a \code{risk} object as returned by \code{calculate_risk}
#'   or a character string giving the path and filename of a SQLite database
#'   of risk data.
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
#' risk.locations <- summarize_location_risk(risk.lines, quantiles = c(0.1, 0.9))
#' }
#'
summarize_location_risk <- function(line.risk, quantiles = c(0.25, 0.75)) {

  if (is.character(line.risk))
    line.risk <- load_line_risk_table(line.risk[1])
  else if (!inherits(line.risk, "risk"))
    stop("Argument line.risk should be a risk data frame or the path to a database file")

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
  loc <- line.risk %>%
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
  pstats <- line.risk %>%
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


  # Set coordinate reference system
  st_crs(pstats) <- st_crs(line.risk)

  pstats
}


#' Summarize scan line risk values for landscape blocks
#'
#' This function takes a set of risk values for scan lines,
#' as returned by \code{calculate_risk}, and summarizes the probability values
#' of lines intersecting with landscape blocks. Blocks are polygons representing
#' management units (e.g. areas for prescribed burning). For each block the
#' function calculates the mean probability of all intersecting lines and,
#' (optionally) specified quantiles.
#'
#' @param line.risk Either a \code{risk} data frame as returned by
#'   \code{calculate_risk} or a character string giving the path and filename
#'   of a SQLite database of risk data.
#'
#' @param blocks Polygons (or multi-polygons) to profile. These can be provided
#'   as an \code{sf} data frame or as the path to a file of vector data
#'   (e.g. ESRI shapefile).
#'
#' @param quantiles A vector of quantiles (probabilities with values in [0,1])
#'   to calculate in addition to the mean. If \code{NULL} or an empty vector
#'   no quantiles are calculated.
#'
#' @param strict.crs If TRUE, the blocks and risk objects must have
#'   exactly the same coordinate reference system defined. If FALSE,
#'   the function will make relaxed assumptions about the reference
#'   systems as described for function \code{lines_in_blocks}.
#'
#' @param quiet If TRUE, suppress messages and warnings about coordinate
#'   reference systems of objects.
#'
#' @return An \code{sf} object with data from the input \code{risk} object
#'   plus columns giving the number of scan lines intersecting each block;
#'   number of locations to which the intersected scan lines belong;
#'   the mean value of scan line probabilities, and (optionally) specified
#'   quantiles.
#'
#' @examples
#' \dontrun{
#' risk.lines <- calculate_risk(scanlines, tsf, forest, 100)
#' blocks <- st_read("c:/foo/bar/data/burning_blocks.shp")
#' risk.blocks <- summarize_block_risk(risk.lines, blocks, quantiles = c(0.1, 0.9))
#' }
#'
#' @importFrom dplyr %>% do group_by mutate row_number select ungroup
#'
#' @seealso \code{\link{lines_in_blocks}}
#'
#' @export
#'
summarize_block_risk <- function(line.risk, blocks, quantiles = c(0.25, 0.75),
                                 strict.crs = FALSE,
                                 quiet = FALSE) {

  if (is.character(line.risk))
    line.risk <- load_line_risk_table(line.risk[1])
  else if (!inherits(line.risk, "risk"))
    stop("Argument line.risk should be a risk data frame or the path to a database file")

  blocks <- .get_sf_object(blocks)

  has.quantiles <- !is.null(quantiles) & length(quantiles) > 0

  if (has.quantiles) {
    qnames <- names(quantile(1, quantiles)) %>% stringr::str_replace("\\%", "")
  }

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

  # Data frame for results when block has no scan lines
  NoData <- cbind(nlines = 0, nlocations = 0, fn(0, "pobs"), fn(0, "pmax"))
  NoData[, 3:ncol(NoData)] <- NA

  intersections <- lines_in_blocks(blocks, line.risk, by = "block",
                                   strict.crs, quiet)

  res <- blocks %>%
    mutate(i__ = row_number()) %>%

    group_by(i__) %>%

    do({
      ii <- intersections[[.$i__]]
      x <- as.data.frame(line.risk[ii, ])

      if (length(ii) == 0)
        NoData
      else
        data.frame(nlines = length(ii),
                   nlocations = n_distinct(line.risk$locationid[ii]),
                   fn(line.risk$pobs[ii], "pobs"),
                   fn(line.risk$pmax[ii], "pmax") )
    }) %>%

    ungroup() %>%
    dplyr::select(-i__)

  st_sf(data.frame(blocks, res))
}



#' Simulate prescribed burning of landscape blocks and update risk values
#'
#' This function takes one or more blocks (polygons representing management
#' units) and simulates prescribed burning by setting the time since fire within
#' each block in turn to zero. Data for scan lines passing through the block
#' being considered are then updated and risk value (probablity of fire
#' traversing the line) recalculated for each line. The line geometries and
#' their prior risk values are taken from the \code{risk} object which is an
#' \code{sf} data frame as returned by \code{calculate_risk}.
#'
#' Sample point locations and prior data for each scan line are retrieved
#' from a database file specified by the \code{"dbname"} attribute of
#' the input \code{risk} object. See \code{calculate_risk} for more details.
#'
#' @note At present this function considers each block singly. Later, an option
#'   could be added to consider blocks in groups.
#'
#' @param line.risk Either a \code{risk} object as returned by \code{calculate_risk}
#'   or a character string giving the path and filename of a SQLite database
#'   of risk data.
#'
#' @param block.risk An \code{sf} object of baseline risk values for blocks as
#'   returned by \code{summarize_block_risk}.
#'
#' @return A data frame of updated scan line risk values with columns:
#'   locationid, lineid, ptreat.
#'
#' @importFrom dplyr %>% group_by summarize
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
#'
#' @examples
#' \dontrun{
#' # Load a polygon shapefile of prescribed burning blocks
#' blocks <- st_read("burnblocks.shp")
#'
#' # Calculate baseline risk values for blocks
#' block.risk <- summarize_block_risk(line.risk, blocks)
#'
#' # Simulate burning each block in turn and report updated
#' # block statistics
#' block.risk.treated <- treat_blocks(block.risk, line.risk)
#' }
#'
#' @export
#'
treat_blocks <- function(block.risk, line.risk) {

  # Ensure the block.risk object has a pobs_mean column
  if (!("pobs_mean" %in% colnames(block.risk)) )
    stop("Argument block.risk should have a column pobs_mean")

  if (is.character(line.risk))
    line.risk <- load_line_risk_table(line.risk[1])
  else if (!inherits(line.risk, "risk"))
    stop("Argument line.risk should be a risk data frame or the path to a database file")

  # FIXME
  line.risk$locationid <- as.character(line.risk$locationid)

  con <- dbConnect(SQLite(), dbname = attr(line.risk, "dbname"))

  intersections <- lines_in_blocks(block.risk, line.risk, by = "block")

  blocks.with.lines <- which(sapply(intersections, length) > 0)

  if (length(blocks.with.lines) == 0) {
    warning("No intersecting scan lines found for any block.\n")
    block.risk$ptreat_mean <- block.risk$pobs_mean
    return(block.risk)
  }


  # For each block intersected by scan lines, identify sample points
  # within the block, set time since fire of those points to zero,
  # re-calculate line risk values, and summarize for the block.

  sql <- "SELECT * FROM pointdata WHERE locationid = :loc and lineid = :line"

  block.risk$ptreat_mean <- block.risk$pobs_mean

  k <- 0
  pb <- txtProgressBar(0, length(blocks.with.lines), style = 3)
  for (iblock in blocks.with.lines) {
    ii <- intersections[[iblock]]
    nlines <- length(ii)

    # Get point data for these scan lines
    pdat <- lapply(ii, function(i) {
      params <- list(loc = line.risk$locationid[i],
                     line = line.risk$lineid[i])

      dbGetQuery(con, sql, params = params)
    })

    pdat <- do.call(rbind, pdat)


    # Set time since fire of points within the block
    # to zero
    pdat <- st_as_sf(pdat, coords = c("x", "y"), crs = st_crs(block.risk))
    ii <- st_intersects(pdat, block.risk[iblock, ])
    ii <- sapply(ii, length) > 0

    pdat$tsf[ii] <- 0

    # Calculate updated line risk values
    ldat <- pdat %>%
      group_by(locationid, lineid) %>%
      summarize(tsf_treat_mean = mean(tsf, na.rm = TRUE)) %>%
      ungroup() %>%

      left_join(as.data.frame(line.risk), by = c("locationid", "lineid")) %>%

      mutate(ptreat = calculate_line_risk(tsf_mean = tsf_treat_mean,
                                          forest_p = forest_p,
                                          distance = distance,
                                          is_west = is_west))

    block.risk$ptreat_mean[iblock] <- mean(ldat$ptreat, na.rm = TRUE)

    k <- k + 1
    setTxtProgressBar(pb, k)
  }
  close(pb)

  dbDisconnect(con)

  block.risk
}


#' Load the table of line data and risk values from a database
#'
#' Reads the table of scan line data and risk values from a SQLite database
#' (as created by function \code{calculate_risk}) and returns it as a spatial
#' data frame (a type of \code{sf} object of class \code{risk}).
#'
#' @importFrom dplyr %>% select
#' @importFrom DBI dbConnect dbDisconnect dbReadTable
#' @importFrom RSQLite SQLite
#'
#' @export
#'
load_line_risk_table <- function(path) {
  con <- dbConnect(SQLite(), path)
  dat <- dbReadTable(con, "linedata")
  dbDisconnect(con)

  # Re-construct line geometries from end-point coordinates
  m <- as.matrix( dat[, c("x0", "y0", "x1", "y1")] )
  lines <- lapply(1:nrow(m), function(i) st_linestring(matrix(m[i,], ncol = 2, byrow = TRUE)))

  # Turn into an 'sf' object with a line geometry column
  dat <- dat %>%
    select(-x0, -y0, -x1, -y1) %>%
    st_sf(., geometry = st_sfc(lines))

  # Add 'risk' class stuff
  class(dat) <- c("risk", class(dat))
  attr(dat, "dbname") <- normalizePath(path)

  dat
}


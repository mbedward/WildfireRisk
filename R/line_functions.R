#' Generate radial scan lines from one or more locations
#'
#' This function generates a set of \code{n} radial scan lines emanating from
#' each of a given set of locations. Line lengths can be constant; individually
#' specified; or set by a user-provided function (e.g. to draw them from a
#' particular distribution). Line angles can be regularly spaced (the default);
#' individually specified; or set by a user-provided function.
#'
#' @param locations A matrix or data frame with either two columns for X-Y
#'   coordinates; or three columns for location identifier (alpha-numeric)
#'   then X-Y coordinates in that order. If location identifiers are not
#'   provided, these will be set as consecutive integers in the returned
#'   scan lines.
#'
#' @param nlines Number of lines to generate for each point.
#'
#' @param lengths Either a single value for uniform line length; or a
#'   vector of \code{n} line lengths; or a function that returns a vector
#'   of lengths.
#'
#' @param angles Either a vector of line directions (degrees); or a
#'   function that returns a vector of directions; or \code{NULL}
#'   for uniform distribution of line angles.
#'
#' @param crs An optional specifier for the coordinate reference system of
#'   the scan lines. This can either be an integer EPSG code, a character
#'   string in proj4 format, or a spatial object from which a map projection
#'   can be taken (e.g. a Raster object or an 'sf' spatial data frame).
#'
#' @return An \code{sf} object (data frame) with columns locationid (integer),
#'   lineid (integer) and geometry (LINESTRING).
#'
#' @examples
#' # Some random point locations in a 20x20km area
#' Npts <- 10
#' pts <- matrix(runif(2*Npts, 0, 20000), ncol = 2)
#'
#' # Add identifying labels for locations
#' pts <- data.frame(id = LETTERS[1:Npts], pts)
#' colnames(pts)[2:3] <- c("x", "y")
#'
#' # Generate 80 lines of uniform length 5km from each point.
#' # Line angles default to regular increments.
#' # Set the coordinate reference system to MGA Zone 55
#' # using its EPSG code 28355
#' #
#' lines <- make_scan_lines(80, pts, lengths = 5000, crs = 28355)
#'
#' # Variable line lengths from an exponential distribution
#' # with mean length 5km (95% bounds: 125m - 18.4km).
#' fun <- function(n) rexp(n, 1 / 5000)
#' lines <- make_scan_lines(80, pts, lengths = fun)
#'
#' @export
#'
make_scan_lines <- function(locations, nlines, lengths, angles = NULL, crs = NULL) {
  
  if (ncol(locations) < 2) stop("Argument locations must have at least 2 columns")
  
  if (is.matrix(locations))
    locations <- as.data.frame(locations)
  
  if (ncol(locations) == 2) {
    # X-Y columns only
    locations <- data.frame(locationid = 1:nrow(locations), locations)
  }
  
  # Ensure identifiers are treated as alpha-numeric
  locations[[1]] <- as.character(locations[[1]])
  
  if (is.vector(lengths)) {
    lengths <- .fixed_length_vector(lengths, nlines)
    length_fun <- function(...) lengths
  }
  else if (is.function(lengths)) {
    .check_nargs(lengths, 1)
    length_fun <- lengths
  } else {
    stop("Argument lengths should be a single value, vector or function")
  }
  
  
  if (is.null(angles)) {
    angles <- seq(0, 2*pi, length.out = nlines+1)[-1]
    angle_fun <- function(...) angles
  }
  else if (is.vector(angles)) {
    angles <- .fixed_length_vector(angles, nlines)
    angle_fun <- function(...) angles
  }
  else if (is.function(angles)) {
    .check_nargs(angles, 1)
    angle_fun <- angles
  } else {
    stop("Argument angles should be a single value, vector or function")
  }
  
  
  lines <- lapply(
    1:nrow(locations),
    
    function(i) {
      x0 <- locations[[i,2]]
      y0 <- locations[[i,3]]
      
      lengths <- length_fun(nlines)
      angles <- angle_fun(nlines)
      
      lapply(1:nlines,
             function(k) {
               x1 <- x0 + lengths[k] * cos(angles[k])
               y1 <- y0 + lengths[k] * sin(angles[k])
               st_linestring(matrix(c(x0, x1, y0, y1), ncol = 2))
             })
    })

  # flatten the list of lists of LINESTRINGS
  lines <- unlist(lines, recursive = FALSE)
  
  dat <- st_sf(locationid = rep(locations[[1]], each = nlines),
               lineid = rep(1:nlines, nrow(locations)),
               geometry = st_sfc(lines),
               stringsAsFactors = FALSE)
  
  if (!is.null(crs)) {
    if (is.integer(crs) | is.character(crs))
      st_crs(dat) <- crs[1]
    else if (inherits(crs, "Raster"))
      st_crs(dat) <- raster::crs(crs, asText = TRUE)
    else if (inherits(crs, "sf"))
      st_crs(dat) <- st_crs(crs)
  }
  
  dat
}




#' Sample one or more raster layers at points along scan lines
#'
#' This function takes one or more raster layers and samples cell values at regularly
#' spaced points along a given set of scan lines. The scan lines are provided
#' as an \code{sf} object as produced by \code{\link{make_scan_lines}}.
#'
#' @note If a raster object has multiple layers (ie. a RasterStack or RasterBrick
#'   object) a warning will be issued and only the first layer will be sampled.
#'
#' @param x Either a single Raster object or a list of Raster objects to sample.
#'   If a named list, the names will be used as the column names for sample values
#'   in the returned data frame.
#'
#' @param lines An \code{sf} object containing scan lines for point locations.
#'
#' @param spacing Spacing between adjacent sample points. If \code{NULL}
#'   (the default), this will be set to half the cell width of the raster layer.
#'
#' @return An \code{sf} object with columns: locationid, lineid (both taken from the input
#'   \code{lines} object), sampleid, geometry (sample point) and a column of
#'   sample values for each of the input rasters.
#'
#' @examples
#' \dontrun{
#'
#' # First example: sample a single raster layer of time since fire values
#' vals <- sample_raster(r.tsf, lines)
#'
#' # Calculate median time since fire value for the central point
#' # of each set of scan lines
#' library(dplyr)
#'
#' vals %>%
#'   as.data.frame() %>%
#'   group_by(locationid) %>%
#'   summarize(tsf = median(value, na.rm = TRUE))
#'
#'
#' # Second example: sample two raster layers for time since fire and
#' # presence of forest cover
#' rr <- list(r.tsf, r.forest)
#' vals <- sample_raster(rr, lines)
#' }
#'
#' @export
#'
sample_raster <- function(x, lines, spacing = NULL) {

  if (inherits(x, "Raster")) {
    x <- list(x)
    names(x) <- deparse(substitute(x))
  }
  else if (is.list(x)) {
    # check all list elements are rasters
    is.r <- sapply(x, function(obj) inherits(obj, "Raster"))
    if (!all(is.r)) stop("All elements in the input list must be Raster objects")

    # check names and set if missing
    nm <- names(x)
    if (is.null(nm)) names(x) <- paste0("layer", 1:length(x))
    else {
      blanks <- nm == ""
      if (any(blanks)) names(x)[blanks] <- paste0("layer", which(blanks))
    }
  }
  else {
    stop("x must be either a Raster object or a list of Raster objects")
  }


  x <- lapply(x, function(r) {
    if (raster::nlayers(r) > 1) {
      warning("Only sampling first layer of multi-layer raster")
      r <- raster::subset(r, 1)
    }
    else r
  })


  cellres <- sapply(x, raster::res)
  if (is.null(spacing)) spacing <- min(cellres)


  # Generate sample points. This will give an sfc object
  # containing MULTIPOINTS, one for each scan line
  mpts <- st_line_sample(lines, density = 1 / spacing)

  # Create a data frame of sample points
  dat <- lapply(1:nrow(lines),
                function(i) {
                  xy <- st_coordinates(mpts[[i]])

                  data.frame(locationid = lines$locationid[i],
                        lineid = lines$lineid[i],
                        sampleid = 1:nrow(xy),
                        x = xy[, 1], 
                        y = xy[, 2])
                })

  dat <- do.call(rbind, dat)


  # Extract values from rasters
  vals <- lapply(x, function(r) {
    raster::extract(r, as.matrix(dat[, c("x", "y")]))
  })
  names(vals) <- names(x)

  # Return result. Note: we don't use cbind here because it doesn't
  # work well with the sf object
  dat <- st_as_sf( data.frame(dat, vals), coords = c("x", "y") )
  st_crs(dat) <- st_crs(lines)
  
  dat
}



#' Calculates line compass bearings
#'
#' Given one or more lines, this function calculates the compass bearing (in
#' degrees) from the start point to the end point of each.
#'
#' @param lines Either a single \code{"LINESTRING"} object or a list of one or more
#'   objects. If the latter, the list can be a geometry column object (class \code{"sfc"})
#'   from an \code{sf} spatial data frame.
#'
#' @return Compass bearings of lines as a numeric vector.
#'
#' @seealso \code{\link[sf]{sfc}}
#'
#' @export
#'
line_bearing <- function(lines) {
  if (inherits(lines, "LINESTRING")) {
    lines <- list(lines)
  }
  else if (inherits(lines, "sfc_LINESTRING")) {
    # no pre-processing required
  }
  else if (is.list(lines)) {
    ok <- all(sapply(X = lines, FUN = inherits, "LINESTRING"))
    if (!ok) stop("All objects in the input list should be class LINESTRING")
  }
  else
    stop("Argument lines should be a LINESTRING (sf geometry) object or a list of LINESTRINGs")


  sapply(lines, function(line) {
    xy <- sf::st_coordinates(line)
    n <- nrow(xy)

    # cartesian angle in degrees
    angle <- 180 * atan2(xy[n, 2] - xy[1, 2], xy[n, 1] - xy[1, 1]) / pi

    # compass bearing
    unname( (450 - angle) %% 360 )
  })
}


#' Tests lines for westerly orientation
#'
#' Given one or more lines, this function checks if the compass bearing from the start
#' to the end point of each lies within the range for westerly orientation.
#'
#' @param lines Either a single \code{"LINESTRING"} object or a list of one or more
#'   objects. If the latter, the list can be a geometry column object (class \code{"sfc"})
#'   from an \code{sf} spatial data frame.
#'
#' @return A logical vector where \code{TRUE} indicates westerly orientation.
#'
#' @export
#'
is_line_west <- function(lines, compass.limits = c(247.5, 292.5)) {
  bearings <- line_bearing(lines)
  bearings >= compass.limits[1] & bearings <= compass.limits[2]
}


#' Generate radial scan lines from one or more locations
#'
#' This function generates a set of \code{n} radial scan lines emanating from
#' each of a given set of locations. Line lengths can be constant; individually
#' specified; or set by a user-provided function (e.g. to draw them from a
#' particular distribution). Line angles can be regularly spaced (the default);
#' individually specified; or set by a user-provided function.
#'
#' @param nlines Number of lines to generate for each point.
#'
#' @param locations One or more locations represented by one of:
#'   \itemize{
#'     \item{A vector of X and Y values for a single location;}
#'     \item{A two-column matrix or data.frame of X and Y values (in that order).}
#'     \item{An sf object with POINT or MULTIPOINT features. In this case the
#'       map projection of the input data (if defined) will be set for the output
#'       scan lines.}
#'   }
#'
#' @param lengths Either a single value for uniform line length; or a
#'   vector of \code{n} line lengths; or a function that returns a vector
#'   of lengths.
#'
#' @param angles Either a vector of line directions (degrees); or a
#'   function that returns a vector of directions; or \code{NULL}
#'   for uniform distribution of line angles.
#'
#' @return An \code{sf} object (data frame) with columns locationid (integer),
#'   lineid (integer) and geometry (LINESTRING).
#'
#' @examples
#' # Some random point locations in a 20x20km area
#' Npts <- 100
#' pts <- matrix(runif(2*Npts, 0, 20000), ncol = 2)
#'
#' # Generate 80 lines of uniform length 5km from each point.
#' # Line angles default to regular increments.
#' lines <- make_scan_lines(80, pts, lengths = 5000)
#'
#' # Variable line lengths from an exponential distribution
#' # with mean length 5km (95% bounds: 125m - 18.4km)
#' fun <- function(n) rexp(n, 1 / 5000)
#' lines <- make_scan_lines(80, pts, lengths = fun)
#'
#' @importFrom sf st_sf st_sfc st_linestring
#'
#' @export
#'
make_scan_lines <- function(nlines, locations, lengths, angles = NULL) {

  if (inherits(locations, "sf")) {
    stop("TODO - add support for sf object")
  }
  else if (is.data.frame(locations))
    centres <- as.matrix(locations[, 1:2])
  else if (is.matrix(locations))
    centres <- locations
  else if (is.vector(locations))
    centres <- matrix(locations[1:2], ncol = 2)

  if (!is.matrix(centres)) {
    stop("Unexpected type of object for argument locations: ", class(locations))
  }


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
    1:nrow(centres),

    function(i) {
      x0 <- centres[i,1]
      y0 <- centres[i,2]
      lengths <- length_fun(nlines)
      angles <- angle_fun(nlines)

      lines <- lapply(1:nlines,
                    function(k) {
                      x1 <- x0 + lengths[k] * cos(angles[k])
                      y1 <- y0 + lengths[k] * sin(angles[k])
                      st_linestring(matrix(c(x0, x1, y0, y1), ncol = 2))
                    })

      st_sf(locationid = i, lineid = 1:nlines, geometry = st_sfc(lines))
    })

  do.call(rbind, lines)
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
#' @importFrom sf st_bind_cols st_cast st_coordinates st_geometry st_sf st_line_sample
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


  # Create an sf object of sample points
  dat <- lapply(1:nrow(lines),
                function(i) {
                  mp <- mpts[[i]]

                  # convert multi-points to a list ('sfc' object) of
                  # simple points
                  p <- st_cast( st_geometry(mp), "POINT" )

                  st_sf(locationid = lines$locationid[i],
                        lineid = lines$lineid[i],
                        sampleid = 1:length(p),
                        geometry = p)
                })

  dat <- do.call(rbind, dat)


  # Extract values from rasters
  mxy <- st_coordinates(dat$geometry)
  vals <- lapply(x, function(r) {
    raster::extract(r, mxy)
  })
  names(vals) <- names(x)

  # Return result. Note: we don't use cbind here because it doesn't
  # work well with the sf object
  st_sf( data.frame(dat, vals) )
}


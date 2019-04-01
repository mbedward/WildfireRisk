
### this function does not need any optimization
###
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

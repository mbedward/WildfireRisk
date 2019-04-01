### investigating to optimize the code
### the fisrt thing is trying to find the bottlenecks of the code
##### for boosting raster lapply we should use snow oackage
#### it seems that library(prioritizr) fast_extract function is usefule for extracting faster 
## fast_extract is not working 
library(here)
library(parallel)
library(doParallel)
library(raster)
library(foreach)
library(microbenchmark)
library(future.apply)
library(prioritizr)

locations <- read.csv( here("par_test/data/ACT_Test_Data/act_cens11_centres.CSV") )
head(locations)
class(locations)
locations[, 1] <- stringr::str_trim( sprintf("%12.0f", locations[, 1]) )



ii <- sample(1:nrow(locations), 2500)
locations <- locations[ii, ]
length_fun <- function(n) 100 + rexp_truncated(n, 1/5000, 20000 - 100)
tsf <- raster( here("par_test/data/act_rasters/actplus_tsf17/hdr.adf") )
dat <- make_scan_line_par(locations, 80, lengths = length_fun, crs = tsf)

## X is the raster objects
## lines is the lines object
sample_raster_endpoints <- function(x, lines) {
  
  if (inherits(x, "Raster")) {
    nm <- deparse(substitute(x))
    x <- list(x)
    names(x) <- nm
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
  
  
  # Coordinates of line distal end-points
  xy <- st_coordinates(lines)[seq(2, 2*nrow(lines), 2), ]

  
  # Extract values from rasters
  vals <- lapply(x, function(r) {
    raster::extract(r, xy[, 1:2])
                    
  })
  str(vals)
  names(vals) <- names(x)
  
  # Return result
  data.frame(locationid = lines$locationid,
             lineid = lines$lineid,
             x = xy[, 1],
             y = xy[, 2],
             vals)
}


sample_raster_endpoints_edited <- function(x, lines) {
  
  if (inherits(x, "Raster")) {
    nm <- deparse(substitute(x))
    x <- list(x)
    names(x) <- nm
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
  
  
  # Coordinates of line distal end-points
  xy <- st_coordinates(lines)[seq(2, 2*nrow(lines), 2), ]
  print(class(xy))
  print(crs(lines,asText = TRUE))
  # Extract values from rasters
  vals <- lapply(x, function(r) {
    prioritizr::fast_extract(r, SpatialPoints(xy[, 1:2]))
  })
  names(vals) <- names(x)
  
  # Return result
  data.frame(locationid = lines$locationid,
             lineid = lines$lineid,
             x = xy[, 1],
             y = xy[, 2],
             vals)
}


profvis::profvis({
  sample_raster_endpoints(tsf,dat)
})

profvis::profvis({
  sample_raster_endpoints_edited(tsf,dat)
})

microbenchmark::microbenchmark(edited = sample_raster_endpoints_edited(tsf,dat),
                               nonedited = sample_raster_endpoints(tsf,dat),
                               times = 10)

sample_raster_endpoints_edited(tsf,dat)


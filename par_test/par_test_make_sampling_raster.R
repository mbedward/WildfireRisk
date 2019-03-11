##This script tries to optimize the code for sample_raster function
## the st_line_string function takes a lot of time to calculate the 
## point from lines. So the source of st_line_sample is obtained from
## github to modify it for parallel processing 
## The problem was calling CLP_gdal_linestring function which is written in c++ is not avaialable to the packages
## finding a way to make it available.
## lapply function which sets coordinate takes less time compared to its parallel version
## although making lapply inside the st_line_sample make the computation far easier
## using data.table::rbindlist seems to be much faster than do.call("rbind",dataframe) so iused rbindlist function
## I should see the result though. Definitley I should use rbindlist function
## running sample_line_point function with parLapplLB is better. 
## the following is the result of the comparison between the edited function and non-edited function
##Unit: seconds
##expr       min         lq       mean     median         uq       max neval
##edited  82.91249   85.13222   90.04062   85.60465   93.83642  102.7173     5
##non_edited 989.49016 1018.26105 1030.91695 1031.38750 1056.47222 1058.9738     5
### the parallel version of the code works better than not parallel but when the number of location is 
## greater than 10. while in this case sample_raster is applied per location in calcualte_risk function and because
## of low number of lines per locationID the parallel processing is not helpful .In fact parallel
## processing just increase the time of calculate_risk function because of overhead. 




library(here)
library(parallel)
library(doParallel)
library(raster)
library(foreach)
library(microbenchmark)
library(data.table)

locations <- read.csv( here("par_test/data/ACT_Test_Data/act_cens11_centres.CSV") )
head(locations)
class(locations)
locations[, 1] <- stringr::str_trim( sprintf("%12.0f", locations[, 1]) )



ii <- sample(1:nrow(locations), 500)
locations <- locations[ii, ]
length_fun <- function(n) 100 + rexp_truncated(n, 1/5000, 20000 - 100)
tsf <- raster( here("par_test/data/act_rasters/actplus_tsf17/hdr.adf") )
dat <- make_scan_line_par(locations, 80, lengths = length_fun, crs = tsf)



### writing the modified version of st_line_sample
## sample_points_line
##using parLapplyLB reduces the computation time alot
## the only problem is that CPL_gdal_linestring_sample is not available
## find a way to make that function available to the package
## this function sampel points from line is an alternative to the 
## st_line_sample function
sample_points_line <- function(x,density){
  
  if (isTRUE(st_is_longlat(x)))
    stop("st_line_sample for longitude/latitude not supported; use st_segmentize?")
  l = st_length(x)
  distList= if(isTRUE(TRUE)){
  n = if (isTRUE(TRUE)) {
    if (!is.na(st_crs(x)) && inherits(density, "units"))
      units(density) = 1/crs_parameters(st_crs(x))$ud_unit # coordinate units
    round(rep(density, length.out = length(l)) * l)
  }
  
  regular = function(n) { (1:n - 0.5)/n }
  c <- parallel::detectCores()
  cl <- parallel::makeCluster(c-2)
  a <- parLapplyLB(cl,seq_along(n), function(i) regular(n[i]) * l[i])
  stopCluster(cl)
  a
  }
x = st_geometry(x)
stopifnot(inherits(x, "sfc_LINESTRING"))
st_sfc(sf:::CPL_gdal_linestring_sample(x, distList), crs = st_crs(x))
  
  
}



sample_raster_edited <- function(x, lines, spacing = NULL) {
  
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
  
  
  cellres <- sapply(x, raster::res)
  if (is.null(spacing)) spacing <- min(cellres)
  
  
  # Generate sample points. This will give an sfc object
  # containing MULTIPOINTS, one for each scan line
  mpts <- sample_points_line(lines, density = 1 / spacing)
  
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

  dat <- data.table::rbindlist(dat,use.names=TRUE, fill=TRUE, idcol=NULL)
  
  
  # Extract values from rasters
  c <- detectCores()
  cl <- makeCluster(c-2)
  vals <- parLapplyLB(cl,x, function(r) {
    raster::extract(r, as.matrix(dat[, c("x", "y")]))
  })
  stopCluster(cl)
  names(vals) <- names(x)
  
  # Return result
  data.frame(dat, vals)
}



sample_raster_edited_notparll <- function(x, lines, spacing = NULL) {
  
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
  
  dat <- data.table::rbindlist(dat,use.names=TRUE, fill=TRUE, idcol=NULL)
  
  
  # Extract values from rasters

  vals <- lapply(x, function(r) {
    raster::extract(r, as.matrix(dat[, c("x", "y")]))
  })
  names(vals) <- names(x)
  
  # Return result
  data.frame(dat, vals)
}




sample_raster <- function(x, lines, spacing = NULL) {
  
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
  
  # Return result
  data.frame(dat, vals)
}







profvis::profvis({
sample_raster_edited(tsf,dat)
})
microbenchmark::microbenchmark(
  edited = sample_raster_edited(tsf,dat),
  non_edited = sample_raster(tsf,dat) ,
  times = 5 )



system.time(sample_raster_edited(tsf,dat))


system.time(sample_raster(tsf,dat))


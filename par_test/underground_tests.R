
library(here)
library(parallel)
library(doParallel)
library(raster)
library(foreach)
library(microbenchmark)

locations <- read.csv( here("par_test/data/ACT_Test_Data/act_cens11_centres.CSV") )
head(locations)
class(locations)
locations[, 1] <- stringr::str_trim( sprintf("%12.0f", locations[, 1]) )



ii <- sample(1:nrow(locations), 500)
locations <- locations[ii, ]
length_fun <- function(n) 100 + rexp_truncated(n, 1/5000, 20000 - 100)
tsf <- raster( here("par_test/data/act_rasters/actplus_tsf17/hdr.adf") )
dat <- make_scan_line_par(locations, 80, lengths = length_fun, crs = tsf)

####sample_raster test code
############
sample_raster_without_docall <- function(x, lines, spacing = NULL) {
  
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
  c <- detectCores()
  cl <- makeCluster(c-2)
  clusterEvalQ(cl,library("sf"))
  dat <- parLapplyLB(cl,1:nrow(lines),
                function(i) {
                  xy <- st_coordinates(mpts[[i]])
                  
                  data.frame(locationid = lines$locationid[i],
                             lineid = lines$lineid[i],
                             sampleid = 1:nrow(xy),
                             x = xy[, 1],
                             y = xy[, 2])
                })
stopCluster(cl)
}





sample_raster_without_docall_lapply <- function(x, lines, spacing = NULL) {
  
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

}




sample_raster_docall_lapply <- function(x, lines, spacing = NULL) {
  
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
  
}



sample_raster_rbindlist_lapply <- function(x, lines, spacing = NULL) {
  
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
  dat <- data.table::rbindlist(dat)
  
}





profvis::profvis({
sample_raster_without_docall(tsf,dat)

})


profvis::profvis({
  sample_raster_without_docall_lapply(tsf,dat)
  
})




profvis::profvis({
  sample_raster_docall_lapply(tsf,dat)
  
})


profvis::profvis({
  sample_raster_rbindlist_lapply(tsf,dat)
  
})










#### TEST for the sample_raster end points

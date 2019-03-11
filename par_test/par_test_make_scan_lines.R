## this file is to investigate the effect of parallel processing in the wildfire risk packagge
## required package is sp, sf, doParallel and doFuture, future, for each, future.apply,
## The only problem that I faced is finding the name of the function that is called inside the length_fun and 
## evaluating that function on the clusters which is the main bottleneck of the code.
## this script contatins the parallelised version of the make_scan_line function. and the result of benchmarking suggest that 
## the function is doing good. The profiling is done by profvis package.
## it seems as the number of observations or number of lines increase parallel version works better while
## for number of lines less than 100 or number of points less than 1000 the non-parallel version is working well
##??????????? I should check that the edited function and the original function have same value

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



ii <- sample(1:nrow(locations), 1500)
locations <- locations[ii, ]
length_fun <- function(n) 100 + rexp_truncated(n, 1/5000, 20000 - 100)
tsf <- raster( here("par_test/data/act_rasters/actplus_tsf17/hdr.adf") )


##### parallel version of make_scan_lines 
## the bottleneck is calling the function inside the length_fun into the clusters

make_scan_line_par <- function(locations, nlines, lengths, angles = NULL, crs = NULL) {
  
  if (ncol(locations) < 2) stop("Argument locations must have at least 2 columns")
  
  if (is.matrix(locations))
    locations <- as.data.frame(locations)
  
  if (ncol(locations) == 2) {
    # X-Y columns only
    locations <- data.frame(locationid = 1:nrow(locations), locations)
  }
  
  # Ensure identifiers are treated as a factor
  locations[[1]] <- factor(locations[[1]])
  
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
  s<- detectCores()
  cl <- makeCluster(s-2)
  registerDoParallel(cl)
  
  lines_1 <- foreach::foreach(x0 = locations[[2]][1:nrow(locations)],y0 = locations[[3]][1:nrow(locations)] ,.packages='sf',.export = "rexp_truncated") %dopar% {
    
    lengths <- length_fun(nlines)
    angles <- angle_fun(nlines) 
    x1<- x0 + lengths*angles
    y1 <- y0 + lengths*angles
    m <- matrix(c(rep(x0,nlines),unlist(x1),rep(y0,nlines),unlist(y1)),ncol = 4)
    res <- lapply(1:nlines, function(i){ st_linestring(matrix(m[i,],ncol = 2))})
  }
  
  
  
  stopCluster(cl)

 # flatten the list of lists of LINESTRINGS
  lines_1 <- unlist(lines_1, recursive = FALSE)
  
  dat_1 <- st_sf(locationid = rep(locations[[1]], each = nlines),
               lineid = rep(1:nlines, nrow(locations)),
               geometry = st_sfc(lines_1),
               stringsAsFactors = FALSE)
  
  if (!is.null(crs)) {
    if (is.integer(crs) | is.character(crs))
      st_crs(dat_1) <- crs[1]
    else if (inherits(crs, "Raster"))
      st_crs(dat_1) <- raster::crs(crs, asText = TRUE)
    else if (inherits(crs, "sf"))
      st_crs(dat_1) <- st_crs(crs)
  }
  
  dat_1
  
}

microbenchmark(
parll = make_scan_line_par(locations, 180, lengths = length_fun, crs = tsf),
notparll = make_scan_lines(locations, 180, lengths = length_fun, crs = tsf),
times = 10
)



##### par Less make scanline
#############
  ##par less make_scan_line


  
 profvis::profvis({
  
  locations
  nlines <- 80
  
  lengths <- length_fun
  
  angles <- NULL 
  
  crs <- tsf
  if (ncol(locations) < 2) stop("Argument locations must have at least 2 columns")
  
  if (is.matrix(locations))
    locations <- as.data.frame(locations)
  
  if (ncol(locations) == 2) {
    # X-Y columns only
    locations <- data.frame(locationid = 1:nrow(locations), locations)
  }
  
  # Ensure identifiers are treated as a factor
  locations[[1]] <- factor(locations[[1]])
  
  if (is.vector(lengths)) {
    lengths <- .fixed_length_vector(lengths, nlines)
    length_fun <- function(...) lengths
  }else if (is.function(lengths)) {
    .check_nargs(lengths, 1)
    length_fun <- lengths
  } else {
    stop("Argument lengths should be a single value, vector or function")
  }
  
  
  if (is.null(angles)) {
    angles <- seq(0, 2*pi, length.out = nlines+1)[-1]
    angle_fun <- function(...) angles
  }else if (is.vector(angles)) {
    angles <- .fixed_length_vector(angles, nlines)
    angle_fun <- function(...) angles
  }else if (is.function(angles)) {
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
})

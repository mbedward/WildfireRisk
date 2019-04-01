### I tyy to avoid writing the calculate risk function using sqllite db to see how things go on 
##
####
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
forest <- raster( here("par_test/data/act_rasters/actfor_comb/hdr.adf") )





### digining body function
lines <- dat
tsf
forest
sample.spacing = 100
dbname = "pointdata.db"
ffdi = 50
kbdi = 100
tsfmax = 50 
  ## foster proportion
  forest <- raster::calc(forest, fun = function(x) !is.na(x) & x != 0)
  
  layers <- list(tsf = tsf, forest = forest)
  
  
  # If the scan lines have a map projection, remove it temporarily
  # so that we don't get pesky errors about incompatible units
  # in the calculations here.
  crs.in <- st_crs(lines)
  st_crs(lines) <- NA_crs_
  
 # dat <-# 
  lines %>%
    group_by(locationid) %>%
    do({
      pdat <- sample_raster_edited_notparll(layers, ., sample.spacing)
      # pdat %>%
      #   group_by(locationid, lineid) %>%
      #   summarize(tsf_mean = mean(tsf, na.rm = TRUE),
      #             forest_p = mean(forest, na.rm = TRUE))
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




system.time(
  risk.lines <- calculate_risk(dat, tsf, forest, sample.spacing = 100)
)

system.time(
  risk.lines_1 <- calculate_risk_edited(dat, tsf, forest, sample.spacing = 100)
  )

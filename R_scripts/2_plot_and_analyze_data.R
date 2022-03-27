#' @author Brian Maitner
#' @description Scripts to wrangle data for the trait shortfalls manuscript

# Load libraries
library(tidyverse)

# Get TRY data


library(sf)
tdwg <- read_sf("manual_downloads/TDWG/level3.shp")
plot(tdwg[1])


#Load in the trait data
  traits <- arrow::open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

#examine the data structure
traits$schema
traits$metadata



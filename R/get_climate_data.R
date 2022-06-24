#climate velocity

library(ClimDatDownloadR)

get_climate_data <- function(out_directory = "data/climate/"){

  # Current climate


  # Future climate for multiple scenarios
    Chelsa.CMIP_5.download(save.location = out_directory,
                           parameter = "bio")


}



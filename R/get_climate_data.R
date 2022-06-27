#climate velocity

library(ClimDatDownloadR)

get_climate_data <- function(out_directory = "data/climate/"){

  # Current climate


  # Future climate for multiple scenarios
    # Chelsa.CMIP_5.download(save.location = out_directory,
    #                        parameter = "bio")

  #This set of scenarios follows Thuiller et al. 2019, Nat. Comm., "Uncertainty in ensembles of global biodiversity scenarios"
  Chelsa.CMIP_5.download(save.location = out_directory,
                         parameter = "bio",
                         model.var = c("CESM1-BGC",
                                       "CMCC-CMS",
                                       "IPSL-CM5A-LR",
                                       "MIROC5",
                                       "MPI-ESM-MR"
                                       ))


}



#next, need a function that iterates through all climate layers/time steps and calcs:
  # country, time, rcp, variable, mean value
  #then we can group and summarize the mean, var, etc.
  # may also need to reproject if these arent equal area (otherwise could bias large counries since high latitude would be disproportionately samples)


MIROC5
ESM-MR

?Chelsa.CMIP_5.download

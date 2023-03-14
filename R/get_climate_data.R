#climate velocity

library(ClimDatDownloadR)
library(arrow)
library(tidyverse)

get_climate_data <- function(out_directory = "data/climate/",
                             out_parquet_directory = "data/climate_averages/",
                             tdwg){


  #ensure output parquet dir exists
    if(!dir.exists(out_parquet_directory)){
      dir.create(out_parquet_directory)
    }

  #skip map downloading if already done


  if(length(list.files(path = out_parquet_directory,pattern = "parquet"))==0){

      crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84"

      #epsgs <- rgdal::make_EPSG()
      crs <- st_crs(6933)



      #Transform tdwg to mollweide
      tdwg <- st_transform(x = tdwg,crs = crs)
      tdwg_sp <- as_Spatial(tdwg)

      # template_rast <-terra::rast(crs = crs$wkt,
      #             resolution =100000,extent = ext(tdwg))


      #plot(tdwg[1])

      #create dir if needed
        if(!dir.exists(out_directory)){
          dir.create(out_directory)
        }


      # Current climate
        Chelsa.Clim.download(parameter = "bio",
                             save.location = out_directory)


      # Future climate for multiple scenarios
        #This set of scenarios follows Thuiller et al. 2019, Nat. Comm., "Uncertainty in ensembles of global biodiversity scenarios"
        Chelsa.CMIP_5.download(save.location = out_directory,
                               parameter = "bio",
                               model.var = c("CESM1-BGC",
                                             "CMCC-CMS",
                                             "IPSL-CM5A-LR",
                                             "MIROC5",
                                             "MPI-ESM-MR"
                                             ))

      #next, need a function that iterates through all climate layers/time steps and calcs:
        # country, time, rcp, variable, mean value
        #then we can group and summarize the mean, var, etc.
        # may also need to reproject if these arent equal area (otherwise could bias large counries since high latitude would be disproportionately samples)




        dirs <- list.dirs(path = file.path(out_directory,"bio"),
                          full.names = FALSE,
                          recursive = FALSE)

        # get relevant layers

          c("bio_V1.2",
            "CESM1-BGC",
            "CMCC-CMS",
            "IPSL-CM5A-LR",
            "MIROC5",
            "MPI-ESM-MR") %>%
            sapply(FUN = function(x){

              dirs[grepl(pattern = x,
                    x = dirs)]

            }) %>%
            unlist() %>%
            as.vector() -> dirs


        # iterate and grab data

          for(i in 1:length(dirs)){

            dir_i <- dirs[i]
            # temp_rasters <- terra::rast(list.files(file.path(out_directory,"bio",dir_i),full.names = TRUE))
            temp_raster_list <- list.files(file.path(out_directory,"bio",dir_i),full.names = TRUE)


          for(j in 1:length(temp_raster_list)){
            message("model ", i)
            message("layer ", j)


            #if output already exists, skip it

            rast_j <- raster::raster(temp_raster_list[j])


            if(file.exists(file.path(out_parquet_directory,paste(names(rast_j),".gz.parquet",sep = "")))){next}



            if(is.na(crs(rast_j))){

              rast_j@crs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

            }


            gc()
            rast_j <- tryCatch(expr = projectRaster(from = rast_j,
                                   res = 10000,
                                   crs = tdwg_sp@proj4string),
                     error = function(e){e}
                     )


            if(inherits(rast_j,"error")){next}
            plot(rast_j)


            vals_j <- extract(x = rast_j,
                              y = tdwg_sp)

            vals_j <- lapply(vals_j, function(x){mean(na.omit(x))})%>%
              unlist()



            out_ij <- data.frame(country = tdwg$LEVEL_3_CO,
                                 variable =names(rast_j),
                                 model = dir_i,
                                 mean_value = vals_j

                                 )

            arrow::write_parquet(x = out_ij,
                                 sink = file.path(out_parquet_directory,paste(names(rast_j),".gz.parquet",sep = "")),
                                 compression = "gzip")
            gc()


          }#j loop



          }# i loop

  }#end stuff that is only done if the parquet files don't exist


    # next load the arrow dataset


    climate_data  <-
      arrow::open_dataset(file.path(out_parquet_directory)) %>%
      collect()%>%
      separate(col = variable,
               sep = "bio",into = c(NA,"bio"),remove = FALSE)%>%
      separate(col = bio,sep = "_",into = c(NA,"bio"),extra = "drop")%>%
      mutate(bio = as.numeric(bio))%>%
        separate(col = model, sep = "_",
                 into = c(NA,"years",NA,NA,NA,NA,"model","rcp"),
                 remove = FALSE)

    climate_data %>%
      mutate(present = case_when(years == "V1.2" ~ "present",
                                 years != "V1.2" ~ "future")) %>%
      filter(present=="present")%>%
      dplyr::select(country,bio,mean_value) %>%
      rename(present_value = mean_value)%>%
      right_join(y = climate_data,by = c("country"="country","bio"="bio"))-> climate_data



    # #summarize % change (or absolute?)

    climate_data %>%
      filter(years!="V1.2")%>%
      mutate(difference = mean_value - present_value,
             pct_difference = (mean_value - present_value)/present_value*100) -> climate_data


    climate_data %>%
      group_by(bio,country,years)%>%
      summarize(mean_difference = mean(difference),
                sd_difference = sd(difference)) -> changes_by_years


    climate_data %>%
      group_by(bio,country)%>%
      summarize(mean_difference = mean(difference),
                sd_difference = sd(difference)) -> changes_overall


    return(list(changes_by_years = changes_by_years,
         changes_overall = changes_overall))




}# end fx






#' @author Brian Maitner
#' @param wcvp wcvp containing appropriate columns
#' @param trait_summary_for_geo_analysis species, traits, # observations
get_georeferenced_trait_coverage <- function(wcvp,tdwg ,trait_summary_for_geo_analysis, temp_file = "temp/temp_georeferenced_trait_coverage.RDS"){

  # Set up output object
    out <- matrix(ncol = 3,nrow = length(unique(wcvp$area_code_l3))*length(unique(trait_summary_for_geo_analysis$TraitName)))
    out <- as.data.frame(out)
    colnames(out) <- c("area","trait", "completeness")


  #Two nested for loops
    counter <- 1
    for(j in 1:length(unique(trait_summary_for_geo_analysis$TraitName))){
    for(i in 1:length(unique(wcvp$area_code_l3))){


          area_i <- unique(wcvp$area_code_l3)[i]

          area_full_i <- tdwg$LEVEL_NAME[which(tdwg$LEVEL_3_CO == area_i)]

          trait_j <- unique(trait_summary_for_geo_analysis$TraitName)[j]

          species_i <- wcvp$taxon_name[which(wcvp$area_code_l3 == area_i)]


          out$area[counter] <- area_i
          out$trait[counter] <- trait_j
          out$completeness[counter] <- length(which(trait_summary_for_geo_analysis$AccSpeciesName %in% species_i &
                                                      trait_summary_for_geo_analysis$TraitName == trait_j &
                                                      trait_summary_for_geo_analysis$LEVEL_NAME == area_full_i))/length(species_i)
          counter <- counter+1
          message(round(counter/nrow(out)*100), " percent done")

          # save the output after each trait is done
            if(i == length(unique(wcvp$area_code_l3))){

              saveRDS(object = out, file = temp_file)

            }


      } # i loop
    } # j loop


  return(out)


}

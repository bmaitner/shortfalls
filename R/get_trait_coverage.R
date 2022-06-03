#' @author Brian Maitner
#' @param wcvp wcvp containing appropriate columns
#' @param trait_summary species, traits, # observations
get_trait_coverage <- function(wcvp, trait_summary, temp_file = "temp/temp_main_trait_coverage.RDS"){

  # Set up output object
    out <- matrix(ncol = 3,nrow = length(unique(wcvp$area_code_l3))*length(unique(trait_summary$TraitName)))
    out <- as.data.frame(out)
    colnames(out) <- c("area", "trait", "completeness")


  #Two nested for loops
    counter <- 1
    for(j in 1:length(unique(trait_summary$TraitName))){
    for(i in 1:length(unique(wcvp$area_code_l3))){


          area_i <- unique(wcvp$area_code_l3)[i]
          trait_j <- unique(trait_summary$TraitName)[j]

          species_i <- wcvp$taxon_name[which(wcvp$area_code_l3 == area_i)]


          out$area[counter] <- area_i
          out$trait[counter] <- trait_j
          out$completeness[counter] <- length(which(trait_summary$AccSpeciesName %in% species_i & trait_summary$TraitName == trait_j))/length(species_i)
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

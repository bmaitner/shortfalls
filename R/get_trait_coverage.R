#' @author Brian Maitner
#' @param wcvp wcvp containing appropriate columns
#' @param trait_summary species, traits, # observations
get_trait_coverage <- function(wcvp, trait_summary){

  # Set up output object
    out <- matrix(ncol = 3,nrow = length(unique(wcvp$area_code_l3))*length(unique(trait_summary$TraitName)))
    out <- as.data.frame(out)
    colnames(out) <- c("area", "trait", "completeness")


  #Two nested for loops
    counter <- 0
    for(i in 1:length(unique(wcvp$area_code_l3))){
      for(j in 1:length(unique(trait_summary$TraitName))){

          area_i <- unique(wcvp$area_code_l3[i])
          trait_j <- unique(trait_summary$TraitName[j])

          species_i <- wcvp$taxon_name[which(wcvp$area_code_l3 == area_i)]


          out$area <- area_i
          out$trait <- trait_j
          out$completeness <- length(which(trait_summary$AccSpeciesName %in% species_i & trait_summary$TraitName == trait_j))/length(species_i)
          counter <- counter+1
          print(counter/nrow(out))


      }
    }


  return(out)


}

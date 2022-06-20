get_family_trait_coverage <- function(wcvp = wcvp,
                                      trait_summary = trait_summary_for_main_analysis,
                                      temp_file = "temp/temp_family_trait_coverage.RDS"){



  # Set up output object
    out <- matrix(ncol = 5,nrow = length(unique(wcvp$family))*length(unique(trait_summary$TraitName)))
    out <- as.data.frame(out)
    colnames(out) <- c("family", "trait", "completeness", "n_species", 'w_data')

  counter <- 1
   for( i in 1:length(unique(wcvp$family))){

     fam_i  <- unique(wcvp$family)[i]
     species_i <- unique(wcvp$taxon_name[which(wcvp$family == fam_i)])

     for(j in 1:length(unique(trait_summary$TraitName))){


        trait_j <- unique(trait_summary$TraitName)[j]

        out$family[counter] <- fam_i
        out$trait[counter] <- trait_j
        out$completeness[counter] <-
          length(unique(trait_summary$AccSpeciesName[which(trait_summary$TraitName == trait_j &
                                                             trait_summary$AccSpeciesName %in% species_i)]))/
          length(species_i)

        out$n_species[counter] <- length(species_i)
        out$w_data[counter] <- length(unique(trait_summary$AccSpeciesName[which(trait_summary$TraitName == trait_j &
                                                                         trait_summary$AccSpeciesName %in% species_i)]))

        counter <- counter+1
        message(round(counter/nrow(out)*100), " percent done")

        # save the output after each trait is done
        if(j == length(unique(trait_summary$TraitName))){

          saveRDS(object = out, file = temp_file)

        }


     }#j
   }#i


  return(out)


}# end fx

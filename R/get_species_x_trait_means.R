

# code to summarize trait data


get_species_x_trait_means <- function(trait_subset){

  # Get metadata

    trait_names <- unique(trait_subset$TraitName)
    species_names <- unique(trait_subset$AccSpeciesName)

  # Iterate over traits, since issues are more likely to be trait-specific

    output <- data.frame(species = species_names)

    traits_to_skip <- c()

    for(i in 1:length(trait_names)){


      trait_i <- trait_names[i]

      data_i <- trait_subset %>% filter(TraitName == trait_i)


    # If data aren't standardized

      if(all(is.na(data_i$StdValue))){

        traits_to_skip <- c(traits_to_skip, trait_i)
        next

      }

    # IF MOST data aren't standardized

      if(length(which(is.na(data_i$StdValue)))/length(data_i$StdValue) > 0.5){

        traits_to_skip <- c(traits_to_skip, trait_i)
        next

      }


    # if data ARE standardized

      if(!all(is.na(data_i$StdValue))){

        # fix multiple units
          if(length(unique(data_i$UnitName)) > 1){

            #toss erroneous data without standardized units
              data_i %>%
                filter(UnitName != "") -> data_i

          }

          if(length(unique(data_i$UnitName)) > 1){

            data_i %>%
              group_by(UnitName) %>%
              count() -> unit_summary_i

            focal_unit_i <- unit_summary_i$UnitName[which.max(unit_summary_i$n)]


            #toss erroneous data without the most common units
            data_i %>%
              filter(UnitName == focal_unit_i) -> data_i

          }



        #check for multiple units
          if(length(unique(data_i$UnitName)) != 1){stop("check units")}

        #get mean values
          data_i %>%
            filter(!is.na(StdValue)) %>%
            group_by(TraitName,AccSpeciesName,UnitName)%>%
            summarise(species_mean = mean(StdValue)) -> data_i


        #re-format
          colnames(data_i)[which(colnames(data_i) == "species_mean")]<-
            paste(unique(data_i$TraitName),"_",unique(data_i$UnitName),sep = "")


        head(data_i)

        data_i %>%
          ungroup()%>%
          select(-TraitName, -UnitName) -> data_i

        merge(x = output,
              y = data_i,
              by.x = "species",
              by.y = "AccSpeciesName",
              all.x = TRUE) -> output


      }







    }# i loop


    return(output)

}




# growth form: select the most common growth forms,do grep, assign the most common growth form
# leaf phenology: select the most common and grep
# dispersal syndrome: huge mess.  probably best ignored
# dispersal unit type: mess
# mycorrhiza type: mess
# ## phyosynthesis pathway: manageable (c3, c4, cam)
# leaf compoundness: fixable, but non-mutually exclusive
# leaf type: potentially fixable, but again, pretty subjective
# nitrogen fixing: potentially fixable
# climate type: fixable
# leaf shape:mess
# leaf margin type: actually pretty good
# leaf distribution: fixable
# Raunkiaer life form: potentially workable
# resprouting: a messs, but could be made to work
# budbank height: fixable
# tolerance to fire: mess
# tolerance to frost: fixable, but may be in days or temperature, or qualitative
# reproductive timing: fixable
# habitat characterization: mess
# ploidy: fixable
# shade tolerance: mess, but could toss the non-standard values if needed
# vegetative regen cap: mess
# root architecture: mess
# plant stem adaptations
# reproduction type: good, but encodes multiple traits
# growth rate: mess

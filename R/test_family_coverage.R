test_family_coverage <- function(family_trait_coverage,trait_list){


  # Set up output object
    out <- matrix(ncol = 6, nrow = length(unique(family_trait_coverage$family))*length(unique(family_trait_coverage$trait)))
    out <- as.data.frame(out)
    colnames(out) <- c("family", "trait","n_spp_in_fam" ,"completeness","exp_completeness","p_value")


  #Two nested for loops
    counter <- 1

    for(i in 1:length(unique(family_trait_coverage$family))){
      for(j in 1:length(unique(family_trait_coverage$trait))){

        trait_j <- unique(family_trait_coverage$trait)[j]
        family_i <- unique(family_trait_coverage$family)[i]

        data_ij <- family_trait_coverage[which(family_trait_coverage$trait == trait_j &
                                      family_trait_coverage$family == family_i),]

        global_ij <- trait_list[which(trait_list$TraitName == trait_j),]


        out$family[counter] <- family_i
        out$trait[counter] <- trait_j
        out$n_spp_in_fam[counter] <- data_ij$n_species
        out$completeness[counter] <- data_ij$completeness
        out$exp_completeness[counter] <- global_ij$pct_coverage_clean/100

        test_i <- binom.test(x = data_ij$w_data,
                   n = data_ij$n_species,
                   p = global_ij$pct_coverage_clean/100,
                   alternative = "two.sided",conf.level = 0.95)

        out$p_value[counter] <- test_i$p.value

        counter <- counter+1
        message(round(counter/nrow(out)*100), " percent done")




      }# j loop
    }# i loop

  return(out)

}

#' @author Brian Maitner

# Load libraries
library(tidyverse)
library(arrow)
source("R/try_to_parquet.R")


##########################################################

# Exclude fern and moss families
library(data.table)
ferns <- fread("manual_downloads/Darwinian_shortfalls/fern_list.txt", quote = "", header = F, col.names = "family") # list of ferns
moss <- fread("manual_downloads/Darwinian_shortfalls/bryophyta_list.txt", quote = "", header = F, col.names = "family") # list of mosses



# Load WCVP
wcvp <- read.table(file = "manual_downloads/WCVP/wcvp_distribution.txt",
                   sep = "|",
                   header = TRUE,
                   quote = "",
                   fill = TRUE,
                   encoding = "UTF-8")

wcvp_names <- read.table(file = "manual_downloads/WCVP/wcvp_names.txt",
                         sep = "|",
                         header = TRUE,
                         quote = "",
                         fill = TRUE,
                         encoding = "UTF-8")


  merge(x = wcvp,
        y = wcvp_names,
        all.x = TRUE) %>%
    filter(taxon_rank == "Species") %>% #include only species
    filter(taxon_status == "Accepted") %>% #only accepted names
    filter(!is.na(accepted_plant_name_id)) %>% #only accepted names
    filter(extinct == 0) %>% #only extant species
    filter(introduced == 0) -> wcvp #exclude introduced


# Remove or correct mistaken Area codes
  shape <- rgdal::readOGR("manual_downloads/Darwinian_shortfalls/level3.shp")
  wcvp$area_code_l3 <- toupper(wcvp$area_code_l3)
  wcvp <- wcvp[wcvp$area_code_l3 %in% shape$LEVEL_3_CO,]


#Toss unneeded columns
wcvp %>%
  select("area_code_l3","introduced","extinct","location_doubtful","taxon_rank","taxon_status","family","genus","species","taxon_name") -> wcvp


rm(wcvp_names)

#Make sure mosses are not present mosses
wcvp <- subset(wcvp, !family %in% moss$family)

#Make a duplicate dataset without ferns (for comparison to the Darwinian paper)
wcvp_no_ferns <- subset(wcvp, !family %in% ferns$family)

saveRDS(object = wcvp, file = "manual_downloads/WCVP/wcvp_cleaned.RDS")
saveRDS(object = wcvp_no_ferns, file = "manual_downloads/WCVP/wcvp_cleaned_no_ferns.RDS",)

n_species <- length(unique(wcvp$taxon_name))


############################################################

# Get TRY data

  # Properly format a list of traits to be pasted into their website
  # This needs to be broken into smaller chunks so that the TRY website can handle it

    # traits <- read.table(file = "manual_downloads/TRY_trait_list_10_03_2022.txt",
    #                      skip = 3,
    #                      sep = "\t",
    #                      stringsAsFactors = FALSE)
    #
    #
    # traits <- traits[,1:5]
    # traits$ObsNum <- as.numeric(traits$ObsNum)
    # colnames(traits) <- traits[1,]
    # traits <- traits[2:nrow(traits),]
    #
    # traits <-  traits[order(traits$ObsNum,decreasing = TRUE),]
    #
    # paste(traits$TraitID[1:10],collapse = ",") # 19811
    # paste(traits$TraitID[11:20],collapse = ",") # 19812
    # paste(traits$TraitID[21:100],collapse = ",") # 19822
    # paste(traits$TraitID[101:1000],collapse = ",") # 19823
    # paste(traits$TraitID[1001:length(traits$TraitID)],collapse = ",") # 19824


  # Convert the dataset to a more useful format

    # 19811
      # try_to_parquet(file = "manual_downloads/TRY/03-11-2022/19811.txt",
      #                output_directory = "manual_downloads/TRY/TRY_parquet/",
      #                batch_size = 70000)

    # 19812
      # try_to_parquet(file = "manual_downloads/TRY/03-11-2022/19812.txt",
      #                output_directory = "manual_downloads/TRY/TRY_parquet/",
      #                batch_size = 80000)



    # 19822
      # try_to_parquet(file = "manual_downloads/TRY/03-11-2022/19822.txt",
      #                output_directory = "manual_downloads/TRY/TRY_parquet/",
      #                batch_size = 80000)

    # 19823
      # try_to_parquet(file = "manual_downloads/TRY/03-11-2022/19823.txt",
      #                output_directory = "manual_downloads/TRY/TRY_parquet/",
      #                batch_size = 80000)
    # 19824
      # try_to_parquet(file = "manual_downloads/TRY/03-11-2022/19824.txt",
      #                output_directory = "manual_downloads/TRY/TRY_parquet/",
      #                batch_size = 80000)

    # "Load" the full dataset (everything in the folder)

      traits <- arrow::open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

    # examine the data structure
      traits$schema
      traits$metadata

      traits %>%
        select(AccSpeciesName) %>%
        collect() %>%
        unique() -> trait_sp_names

    library(TNRS)

    #matching the trait names using WCVP for consistency
    # TNRS(taxonomic_names = trait_sp_names,
    #      sources = "wcvp") -> tnrsed_trait_sp_names
    #saveRDS(object = tnrsed_trait_sp_names,file = "data/tnrsed_names.RDS")
      tnrsed_trait_sp_names <- readRDS("data/tnrsed_names.RDS")

    #generate a list of traits for annotation
      # we want to know the number of species with data

      # traits %>%
      #   filter(!is.na(TraitID)) %>%
      #   group_by(AccSpeciesName, TraitName) %>%
      #   count() %>%
      #   ungroup() %>%
      #   group_by(TraitName) %>%
      #   count() %>%
      #   collect() %>%
      #   mutate(pct_coverage = round(x = (n/n_species) *100,digits = 2)) %>%
      #   write.csv(file = "manual_downloads/trait_annotation/trait_list.csv",
      #             row.names = FALSE)

      trait_list <- read.csv("manual_downloads/trait_annotation/trait_list.csv")
      hist(log10(trait_list$n), breaks = 100, xlab = "log10(n species with data)", main = "")

    # Count the number of observations per species x trait
      traits %>%
        filter(!is.na(TraitID))%>%
        group_by(AccSpeciesName, TraitName) %>%
        count() %>%
        collect() -> trait_summary #1 636 705

      #correct or toss names that cannot be matched

      trait_summary[which(trait_summary$AccSpeciesName %in% wcvp$taxon_name),] -> good_names
      trait_summary[which(!trait_summary$AccSpeciesName %in% wcvp$taxon_name),] -> bad_names

      merge(x = bad_names,
            y = tnrsed_trait_sp_names,
            by.x = "AccSpeciesName",
            by.y = "Name_submitted",
            all.x = TRUE,
            all.y = FALSE) %>%
        filter(Accepted_name_rank == "species") %>% #toss names that couldn't be matched to species
        select(AccSpeciesName, TraitName, n, Accepted_species) %>%
        mutate(AccSpeciesName = Accepted_species) %>%
        select(-Accepted_species) -> bad_names

      #Only keep TNRSed names that match to our list of accepted taxa
      bad_names[which(bad_names$AccSpeciesName %in% wcvp$taxon_name),] -> bad_names

      #Add the totals for the good and bad names together  (needed in case a bad name matched to a good name that was already present)
      rbind(good_names, bad_names) %>%
        group_by(AccSpeciesName, TraitName) %>%
        summarize(n = sum(n)) -> trait_summary # 1 446 171




##########################################################

  # Selecting focal traits

    #update trait list per the new coverage
      trait_summary %>%
        group_by(AccSpeciesName, TraitName) %>%
        count() %>%
        ungroup() %>%
        group_by(TraitName) %>%
        count() %>%
        collect() %>%
        mutate(pct_coverage_clean = (n/n_species)*100) %>%
        rename(n_clean = n)-> trait_list_v2

      merge(x = trait_list, y = trait_list_v2, by = "TraitName", all = TRUE) ->
        trait_list

      rm(trait_list_v2)

      trait_list$n_clean[which(is.na(trait_list$n_clean))] <- 0
      trait_list$pct_coverage_clean[which(is.na(trait_list$pct_coverage_clean))] <- 0


    # How many traits with observation of at least 1%, and are general?
        length(which(trait_list$pct_coverage_clean >= 1 &
                       trait_list$general == 1)) #55

  # Completeness per trait

    #Best coverage
      trait_list[which.max(trait_list$pct_coverage_clean),]

    #Worst coverage
      trait_list[which.min(trait_list$n_clean),]
      length(which(trait_list$n_clean==1))
      length(which(trait_list$n_clean==0))
      length(which(trait_list$n_clean <= 1))

      hist(trait_list$pct_coverage_clean,breaks = 100,main = "Histogram of Trait Coverage",xlab = "Percent Coverage")
      hist(log10(trait_list$pct_coverage_clean),breaks = 100,main = "Histogram of Trait Coverage",xlab = "log(Percent Coverage)")

    #averages

      mean(na.omit(trait_list$pct_coverage_clean))# 0.20 %
      median(na.omit(trait_list$pct_coverage_clean))# 0.0049 %

      Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      }

      Mode(trait_list$pct_coverage_clean)
      Mode(trait_list$n_clean)

    #subset to traits with 1% coverage or more

      trait_list[which(trait_list$pct_coverage_clean >= 1),] %>%
        arrange(pct_coverage_clean)

      trait_list[which(trait_list$pct_coverage_clean >= 1),] %>%
        summarise(mean_coverage = mean(pct_coverage_clean),
                  median_coverage = median(pct_coverage_clean),
                  Mode_coverage = Mode(pct_coverage_clean))

      traits_for_main_analysis <-
        trait_list %>%
        filter(pct_coverage_clean >= 1 & general == 1)

      trait_summary_for_main_analysis <-
        trait_summary %>%
        filter(TraitName %in% traits_for_main_analysis$TraitName)

      saveRDS(object = trait_summary_for_main_analysis,file = "data/trait_summary_for_main_analysis.RDS")
      saveRDS(object = trait_list,file = "data/trait_list_w_coverage.RDS")


###########################################################
  # We need a dataset that lists completeness as a function of country x trait, which we can join to the shapefile for plotting or use in other analyses


  # try: species, trait, number of observations
  # wcvp: species, country


source("R/get_trait_coverage.R")




  trait_coverage <-
    get_trait_coverage(wcvp = wcvp,
                       trait_summary = trait_summary_for_main_analysis)

  saveRDS(object = trait_coverage,
          file = "data/focal_trait_coverage.rds")



  #family_trait_coverage

    # family_trait_coverage <-
    # get_family_trait_coverage(wcvp = wcvp,
    #                           trait_summary = trait_summary_for_main_analysis,
    #                           temp_file = "temp/temp_family_trait_coverage.RDS")
    #
    # saveRDS(object = family_trait_coverage,
    #         file = "data/focal_trait_coverage_family.rds")

  trait_coverage <- read_rds("data/focal_trait_coverage.rds")


  # TDWG polygons from https://github.com/tdwg/wgsrpd/tree/master/level3 on 3/25/2022

  colnames(trait_coverage)

  coverage_wide <-
    trait_coverage %>%
      pivot_wider(id_cols = area,
                  names_from = trait,
                  values_from = completeness)


countries <- sf::read_sf("manual_downloads/TDWG/level3.shp")


plot(countries[1])

countries <-merge(x = countries,
                  y = coverage_wide,
                  by.x = "LEVEL3_COD",
                  by.y = "area")


ggplot(data = countries)+
  geom_sf(aes(fill = 100 * `Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)`))+
  scale_fill_viridis_c(option = "plasma")+
  labs(fill = "%")+
  ggtitle("Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)")


ggplot(data = countries)+
  geom_sf(aes(fill = 100 * `Plant growth form`))+
  scale_fill_viridis_c(option = "plasma")+
  labs(fill = "%")+
  ggtitle("Plant growth form")

ggplot(data = countries)+
  geom_sf(aes(fill = 100 * `Plant height vegetative`))+
  scale_fill_viridis_c(option = "plasma")+
  labs(fill = "%")+
  ggtitle("Plant height vegetative")

ggplot(data = countries)+
  geom_sf(aes(fill = 100 * `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded`))+
  scale_fill_viridis_c(option = "plasma")+
  labs(fill = "%")+
  ggtitle("Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded")

##############################################################################

# Generate a dataset of traits measured within given botanical countries

# First, we'll wrangle the metadata of 3 types: country, state, latitude and longitude


#ObservationID is a primary key
#


traits <- arrow::open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

traits %>%
  select(DataName,TraitID)%>%
  filter(is.na(TraitID)) %>%
  select(DataName)%>%
  collect()%>%
  unique() -> trait_md_options

# Get MD useful for inferring trait location (country or state or lat/long)
  #this will include anything with country, lat/long, or state, EXCEPT where they include qualifiers that cast doubt on the location

  traits %>%
    filter(is.na(TraitID)) %>%
    select(ObservationID, DataName, OrigValueStr)%>%
    filter((grepl(pattern = "country",ignore.case = TRUE,x = DataName)|
              grepl(pattern = "state",ignore.case = TRUE,x = DataName)|
              grepl(pattern = "latitude",ignore.case = TRUE,x = DataName)|
              grepl(pattern = "longitude",ignore.case = TRUE,x = DataName)) &
             !grepl(pattern = "provenance",ignore.case = TRUE,x = DataName) &
             !grepl(pattern = "origin",ignore.case = TRUE,x = DataName) &
             !grepl(pattern = "maximum",ignore.case = TRUE,x = DataName)&
             !grepl(pattern = "minimum",ignore.case = TRUE,x = DataName)) %>%
    collect() %>%
    unique() %>%
    pivot_wider(id_cols = ObservationID,
                names_from = DataName,
                values_from = OrigValueStr) -> useful_md #668k

  source("R/get_countries.R")
  tdwg <- read_sf("manual_downloads/TDWG/old_lv3/level3.shp")

  useful_md <- get_countries(useful_md = useful_md, tdwg = tdwg)

  #How many of the records couldn't be assigned to a bot country?

  stop("code") #og useful m is 668 759 if same number isn't returned, modify to record og number

  useful_md %>%
    select(ObservationID, LEVEL_NAME) %>%
    filter(is.na(LEVEL_NAME)) -> useless_md

  #Toss anything without a LEVEL3 name

    useful_md %>%
      select(ObservationID, LEVEL_NAME) %>%
      filter(!is.na(LEVEL_NAME)) -> useful_md

    message(nrow(useless_md)/(nrow(useless_md)+nrow(useful_md))*100,"% of metadata cannot be used due to errors, etc.")

  #Now, pull only the trait observations with an ObservationID in the useful_md dataset
    #and append country to traits where possible

    traits %>%
      filter(!is.na(TraitID)) %>%
      select(ObservationID, AccSpeciesName, TraitName) %>%
      filter(ObservationID %in% useful_md$ObservationID) %>%
      collect() %>%
        inner_join(y = useful_md,by = "ObservationID")%>%
      group_by(AccSpeciesName, TraitName, LEVEL_NAME) %>%
      count() -> stc #stc = species + trait + country


    #matching the trait names using WCVP for consistency
      TNRS(taxonomic_names = unique(stc$AccSpeciesName),
            sources = "wcvp") -> tnrsed_stc_names

      stc[which(stc$AccSpeciesName %in% wcvp$taxon_name),] -> good_stc_names
      stc[which(!stc$AccSpeciesName %in% wcvp$taxon_name),] -> bad_stc_names

      merge(x = bad_stc_names,
            y = tnrsed_stc_names,
            by.x = "AccSpeciesName",
            by.y = "Name_submitted",
            all.x = TRUE,
            all.y = FALSE) %>%
        filter(Accepted_name_rank == "species") %>% #toss names that couldn't be matched to species
        select(AccSpeciesName, TraitName,LEVEL_NAME, n, Accepted_species) %>%
        mutate(AccSpeciesName = Accepted_species) %>%
        select(-Accepted_species) -> bad_stc_names

      #Toss names that don't match to our list of accepted taxa
        bad_stc_names[which(bad_stc_names$AccSpeciesName %in% wcvp$taxon_name),] -> bad_stc_names

      #Add the totals for the good and bad names together  (needed in case a bad name matched to a good name that was already present)
        rbind(good_stc_names, bad_stc_names) %>%
          group_by(AccSpeciesName, TraitName, LEVEL_NAME) %>%
          summarize(n = sum(n)) -> stc

      #We need to know the number of fraction of trait x species combinations overall
        trait_list_stc <- read.csv("manual_downloads/trait_annotation/trait_list.csv")


        stc %>%
          group_by(AccSpeciesName, TraitName) %>%
          count() %>%
          ungroup() %>%
          group_by(TraitName) %>%
          count() %>%
          collect() %>%
          mutate(pct_coverage_clean = round(x = (n/n_species) *100,digits = 2))%>%
          rename(n_clean = n) -> stc_trait_list_v2

        merge(x = trait_list_stc,
              y = stc_trait_list_v2,
              by = "TraitName", all = TRUE) ->
          stc_trait_list

        rm(stc_trait_list_v2)

        # How many traits with observation of at least 1%, and are general?
        length(which(stc_trait_list$pct_coverage_clean >= 1 &
                       stc_trait_list$general == 1)) #29


        traits_for_geo_analysis <-
          stc_trait_list %>%
          filter(pct_coverage_clean >= 1 & general == 1)

        trait_summary_for_geo_analysis <-
          stc %>%
          filter(TraitName %in% traits_for_geo_analysis$TraitName)

        nrow(trait_summary_for_geo_analysis)


  #run through a country-specific version of get_trait_coverage()


        # georeferenced_trait_coverage <-
        #   get_georeferenced_trait_coverage(wcvp = wcvp,
        #                      trait_summary_for_geo_analysis = trait_summary_for_geo_analysis,
        #                      tdwg = tdwg,
        #                      temp_file = "temp/temp_georeferenced_trait_coverage.RDS")
        #
        # saveRDS(object = georeferenced_trait_coverage,
        #         file = "data/focal_georeferenced_trait_coverage.rds")

        georeferenced_trait_coverage <- read_rds("data/focal_georeferenced_trait_coverage.rds")


        # TDWG polygons from https://github.com/tdwg/wgsrpd/tree/master/level3 on 3/25/2022


        georeferenced_coverage_wide <-
          georeferenced_trait_coverage %>%
          pivot_wider(id_cols = area,
                      names_from = trait,
                      values_from = completeness)


        countries <- sf::read_sf("manual_downloads/TDWG/level3.shp")


        plot(countries[1])

        geo_countries <-merge(x = countries,
                          y = georeferenced_coverage_wide,
                          by.x = "LEVEL3_COD",
                          by.y = "area")


        ggplot(data = geo_countries)+
          geom_sf(aes(fill = 100 * `Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)`))+
          scale_fill_viridis_c(option = "plasma")+
          labs(fill = "%")+
          ggtitle("Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)")


        #get overall trait coverage
          # wcvp x n traits =  expected
          stc %>%
            group_by(TraitName)%>%
            summarise(n_obs = n())%>%
            mutate(frac_cov = n_obs/(nrow(wcvp)) )-> stc_overall_coverage
          mean(stc_overall_coverage$frac_cov)*100
          max(stc_overall_coverage$frac_cov)
          stc_overall_coverage %>% slice_max(order_by = frac_cov)

          #How many trais even had georeferenced data?
            length(unique(stc_overall_coverage$TraitName))

        # Geo focal trait coverage
        georeferenced_trait_coverage %>%
          slice_max(order_by = completeness) #66.6 % coverage (antarctica)

        mean(georeferenced_trait_coverage$completeness)*100

############################################

  #Wood traits: would need to figure out what the woody species are by

    # any species listed as having a trait measurement for wood
    # any species with growth form implying woodiness
    # species noted as woody

        traits <- arrow::open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

        traits %>%
          group_by(TraitName) %>%
          summarize(count = n()) %>%
          collect() -> trait_counts

        trait_counts%>%
          filter(grepl(pattern = "wood", ignore.case = TRUE,x = TraitName)|
                   grepl(pattern = "tree", ignore.case = TRUE,x = TraitName))%>%
          filter(TraitName != "Stem specific density (SSD) or wood density (stem dry mass per stem fresh volume)")%>%
          filter(TraitName != "Plant woodiness") -> wood_traits


        traits %>%
          select(AccSpeciesName,TraitName) %>%
          filter(TraitName %in% wood_traits$TraitName) %>%
          select(AccSpeciesName) %>%
          collect() %>%
          unique() -> trs #these are species whose traits imply woodiness


        traits %>%
          filter(TraitName == "Plant growth form") %>%
          select(AccSpeciesName, TraitName, OrigValueStr) %>%
          collect() %>%
          filter(grepl(pattern = "tree",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "shrub",x = OrigValueStr,ignore.case = TRUE)|
                 grepl(pattern = "liana",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "woody",x = OrigValueStr,ignore.case = TRUE)) %>%
          select(AccSpeciesName) %>%
          unique() -> gfs #species explicitlt stated as having a woody growth form

        traits %>%
          filter(TraitName == "Plant woodiness") %>%
          select(AccSpeciesName, TraitName, OrigValueStr) %>%
          collect()%>%
          filter(grepl(pattern = "woody",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "woody/nonwoody",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "w",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "W",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "wood at base",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "woody at base",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "woody rootstock" ,x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "True" ,x = OrigValueStr,ignore.case = TRUE)) %>%
          select(AccSpeciesName)-> wps #species stated as woody



      putative_wood <- bind_rows(gfs,trs,wps) %>% unique()

        rm(gfs,trs,wps)

      putative_wood <- TNRS::TNRS(taxonomic_names = putative_wood$AccSpeciesName,
                                  sources = "wcvp") %>%
                        filter(Accepted_species != "")


        # Load WCVP
        wcvp <- read.table(file = "manual_downloads/WCVP/wcvp_distribution.txt",
                           sep = "|",
                           header = TRUE,
                           quote = "",
                           fill = TRUE,
                           encoding = "UTF-8")

        wcvp_names <- read.table(file = "manual_downloads/WCVP/wcvp_names.txt",
                                 sep = "|",
                                 header = TRUE,
                                 quote = "",
                                 fill = TRUE,
                                 encoding = "UTF-8")


        merge(x = wcvp,
              y = wcvp_names,
              all.x = TRUE) %>%
          filter(taxon_rank == "Species") %>% #include only species
          filter(taxon_status == "Accepted") %>% #only accepted names
          filter(!is.na(accepted_plant_name_id)) %>% #only accepted names
          filter(extinct == 0) %>% #only extant species
          filter(introduced == 0) -> wcvp #exclude introduced


        # Remove or correct mistaken Area codes
        shape <- rgdal::readOGR("manual_downloads/Darwinian_shortfalls/level3.shp")
        wcvp$area_code_l3 <- toupper(wcvp$area_code_l3)
        wcvp <- wcvp[wcvp$area_code_l3 %in% shape$LEVEL_3_CO,]

          rm(wcvp_names, shape)
          gc()

      # Subset WCVP to woody species

        wcvp_wood <-
        wcvp %>%
          filter(taxon_name %in% putative_wood$Accepted_species)

        # saveRDS(object = wcvp_wood,
        #         file = "manual_downloads/WCVP/wcvp_cleaned_woody.RDS")

        wcvp_wood <- readRDS(file = "manual_downloads/WCVP/wcvp_cleaned_woody.RDS")
        rm(wcvp)

      #need trait summary for wood traits

        traits %>%
          filter(grepl(pattern = "wood", ignore.case = TRUE,x = TraitName)) %>% #filter to only traits pertaining to wood
          filter(TraitName != "Stem specific density (SSD) or wood density (stem dry mass per stem fresh volume)") %>% #toss this trait which could also be applied to non-wood
          filter(TraitName != "Plant woodiness") %>% #toss this trait which could also be applied to non-wood
          select(AccSpeciesName,TraitName) %>%
          group_by(AccSpeciesName, TraitName)%>%
          summarize(n = n())%>%
          collect() -> wood_trait_summary

        #check/fix wood species names
          good_wood_summary <- wood_trait_summary[which(wood_trait_summary$AccSpeciesName %in% wcvp_wood$taxon_name),]
          bad_wood_summary <- wood_trait_summary[which(!wood_trait_summary$AccSpeciesName %in% wcvp_wood$taxon_name),]


          bad_wood_summary%>%
            dplyr::select(AccSpeciesName)%>%
            ungroup()%>%
            unique()%>%
            mutate(ID=row_number())%>%
            select(ID,AccSpeciesName)-> bad_wood_names

          TNRSed_bad_wood_names <- TNRS::TNRS(taxonomic_names = bad_wood_names,
                                              sources = "wcvp")

          bad_wood_names <-
          bad_wood_names %>%
            inner_join(y = TNRSed_bad_wood_names%>%mutate(ID=as.numeric(ID)),by = "ID")

          bad_wood_summary <-
          bad_wood_names %>%
            select(AccSpeciesName, Accepted_species) %>%
            right_join(bad_wood_summary)%>%
            dplyr::select(-AccSpeciesName)%>%
            rename("AccSpeciesName"  = Accepted_species)


          wood_trait_summary <-
          good_wood_summary%>%
            bind_rows(bad_wood_summary)%>%
            group_by(AccSpeciesName,TraitName)%>%
            summarise(n=sum(n))%>%
            filter(AccSpeciesName != "")

          rm(good_wood_summary,bad_wood_names,bad_wood_names)


          source("R/get_trait_coverage.R")


        # wood_trait_coverage <-
        #   get_trait_coverage(wcvp = wcvp_wood,
        #                      trait_summary = wood_trait_summary)
        #
        # saveRDS(object = wood_trait_coverage,
        #         file = "data/wood_trait_coverage.rds")

          wood_trait_coverage <- readRDS(file = "data/wood_trait_coverage.rds")

        # Need to also look how good any of the coverage is for wood traits

          n_woody_species <- length(unique(wcvp_wood$taxon_name))


          wood_trait_summary %>%
            group_by(AccSpeciesName, TraitName) %>%
            count() %>%
            ungroup() %>%
            group_by(TraitName) %>%
            count() %>%
            collect() %>%
            mutate(pct_coverage_clean = (n/n_woody_species)*100) %>%
            rename(n_clean = n) -> woody_trait_list

          n_woody_species

        # focal dataset as per the others

          wood_traits_focal <-
          wood_trait_coverage %>%
            filter(trait %in% woody_trait_list$TraitName[which(woody_trait_list$pct_coverage_clean >= 1)])


          # saveRDS(object = wood_traits_focal,
          #         file = "data/focal_wood_trait_coverage.rds")


###############################################

  #Seed traits and flower traits - by taxonomy

        # Load WCVP
          wcvp <- read.table(file = "manual_downloads/WCVP/wcvp_distribution.txt",
                             sep = "|",
                             header = TRUE,
                             quote = "",
                             fill = TRUE,
                             encoding = "UTF-8")

          wcvp_names <- read.table(file = "manual_downloads/WCVP/wcvp_names.txt",
                                   sep = "|",
                                   header = TRUE,
                                   quote = "",
                                   fill = TRUE,
                                   encoding = "UTF-8")


          merge(x = wcvp,
                y = wcvp_names,
                all.x = TRUE) %>%
            filter(taxon_rank == "Species") %>% #include only species
            filter(taxon_status == "Accepted") %>% #only accepted names
            filter(!is.na(accepted_plant_name_id)) %>% #only accepted names
            filter(extinct == 0) %>% #only extant species
            filter(introduced == 0) -> wcvp #exclude introduced


        # Remove or correct mistaken Area codes
          shape <- rgdal::readOGR("manual_downloads/Darwinian_shortfalls/level3.shp")
          wcvp$area_code_l3 <- toupper(wcvp$area_code_l3)
          wcvp <- wcvp[wcvp$area_code_l3 %in% shape$LEVEL_3_CO,]

          rm(wcvp_names, shape)
          gc()


        # Remove non-seed plants

          fams <- unique(wcvp$family)

          BIEN_fams <- BIEN::BIEN_taxonomy_family(fams)
          BIEN_fams <- BIEN_fams %>%
            select(order,scrubbed_family,higher_plant_group)%>%
            unique()


          # A few manual corrections
              BIEN_fams$higher_plant_group[which(is.na(BIEN_fams$higher_plant_group)&
                                                   BIEN_fams$scrubbed_family == "Ranunculaceae")] <- "flowering plants"

              BIEN_fams$higher_plant_group[which(is.na(BIEN_fams$higher_plant_group)&
                                                   BIEN_fams$scrubbed_family == "Mazaceae")] <- "flowering plants"

              BIEN_fams$higher_plant_group[which(is.na(BIEN_fams$higher_plant_group)&
                                                   BIEN_fams$scrubbed_family == "Lindsaeaceae")] <- "ferns and allies"

              BIEN_fams$higher_plant_group[which(BIEN_fams$order=="Malpighiales"&
                                                   BIEN_fams$scrubbed_family == "Pteridaceae")] <- "ferns and allies"


          seed_fams <-
          BIEN_fams$scrubbed_family[which(BIEN_fams$higher_plant_group %in%
                                            c("flowering plants","gymnosperms (conifers)","gymnosperms (non-conifer)"))]
          seed_fams <- unique(seed_fams)


          wcvp_seed <-
            wcvp %>%
            filter(family %in% seed_fams)

          # saveRDS(object = wcvp_seed,
          #         file = "manual_downloads/WCVP/wcvp_cleaned_seed.RDS")

          wcvp_seed <- readRDS(file = "manual_downloads/WCVP/wcvp_cleaned_seed.RDS")


          #we'll also do flowers while we're at it

          flower_fams <-
            BIEN_fams$scrubbed_family[which(BIEN_fams$higher_plant_group %in%
                                              c("flowering plants"))]
          flower_fams <- unique(flower_fams)

          wcvp_flower <-
            wcvp %>%
            filter(family %in% flower_fams)

          saveRDS(object = wcvp_flower,
                  file = "manual_downloads/WCVP/wcvp_cleaned_flower.RDS")

          wcvp_flower <- readRDS(file = "manual_downloads/WCVP/wcvp_cleaned_flower.RDS")

          rm(wcvp,seed_fams,flower_fams,fams,BIEN_fams)



      # Identify traits

          traits <- arrow::open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

          traits %>%
            group_by(TraitName) %>%
            summarize(count = n()) %>%
            collect() -> trait_counts



          traits %>%
            filter(grepl(pattern = "seed", ignore.case = TRUE,x = TraitName)|
                     grepl(pattern = "flower", ignore.case = TRUE,x = TraitName)|
                     grepl(pattern = "inflorescence", ignore.case = TRUE,x = TraitName)) %>% #filter to only traits pertaining to wood
            select(AccSpeciesName,TraitName) %>%
            group_by(AccSpeciesName, TraitName)%>%
            summarize(n = n())%>%
            collect() -> flower_and_seed_trait_summary

          #check/fix wood species names
          good_flower_and_seed_summary <- flower_and_seed_trait_summary[which(flower_and_seed_trait_summary$AccSpeciesName %in% wcvp_flower$taxon_name|
                                                                                flower_and_seed_trait_summary$AccSpeciesName %in% wcvp_seed$taxon_name  ),]

          bad_flower_and_seed_summary <- flower_and_seed_trait_summary[which(!(flower_and_seed_trait_summary$AccSpeciesName %in% wcvp_flower$taxon_name|
                                                                               flower_and_seed_trait_summary$AccSpeciesName %in% wcvp_seed$taxon_name)),]



          #nrow(good_flower_and_seed_summary)+nrow(bad_flower_and_seed_summary)==nrow(flower_and_seed_trait_summary)


          bad_flower_and_seed_names <-
          bad_flower_and_seed_summary %>%
            dplyr::select(AccSpeciesName) %>%
            ungroup() %>%
            unique() %>%
            mutate(ID=row_number()) %>%
            select(ID,AccSpeciesName)

          TNRSed_bad_flower_and_seed_names <-
            TNRS::TNRS(taxonomic_names = bad_flower_and_seed_names,
                       sources = "wcvp")

          bad_flower_and_seed_names <-
            bad_flower_and_seed_names %>%
            inner_join(y = TNRSed_bad_flower_and_seed_names %>%
                         mutate(ID=as.numeric(ID)), by = "ID")

          bad_flower_and_seed_summary <-
            bad_flower_and_seed_names %>%
            select(AccSpeciesName, Accepted_species) %>%
            right_join(bad_flower_and_seed_summary) %>%
            dplyr::select(-AccSpeciesName) %>%
            rename("AccSpeciesName"  = Accepted_species)

          flower_and_seed_trait_summary <-
            good_flower_and_seed_summary %>%
            bind_rows(bad_flower_and_seed_summary) %>%
            group_by(AccSpeciesName,TraitName) %>%
            summarise(n=sum(n)) %>%
            filter(AccSpeciesName != "")

          rm(good_flower_and_seed_summary,
             bad_flower_and_seed_names,
             bad_flower_and_seed_summary)


        #split into flower summary and seed summary

          flower_and_seed_trait_summary %>%
          filter(grepl(pattern = "seed", ignore.case = TRUE,x = TraitName)) %>% #filter to only traits pertaining to seeds
            select(AccSpeciesName,TraitName) %>%
            group_by(AccSpeciesName, TraitName)%>%
            summarize(n = n())%>%
            collect() -> seed_trait_summary


          seed_traits_to_omit <-
            c("Fruit/seed conspicuous",
              "Fruit/seed color",
              "Fruit/seed abundance",
              "Seedling vigor",
              "Plant morphological adaptations: seed or dispersal unit metamorphoses",
              "Seed number per inflorescence (total, fertile, infertile)",
              "Seed number per flower",
              "Seed number per ramet",
              "Seed or fruit biomass per ground area",
              "Seed mass per inflorescence")


          seed_trait_summary <-
          seed_trait_summary %>%
            filter(!TraitName %in% seed_traits_to_omit)


          n_seed_species <- length(unique(wcvp_seed$taxon_name))


          seed_trait_summary %>%
            group_by(AccSpeciesName, TraitName) %>%
            count() %>%
            ungroup() %>%
            group_by(TraitName) %>%
            count() %>%
            collect() %>%
            mutate(pct_coverage_clean = (n/n_seed_species)*100) %>%
            rename(n_clean = n) -> seed_trait_list

          n_seed_species

          seed_traits_focal <-
            seed_trait_coverage %>%
            filter(trait %in% seed_trait_list$TraitName[which(seed_trait_list$pct_coverage_clean >= 1)])


          # seed_trait_coverage <-
          #   get_trait_coverage(wcvp = wcvp_seed,
          #                      trait_summary = seed_trait_summary)

          # saveRDS(object = seed_trait_coverage,
          #         file = "data/seed_trait_coverage.rds")


          # focal dataset as per the others

          seed_traits_focal <-
            seed_trait_coverage %>%
            filter(trait %in% seed_trait_list$TraitName[which(seed_trait_list$pct_coverage_clean >= 1)])


          # saveRDS(object = seed_traits_focal,
          #         file = "data/focal_seed_trait_coverage.rds")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


            #split into flower summary and seed summary

            flower_and_seed_trait_summary %>%
            filter(grepl(pattern = "flower", ignore.case = TRUE,x = TraitName)|
                     grepl(pattern = "inflorescence", ignore.case = TRUE,x = TraitName) ) %>% #filter to only traits pertaining to seeds
            filter(!grepl(pattern = "litter",ignore.case = TRUE,x = TraitName)) %>% #ignore litter traits that might include flowers
              filter(!grepl(pattern = "Plant dead organ retainment: leaves, branches, flowers",ignore.case = TRUE,x = TraitName)) %>%
            select(AccSpeciesName,TraitName) %>%
            group_by(AccSpeciesName, TraitName)%>%
            summarize(n = n()) -> flower_trait_summary

          # flower_trait_summary%>%
          #   group_by(TraitName)%>%
          #   summarise(n=n())-> flower_trait_counts


          n_flower_species <- length(unique(wcvp_flower$taxon_name))


          flower_trait_summary %>%
            group_by(AccSpeciesName, TraitName) %>%
            count() %>%
            ungroup() %>%
            group_by(TraitName) %>%
            count() %>%
            collect() %>%
            mutate(pct_coverage_clean = (n/n_flower_species)*100) %>%
            rename(n_clean = n) -> flower_trait_list

          n_flower_species


          flower_trait_coverage <-
            get_trait_coverage(wcvp = wcvp_flower,
                               trait_summary = flower_trait_summary)

          saveRDS(object = flower_trait_coverage,
                  file = "data/flower_trait_coverage.rds")

          # focal dataset as per the others

          flower_traits_focal <-
            flower_trait_coverage %>%
            filter(trait %in% flower_trait_list$TraitName[which(flower_trait_list$pct_coverage_clean >= 1)])

          length(unique(flower_traits_focal$trait))


          # saveRDS(object = flower_traits_focal,
          #         file = "data/focal_flower_trait_coverage.rds")
          #















############################################


  #Generating predictor variables

    # Biological variables included
        # species richness,
        # mean species range size, and
        # endemism.

  # Plant species richness was determined as the total number of species for each botanical country.
    # Because species richness was used as an explanatory variable, we did not area-standardize it despite the fact that species richness is clearly affected by the size of the botanical country.

    wcvp %>%
      group_by(area_code_l3) %>%
      summarise(richness = n())


  # Mean species range was calculated as the average of the range size in km2 of the species occurring in each botanical country, where range size was defined as the sum of the area of the countries in which each species occurs.

    socio_vars %>%
      select(LEVEL_3_CO, AREA_SQKM) %>%
      merge(x = wcvp,
            y = .,
            by.x = "area_code_l3",
            by.y = "LEVEL_3_CO",
            all.x = TRUE) %>%
      group_by(taxon_name) %>%
      summarise(range_size = sum(AREA_SQKM)) %>%
      merge(x = wcvp,
            y = .,by = "taxon_name",
            all.x = TRUE) %>%
      group_by(area_code_l3) %>%
      summarise(mean_species_range = mean(na.omit(range_size)))

  # Endemism was calculated as the proportion of species in a botanical country that were not found in any other botanical country.

    wcvp %>%
      group_by(taxon_name) %>%
      summarize(n_countries = n()) %>%
      merge(x = wcvp,
            y = .) %>%
      group_by(area_code_l3) %>%
      filter(n_countries == 1) %>%
    summarise(endemics = n())












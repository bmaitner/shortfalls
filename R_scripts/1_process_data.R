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

      #Toss names that don't match to our list of accepted taxa
      bad_names[which(!bad_names$AccSpeciesName %in% wcvp$taxon_name),] -> bad_names

      #Add the totals for the good and bad names together  (needed in case a bad name matched to a good name that was already present)
      rbind(good_names, bad_names) %>%
        group_by(AccSpeciesName, TraitName) %>%
        summarize(n = sum(n)) -> trait_summary # 1 320 574




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
        mutate(pct_coverage_clean = round(x = (n/n_species) *100,digits = 2))%>%
        rename(n_clean = n)-> trait_list_v2

      merge(x = trait_list,y = trait_list_v2, by = "TraitName", all = TRUE) ->
        trait_list

      rm(trait_list_v2)

    # How many traits with observation of at least 1%, and are general?
        length(which(trait_list$pct_coverage_clean >= 1 &
                       trait_list$general == 1)) #54

  # Completeness per trait

    #Best coverage
      trait_list[which.max(trait_list$pct_coverage_clean),]

    #Worst coverage
      trait_list[which.min(trait_list$n),]
      length(which(trait_list$n==1))

    #averages
      mean(na.omit(trait_list$pct_coverage_clean))# 0.22 %
      median(na.omit(trait_list$pct_coverage_clean))# 0.01 %


    #subset to traits with 1% coverage or more
    trait_list[which(trait_list$pct_coverage_clean >= 1),]

###########################################################
  # We need a dataset that lists completeness as a function of country x trait, which we can join to the shapefile for plotting or use in other analyses


  # try: species, trait, number of observations
  # wcvp: species, country


source("R/get_trait_coverage.R")

    length(which(trait_list$pct_coverage_clean >= 1 &
                   trait_list$general == 1)) #54

    traits_for_main_analysis <-
      trait_list %>%
        filter(pct_coverage_clean >= 1 & general == 1)

    trait_summary_for_main_analysis <-
      trait_summary %>%
        filter(TraitName %in% traits_for_main_analysis$TraitName)


  trait_coverage <-
    get_trait_coverage(wcvp = wcvp,
                       trait_summary = trait_summary_for_main_analysis)

  saveRDS(object = trait_coverage,
          file = "data/focal_trait_coverage.rds")

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

traits %>% head() %>% collect() -> thead


# Get MD useful for inferring trait location (country or state or lat/long)

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
                values_from = OrigValueStr) -> useful_md

  #Next, need to iterate through md and assign a tdwg polygon based on lat/lon if possible, and other metadata if not
    #may need to pull in GNRS here actually

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












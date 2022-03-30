#' @author Brian Maitner

# Load libraries
library(tidyverse)
source("R/try_to_parquet.R")

# Get TRY data

  # Properly format a list of traits to be pasted into their website
  # This needs to be broken into smaller chunks so that the TRY website can handle it

    traits <- read.table(file = "manual_downloads/TRY_trait_list_10_03_2022.txt",
                         skip = 3,
                         sep = "\t",
                         stringsAsFactors = FALSE)


    traits <- traits[,1:5]
    traits$ObsNum <- as.numeric(traits$ObsNum)
    colnames(traits) <- traits[1,]
    traits <- traits[2:nrow(traits),]

    traits <-  traits[order(traits$ObsNum,decreasing = TRUE),]

    paste(traits$TraitID[1:10],collapse = ",") # 19811
    paste(traits$TraitID[11:20],collapse = ",") # 19812
    paste(traits$TraitID[21:100],collapse = ",") # 19822
    paste(traits$TraitID[101:1000],collapse = ",") # 19823
    paste(traits$TraitID[1001:length(traits$TraitID)],collapse = ",") # 19824


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

    #examine the data structure
      traits$schema
      traits$metadata

    # Count the number of observations per species x trait
      traits %>%
        filter(!is.na(TraitID))%>%
        group_by(AccSpeciesName, TraitName) %>%
        count() %>%
        collect() -> trait_summary #1 029 535

    # Take a look at available types of metadata
      traits %>%
        filter(is.na(TraitID)) %>%
        group_by(DataName) %>%
        count() %>%
        collect() -> nontrait_summary

  # Need to subset the main dataset to only observations with country or lat/long data, then use this to fill in country information
    # Once I've got that, I can calculate observations x species x traits
      c("Location Country",
        "Latitude (decimal degrees)",
        "Latitude",
        "Latitude estimated")

##########################################################
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
    select("area_code_l3","introduced","extinct","location_doubtful","taxon_rank","taxon_status","family","genus","species","taxon_name")%>%
    filter(taxon_rank == "Species") -> wcvp

  rm(wcvp_names)



##########################################################

  # We need a dataset that lists completeness as a function of country x trait, which we can join to the shapefile for plotting or use in other analyses


  # try: species, trait, number of observations
  # wcvp: species, country

      # Count the number of observations per species x trait
      traits %>%
        filter(!is.na(TraitID))%>%
        group_by(AccSpeciesName, TraitName) %>%
        count() %>%
        collect() -> trait_summary #1 029 535

source("R/get_trait_coverage.R")

  trait_coverage <-
    get_trait_coverage(wcvp = wcvp,
                       trait_summary = trait_summary)

  saveRDS(object = trait_coverage,
          file = "data/trait_coverage.rds")

  trait_coverage <- read_rds("data/trait_coverage.rds")


  # TDWG polygons from https://github.com/tdwg/wgsrpd/tree/master/level3 on 3/25/2022

  colnames(trait_coverage)

  trait_coverage %>%
    pivot_wider(id_cols = area,
                names_from = trait,
                values_from = completeness)

  trait_coverage %>%
    filter(trait == "Plant growth form")




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
      try_to_parquet(file = "manual_downloads/TRY/03-11-2022/19811.txt",
                     output_directory = "manual_downloads/TRY/TRY_parquet/",
                     batch_size = 80000)

    # 19812
      try_to_parquet(file = "manual_downloads/TRY/03-11-2022/19812.txt",
                     output_directory = "manual_downloads/TRY/TRY_parquet/",
                     batch_size = 80000)

    # 19822
      try_to_parquet(file = "manual_downloads/TRY/03-11-2022/19822.txt",
                     output_directory = "manual_downloads/TRY/TRY_parquet/",
                     batch_size = 80000)

    # 19823
      try_to_parquet(file = "manual_downloads/TRY/03-11-2022/19823.txt",
                     output_directory = "manual_downloads/TRY/TRY_parquet/",
                     batch_size = 80000)

    x <- open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

    x$schema
    x$metadata

    library(tidyverse)

    x %>%
      select(AccSpeciesName) %>%
      collect() -> y






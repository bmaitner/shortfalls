#' @author Brian Maitner

# Load libraries
library(tidyverse)
library(arrow)
source("R/try_to_parquet.R")


##########################################################

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

  # How much of the world is covered by the wcvp?

    # tdwg <- read_sf("manual_downloads/TDWG/old_lv3/level3.shp")
    #
    # tdwg %>%
    #   st_transform(crs = st_crs(6933)) %>%
    #   mutate(area = st_area(geometry)) %>%
    #   mutate(in_wcvp = LEVEL_3_CO %in% wcvp$area_code_l3) -> tdwg
    #
    # tdwg$LEVEL_NAME[which(tdwg$in_wcvp == FALSE)] # Only Bouvet Island doesn't have records.  They also have no vascular plants.
    #
    # tdwg %>%
    #   group_by(in_wcvp) %>%
    #   summarise(total_area = sum(area)) %>%
    #   st_drop_geometry() %>%
    #   ungroup() %>%
    #   mutate(all_area = sum(total_area))%>%
    #   mutate(pct_area = (total_area/all_area) * 100)-> tdwg
    #
    #   tdwg
    #
    #   rm(tdwg)


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
  dplyr::select("area_code_l3","introduced","extinct","location_doubtful","taxon_rank","taxon_status","family","genus","species","taxon_name") -> wcvp


rm(wcvp_names)

#Make sure mosses are not present mosses
wcvp <- subset(wcvp, !family %in% moss$family)

#Make a duplicate dataset without ferns (for comparison to the Darwinian paper)
wcvp_no_ferns <- subset(wcvp, !family %in% ferns$family)

# saveRDS(object = wcvp, file = "manual_downloads/WCVP/wcvp_cleaned.RDS")
# saveRDS(object = wcvp_no_ferns, file = "manual_downloads/WCVP/wcvp_cleaned_no_ferns.RDS",)

n_species <- length(unique(wcvp$taxon_name)) #349 928


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
        dplyr::select(AccSpeciesName,SpeciesName) %>%
        collect() %>%
        unique() -> trait_sp_names

      trait_sp_names %>%
      dplyr::filter(!grepl("unidentified",SpeciesName,ignore.case = TRUE))%>%
        dplyr::filter(!grepl("unidentified",AccSpeciesName,ignore.case = TRUE))-> trait_sp_names

      trait_sp_names %>%
        group_by(SpeciesName) %>%
        mutate(SpeciesName_id = cur_group_id()) %>%
        ungroup()%>%
        group_by(AccSpeciesName)%>%
        mutate(AccSpeciesName_id = cur_group_id()) %>%
        ungroup() -> trait_sp_names


      # traits %>%
      #   select() %>%
      #   collect() %>%
      #   unique() -> trait_sp_names

      traits %>%
        dplyr::select(TraitName) %>%
        collect() %>%
        unique() -> traits_names_og

    library(TNRS)

    #matching the trait names using WCVP for consistency

    # trait_sp_names %>%
    #   select(SpeciesName_id,SpeciesName)%>%
    #   unique() %>%
    # TNRS(sources = "wcvp") -> tnrsed_raw_trait_sp_names

    #saveRDS(object = tnrsed_raw_trait_sp_names,file = "data/tnrsed_raw_names.RDS")

    tnrsed_raw_trait_sp_names <- readRDS("data/tnrsed_raw_names.RDS")

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

      #fix the weird characters that break
        trait_list %>%
          mutate(TraitName = gsub(pattern = "â€¦", replacement = "…",x=TraitName))%>%
          mutate(TraitName = gsub(pattern = "â€\u009d", replacement = "”",x=TraitName))%>%
          mutate(TraitName = gsub(pattern = "â€œ", replacement = "“",x=TraitName))%>%
          mutate(TraitName = gsub(pattern = "â€œ", replacement = "“",x=TraitName))%>%
          mutate(TraitName = gsub(pattern = "Â",
                                  replacement = "",
                                  x=TraitName,fixed = TRUE))-> trait_list


      hist(log10(trait_list$n), breaks = 100, xlab = "log10(n species with data)", main = "")

    # Count the number of observations per species x trait

      #this bit of code makes the first pass at the estimating the number of species per trait

      # this gives us the total number of raw species x trait combinations

        traits %>%
          filter(!is.na(TraitID)) %>%
          group_by(SpeciesName, TraitName, TraitID) %>%
          count() %>%
          collect() %>%
          nrow() # 1 905 810


        traits %>%
          filter(!is.na(TraitID)) %>%
          dplyr::select(SpeciesName) %>%
          collect() %>%
          unique() %>%
          nrow() # 293 045 taxon names (including a lot of junk)

        traits %>%
          filter(!is.na(TraitID)) %>%
          dplyr::select(TraitID) %>%
          collect() %>%
          nrow() # 10 465 028


      traits %>%
        filter(!is.na(TraitID)) %>%
        dplyr::filter(!grepl("unidentified",SpeciesName,ignore.case = TRUE))%>%
        dplyr::filter(!grepl("unidentified",AccSpeciesName,ignore.case = TRUE))%>%
        group_by(SpeciesName, TraitName, TraitID) %>%
        count() %>%
        collect() -> trait_summary # 1 901 820
          #1905810 - 1901820 = 3990 taxon x trait combinations lost here

        #293045 - length(unique(trait_summary$SpeciesName)) #418 names lost by filering "unidentified"


      # this bit of code "corrects" erroneous TraitIDs if multiple exist for the same trait name.  It simply assumes the more common TraitID is correct

      trait_summary %>%
        group_by(TraitName,TraitID) %>%
        summarise(total_n = sum(n)) %>%
        mutate(max_n = max(total_n)) %>%
        filter(total_n == max_n) %>%
        dplyr::select(TraitName,TraitID) %>%
        inner_join(trait_summary %>%
                     rename(raw_TraitID = TraitID)) %>%
        group_by(TraitName,TraitID,SpeciesName)%>%
        summarize(n = sum(n)) -> trait_summary #1 901 815 (lost 5 duplicated taxon x trait combos)

      trait_summary %>%
        dplyr::select(TraitName,TraitID)%>%
        unique() -> tID_lookup

      #make sure names are still all good
        if(any(which(!tID_lookup$TraitName %in% traits_names_og$TraitName))){
         stop("check trait names")
        }

      # connect the TNRSed accepted names to the trait summary

        tnrsed_raw_trait_sp_names %>%
          dplyr::select(Name_submitted,Accepted_species,Name_score,Genus_score,Specific_epithet_score)%>%
          unique()%>%
          right_join(trait_summary,
                     by = c("Name_submitted"="SpeciesName"))-> trait_summary

        #Fix hybrid notation issue
        # wcvp uses "×" for hybrids
        # TNRS uses "x"
        trait_summary$Accepted_species <- gsub(pattern = " x ",
                                               replacement = " × ",
                                               x = trait_summary$Accepted_species)

      #How many raw names match to the WCVP vs TNRSed names?

      trait_summary[which(trait_summary$Name_submitted %in% wcvp$taxon_name),] -> correct_initially #905 120
      trait_summary[which(trait_summary$Accepted_species %in% wcvp$taxon_name),] -> correct_after_TNRS
      correct_after_TNRS %>% filter(Name_score > 0.53) -> correct_after_TNRS #1 822 676

      length(unique(correct_initially$Name_submitted)) #100 154 species initially correct
      length(unique(correct_after_TNRS$Name_submitted)) #271 068 species names correct after TNRS

      length(unique(correct_after_TNRS$Accepted_species)) #131 052 species names match after TNRS
      length(unique(correct_after_TNRS$Accepted_species))/n_species #37.5% of wcvp species


        #for comparison, how many of the "AccSpeciesNames" match directly?
        traits %>%
          filter(!is.na(TraitID)) %>%
          dplyr::select(AccSpeciesName) %>%
          collect() %>%
          unique() -> try_corrected_names

          length(which(try_corrected_names$AccSpeciesName %in% wcvp$taxon_name)) #116 896

          #131052/116896 112%
          #131052-116896 = 14 156 names rescued

          # How many species could not be assigned an accepted species?
              trait_summary %>%
                filter(Accepted_species == "") %>%
                dplyr::select(Name_submitted) %>%
                unique()%>%
                nrow() #28 550 of the raw names could not be matched to a species with TNRS (lots of sp. in here)

          # How many species could not be assigned an accepted species OR matched poorly?
              trait_summary %>%
                filter(Accepted_species == "" | Name_score <= 0.53) %>%
                dplyr::select(Name_submitted) %>%
                unique()%>%
                nrow() #28 749 of the names submitted could not be matched to a species (lots of sp. in here)


      # Filter by name match threshold.  Using .53 based on correspondence with Brad Boyle
        trait_summary %>%
          filter(Accepted_species != "") %>%
          filter(Name_score > 0.53) %>%
          filter(Accepted_species %in% wcvp$taxon_name)-> trait_summary

        #1 752 602

      #28749 + 418 = 29 167 names lost because they couldn't be matched to a species

      # combine the corrected names and toss unneeded columns

        trait_summary %>%
          dplyr::select(Accepted_species, TraitName,n) %>%
          group_by(Accepted_species, TraitName) %>%
          summarise(n = sum(n)) -> trait_summary #1 459 389

      #get counts of # of traits, # of species
        length(unique(trait_summary$Accepted_species)) #131 051
        length(unique(trait_summary$TraitName)) #1916

      trait_summary %>%
        dplyr::select(-n)%>%
        unique()%>%
        nrow() #1 459 389 unique species x trait combinations

      #make sure names are still all good
      if(any(which(!trait_summary$TraitName %in% traits_names_og$TraitName))){
        stop("check trait names")
      }

##########################################################

  # Selecting focal traits

    #update trait list per the new coverage
      trait_summary %>%
        group_by(Accepted_species, TraitName) %>%
        count() %>%
        ungroup() %>%
        group_by(TraitName) %>%
        count() %>%
        collect() %>%
        mutate(pct_coverage_clean = (n/n_species)*100) %>%
        rename(n_clean = n)-> trait_list_v2

      merge(x = trait_list,
            y = trait_list_v2,
            by = "TraitName",
            all = TRUE) ->
        trait_list

      rm(trait_list_v2)

      trait_list$n_clean[which(is.na(trait_list$n_clean))] <- 0
      trait_list$pct_coverage_clean[which(is.na(trait_list$pct_coverage_clean))] <- 0

    # How many traits have at least 1% coverage
      length(which(trait_list$pct_coverage_clean >= 1)) #78
        #2027-82 = 1945 excluded by this

    # How many traits with observation of at least 1%, and are general?
        length(which(trait_list$pct_coverage_clean >= 1 &
                       trait_list$general == 1)) #53

        #78 - 53 = 25 traits excluded by these criteria

  # Completeness per trait

    #Best coverage
      trait_list[which.max(trait_list$pct_coverage_clean),] #32.47

    #Worst coverage
      trait_list[which.min(trait_list$n_clean),] #0
      length(which(trait_list$n_clean==1)) #268
      length(which(trait_list$n_clean==0)) #111
      length(which(trait_list$n_clean <= 1)) #379

      hist(trait_list$pct_coverage_clean,breaks = 100,main = "Histogram of Trait Coverage",xlab = "Percent Coverage")
      hist(log10(trait_list$pct_coverage_clean),breaks = 100,main = "Histogram of Trait Coverage",xlab = "log(Percent Coverage)")

    # Number of SLA vs constituent parts
      trait_summary %>%
        dplyr::filter(TraitName == "Leaf dry mass (single leaf)")%>%
        group_by(Accepted_species)%>%
        count()%>%
        nrow() #6064

      trait_summary %>%
        dplyr::filter(grepl(pattern = "Leaf area per leaf dry mass",
                            x = TraitName))%>%
        group_by(Accepted_species)%>%
        count()%>%
        nrow() #14013

      trait_summary %>%
        dplyr::filter(grepl(pattern = "Leaf area per leaf dry mass",
                            x = TraitName))%>%
        group_by(TraitName)%>%
        summarise(species = length(unique(Accepted_species)))
      #388 to 11256 for each type



      trait_list %>%
        ggplot(mapping = aes(pct_coverage_clean))+
        geom_histogram(fill = "magenta")+
        #xlim(c(0,10000))+
      scale_x_log10(#labels = scales::label_number(accuracy = .0001),
                    labels = c(0,"0.0001",0.001,0.01,.1,1,10,100),
                    breaks=c(0,0.0001,0.001,0.01,.1,1,10,100),
                    limits=c(0.0001,100),
                    expand=expansion(mult = c(0,.01)))+
        #scale_y_continuous(expand = c(0,0))+
        #scale_y_continuous(expand = expansion(add = c(0,5)))+
        scale_y_continuous(expand = expansion(mult = c(0,0.05)))+
        xlab("Percent Completeness (log scale)")+
        ylab("Count")+
        ggtitle("Histogram of Trait Completeness")+
        geom_vline(xintercept = 1,lty=2)+
        theme_bw()-> trait_hist

      ggsave(filename = "plots/trait_histogram.jpg",
             plot = trait_hist,height = 3,width = 5,units = "in")


    #averages

      mean(na.omit(trait_list$pct_coverage_clean))# 0.21 %
      median(na.omit(trait_list$pct_coverage_clean))# 0.0051 %

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

      # mean_coverage median_coverage Mode_coverage
      # 1      3.477565        1.884102      1.347706

      traits_for_main_analysis <-
        trait_list %>%
        filter(pct_coverage_clean >= 1 & general == 1)

      trait_summary_for_main_analysis <-
        trait_summary %>%
        filter(TraitName %in% traits_for_main_analysis$TraitName)


      #How many trait observations?
        sum(trait_summary_for_main_analysis$n) #5 100 952
        length(unique(trait_summary_for_main_analysis$Accepted_species))#122 230
        length(unique(trait_summary_for_main_analysis$Accepted_species))/n_species #34.9%


      # saveRDS(object = trait_summary,file = "data/cleaning_raw_names/trait_summary_overall.RDS")
      # saveRDS(object = trait_summary_for_main_analysis,file = "data/cleaning_raw_names/trait_summary_for_main_analysis.RDS")
      # saveRDS(object = trait_list,file = "data/cleaning_raw_names/trait_list_w_coverage.RDS")
      # saveRDS(object = tID_lookup,file = "data/cleaning_raw_names/tID_lookup.RDS")

      trait_summary_for_main_analysis <-  readRDS("data/cleaning_raw_names/trait_summary_for_main_analysis.RDS")
      trait_list_w_coverage <- readRDS("data/cleaning_raw_names/trait_list_w_coverage.RDS")


      # Per reviewer 2's suggestion, create a table containing traits and coverage for the focal traits

        trait_summary_for_main_analysis %>%
          group_by(TraitName) %>%
          summarise(species_with_data = n()) %>%
          inner_join(y = tID_lookup ) %>%
          arrange(-species_with_data) -> focal_trait_coverage

        # write.csv(x = focal_trait_coverage,
        #           file = "tables/cleaning_raw_names/focal_trait_coverage.csv",row.names = FALSE)

        focal_trait_coverage %>%
          mutate(pct_species_with_data = (species_with_data/n_species)*100) -> focal_trait_coverage

        min(focal_trait_coverage$pct_species_with_data) #1.01 %
        max(focal_trait_coverage$pct_species_with_data) #32.47 %
        mean(focal_trait_coverage$pct_species_with_data) #3.68 %
        median(focal_trait_coverage$pct_species_with_data) #1.89

###########################################################
  # We need a dataset that lists completeness as a function of country x trait, which we can join to the shapefile for plotting or use in other analyses


  # try: species, trait, number of observations
  # wcvp: species, country


# source("R/get_trait_coverage.R")
#
#   trait_coverage <-
#     get_trait_coverage(wcvp = wcvp,
#                        trait_summary = trait_summary_for_main_analysis%>%
#                          rename(AccSpeciesName = Accepted_species))
#
#   saveRDS(object = trait_coverage,
#           file = "data/cleaning_raw_names/focal_trait_coverage.rds")



  # family_trait_coverage

     source("R/get_family_trait_coverage.R")

      # family_trait_coverage <-
      # get_family_trait_coverage(wcvp = wcvp,
      #                           trait_summary = trait_summary_for_main_analysis%>%
      #                             rename(AccSpeciesName = Accepted_species),
      #                           temp_file = "temp/cleaning_raw_names/temp_family_trait_coverage.RDS")
      #
      # saveRDS(object = family_trait_coverage,
      #         file = "data/cleaning_raw_names/focal_trait_coverage_family.rds")

  family_trait_coverage <- readRDS("data/cleaning_raw_names/focal_trait_coverage_family.rds")

  trait_coverage <- read_rds("data/cleaning_raw_names/focal_trait_coverage.rds")


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

# Numbers for trait coverage

  min(trait_coverage$completeness) # 0 pct completeness
  max(trait_coverage$completeness) # 100 pct completeness
  mean(trait_coverage$completeness) #19.4 %
  median(trait_coverage$completeness) # 12.6 %

  mean_completeness_across_country_trait <-
  trait_coverage %>%
    group_by(trait) %>%
    summarise(mean_coverage = mean(completeness))
  min(mean_completeness_across_country_trait$mean_coverage)*100 #6.02
  max(mean_completeness_across_country_trait$mean_coverage)*100 # 71.52
  mean(mean_completeness_across_country_trait$mean_coverage)*100 # 19.36
  median(mean_completeness_across_country_trait$mean_coverage)*100 #15.22

##############################################################################

#Trait coverage across all traits

  source("R/get_trait_coverage.R")
  wcvp <- readRDS(file = "manual_downloads/WCVP/wcvp_cleaned.RDS")
  trait_summary <- readRDS(file = "data/cleaning_raw_names/trait_summary_overall.RDS")

  # total_trait_coverage <-
  # get_trait_coverage(wcvp = wcvp,
  #                    trait_summary = trait_summary %>%
  #                      rename(AccSpeciesName = Accepted_species))
  #
  #
  # saveRDS(object = total_trait_coverage,file = "data/cleaning_raw_names/total_trait_coverage.RDS")
  total_trait_coverage <- readRDS("data/cleaning_raw_names/total_trait_coverage.RDS")



##############################################################################

# Generate a dataset of traits measured within given botanical countries

# First, we'll wrangle the metadata of 3 types: country, state, latitude and longitude


#ObservationID is a primary key
#


traits <- arrow::open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

traits %>%
  dplyr::select(DataName,TraitID)%>%
  filter(is.na(TraitID)) %>%
  dplyr:::select(DataName)%>%
  collect()%>%
  unique() -> trait_md_options

# Get MD useful for inferring trait location (country or state or lat/long)
  #this will include anything with country, lat/long, or state, EXCEPT where they include qualifiers that cast doubt on the location

  traits %>%
    filter(is.na(TraitID)) %>%
    dplyr::select(ObservationID, DataName, OrigValueStr)%>%
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
                values_from = OrigValueStr) -> useful_md #707k

  source("R/get_countries.R")
  tdwg <- read_sf("manual_downloads/TDWG/old_lv3/level3.shp")

  useful_md <- get_countries(useful_md = useful_md, tdwg = tdwg)

  #How many of the records couldn't be assigned to a bot country?

    #og useful md is 707 402 if same number isn't returned, modify to record og number

  useful_md %>%
    dplyr::select(ObservationID, LEVEL_NAME) %>%
    filter(is.na(LEVEL_NAME)) -> useless_md #56 854

  #Toss anything without a LEVEL3 name

    useful_md %>%
      dplyr::select(ObservationID, LEVEL_NAME) %>%
      filter(!is.na(LEVEL_NAME)) -> useful_md #650 548

    message(nrow(useless_md)/(nrow(useless_md)+nrow(useful_md))*100,"% of metadata cannot be used due to errors, etc.")
      #8% of md can't be used

  #Now, pull only the trait observations with an ObservationID in the useful_md dataset
    #and append country to traits where possible

    traits %>%
      filter(!is.na(TraitID)) %>%
      dplyr::select(ObservationID, SpeciesName, TraitName) %>%
      filter(ObservationID %in% useful_md$ObservationID) %>%
      collect() %>%
        inner_join(y = useful_md,by = "ObservationID") -> stc

    #correct names using previous tnrs results

      tnrsed_raw_trait_sp_names <- readRDS("data/tnrsed_raw_names.RDS")

      stc %>%
        inner_join(tnrsed_raw_trait_sp_names %>%
                     dplyr::select(Name_submitted,Accepted_species,Name_score,Genus_score,Specific_epithet_score)%>%
                     unique(),
                   by = c("SpeciesName" = "Name_submitted"))-> stc


      #Fix hybrid notation issue
      # wcvp uses "×" for hybrids
      # TNRS uses "x"
      stc$Accepted_species <- gsub(pattern = " x ",
                                             replacement = " × ",
                                             x = stc$Accepted_species)


      # Filter poor quality names
      # Filter by name match threshold.  Using .53 based on correspondence with Brad Boyle

        stc %>%
          filter(Accepted_species != "") %>%
          filter(Name_score > 0.53) %>%
          filter(Accepted_species %in% wcvp$taxon_name) -> stc

      stc %>%
        dplyr::select(AccSpeciesName = Accepted_species, TraitName, LEVEL_NAME)%>%
      group_by(AccSpeciesName, TraitName, LEVEL_NAME) %>%
      count() -> stc #stc = species + trait + country

      sum(stc$n) #4 191 393 validly georeferenced traits with good names
      nrow(stc) #621 596 unique species x trait x country combinations

      #We need to know the number of fraction of trait x species combinations overall
        trait_list_stc <- read.csv("manual_downloads/trait_annotation/trait_list.csv")

        trait_list_stc %>%
          mutate(TraitName = gsub(pattern = "â€¦", replacement = "…",x=TraitName))%>%
          mutate(TraitName = gsub(pattern = "â€\u009d", replacement = "”",x=TraitName))%>%
          mutate(TraitName = gsub(pattern = "â€œ", replacement = "“",x=TraitName))%>%
          mutate(TraitName = gsub(pattern = "â€œ", replacement = "“",x=TraitName))%>%
          mutate(TraitName = gsub(pattern = "Â",
                                  replacement = "",
                                  x=TraitName,fixed = TRUE))-> trait_list_stc


        #get number of observations retained

          stc %>%
            left_join(trait_list_stc%>%
                        dplyr::select(TraitName,general),by = "TraitName")-> geo_obs_pre_filtering



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

        stc_trait_list$pct_coverage_clean[which(is.na(stc_trait_list$pct_coverage_clean))] <- 0
        stc_trait_list$n_clean[which(is.na(stc_trait_list$n_clean))] <- 0

        # What is the mean coverage across all georeferenced traits?

          mean(stc_trait_list$pct_coverage_clean) #0.067%
          max(stc_trait_list$pct_coverage_clean) #5.57%

        # How many traits have at least one georef value?
          stc_trait_list$TraitName[which(stc_trait_list$n_clean>0)]%>%
            unique()%>%
            length() #1550

          (stc_trait_list$TraitName[which(stc_trait_list$n_clean>0)]%>%
            unique()%>%
            length())/nrow(stc_trait_list) #76.46% of traits have any georeferenced values


        # How many traits with observation of at least 1%, and are general?
        length(which(stc_trait_list$pct_coverage_clean >= 1 &
                       stc_trait_list$general == 1)) #28

        traits_for_geo_analysis <-
          stc_trait_list %>%
          filter(pct_coverage_clean >= 1 & general == 1)

        # coverage stats for the focal geo dataset
        traits_for_geo_analysis[which.max(traits_for_geo_analysis$pct_coverage_clean),]#5.57
        traits_for_geo_analysis[which.min(traits_for_geo_analysis$pct_coverage_clean),]#1.10
        mean(traits_for_geo_analysis$pct_coverage_clean)#1.85
        median(traits_for_geo_analysis$pct_coverage_clean)#1.65


        trait_summary_for_geo_analysis <-
          stc %>%
          filter(TraitName %in% traits_for_geo_analysis$TraitName)

        nrow(trait_summary_for_geo_analysis) #266 583 trait x species x country combos in focal geo ref dataset
        sum(trait_summary_for_geo_analysis$n) #2 017 606

  #run through a country-specific version of get_trait_coverage()
#
#         source("R/get_georeference_trait_coverage.R")
#         georeferenced_trait_coverage <-
#           get_georeferenced_trait_coverage(wcvp = wcvp,
#                              trait_summary_for_geo_analysis = trait_summary_for_geo_analysis,
#                              tdwg = tdwg,
#                              temp_file = "temp/temp_georeferenced_trait_coverage.RDS")
#
        # saveRDS(object = georeferenced_trait_coverage,
        #         file = "data/cleaning_raw_names/focal_georeferenced_trait_coverage.rds")

        georeferenced_trait_coverage <- read_rds("data/cleaning_raw_names/focal_georeferenced_trait_coverage.rds")

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
          select(SpeciesName,TraitName) %>%
          filter(TraitName %in% wood_traits$TraitName) %>%
          select(SpeciesName) %>%
          collect() %>%
          unique() -> trs #these are species whose traits imply woodiness


        traits %>%
          filter(TraitName == "Plant growth form") %>%
          select(SpeciesName, TraitName, OrigValueStr) %>%
          collect() %>%
          filter(grepl(pattern = "tree",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "shrub",x = OrigValueStr,ignore.case = TRUE)|
                 grepl(pattern = "liana",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "woody",x = OrigValueStr,ignore.case = TRUE)) %>%
          select(SpeciesName) %>%
          unique() -> gfs #species explicitlt stated as having a woody growth form

        traits %>%
          filter(TraitName == "Plant woodiness") %>%
          select(SpeciesName, TraitName, OrigValueStr) %>%
          collect()%>%
          filter(grepl(pattern = "woody",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "woody/nonwoody",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "w",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "W",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "wood at base",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "woody at base",x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "woody rootstock" ,x = OrigValueStr,ignore.case = TRUE)|
                   grepl(pattern = "True" ,x = OrigValueStr,ignore.case = TRUE)) %>%
          select(SpeciesName)-> wps #species stated as woody



      putative_wood <- bind_rows(gfs,trs,wps) %>% unique()

        rm(gfs,trs,wps)

      # TNRS using the earlier data

        tnrsed_raw_trait_sp_names <- readRDS("data/tnrsed_raw_names.RDS")

        tnrsed_raw_trait_sp_names %>%
          filter(Name_submitted %in% putative_wood$SpeciesName) -> putative_wood


        #Fix hybrid notation issue
        # wcvp uses "×" for hybrids
        # TNRS uses "x"
        putative_wood$Accepted_species <- gsub(pattern = " x ",
                                     replacement = " × ",
                                     x = putative_wood$Accepted_species)


        # Filter poor quality names
        # Filter by name match threshold.  Using .53 based on correspondence with Brad Boyle

        putative_wood %>%
          filter(Accepted_species != "") %>%
          filter(Name_score > 0.53) %>%
          filter(Accepted_species %in% wcvp$taxon_name) -> putative_wood

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
          filter(SpeciesName %in% putative_wood$Name_submitted)%>%
          filter(grepl(pattern = "wood", ignore.case = TRUE,x = TraitName)) %>% #filter to only traits pertaining to wood
          filter(TraitName != "Stem specific density (SSD) or wood density (stem dry mass per stem fresh volume)") %>% #toss this trait which could also be applied to non-wood
          filter(TraitName != "Plant woodiness") %>% #toss this trait which could also be applied to non-wood
          select(SpeciesName,TraitName) %>%
          group_by(SpeciesName, TraitName)%>%
          summarize(n = n())%>%
          collect() -> wood_trait_summary

        #check/fix wood species names

          wood_trait_summary %>%
            left_join(putative_wood %>%
                        select(Name_submitted, Accepted_species)%>%
                        unique(),
                      by = c("SpeciesName" = "Name_submitted")) -> wood_trait_summary

          wood_trait_summary %>%
            ungroup()%>%
            rename(AccSpeciesName = Accepted_species)%>%
            select(-SpeciesName)%>%
            group_by(TraitName,AccSpeciesName)%>%
            summarize(n = sum(n)) -> wood_trait_summary


        #   source("R/get_trait_coverage.R")
        #
        # wood_trait_coverage <-
        #   get_trait_coverage(wcvp = wcvp_wood,
        #                      trait_summary = wood_trait_summary)
        #
        # saveRDS(object = wood_trait_coverage,
        #         file = "data/cleaning_raw_names/wood_trait_coverage.rds")

          wood_trait_coverage <- readRDS(file = "data/cleaning_raw_names/wood_trait_coverage.rds")

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

          n_woody_species #91552

        # focal dataset as per the others

          wood_traits_focal <-
          wood_trait_coverage %>%
            filter(trait %in% woody_trait_list$TraitName[which(woody_trait_list$pct_coverage_clean >= 1)])


          # saveRDS(object = wood_traits_focal,
          #         file = "data/cleaning_raw_names/focal_wood_trait_coverage.rds")


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

          # saveRDS(object = wcvp_flower,
          #         file = "manual_downloads/WCVP/wcvp_cleaned_flower.RDS")

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
                     grepl(pattern = "inflorescence", ignore.case = TRUE,x = TraitName)) %>%
            select(SpeciesName,TraitName) %>%
            group_by(SpeciesName, TraitName)%>%
            summarize(n = n())%>%
            collect() -> flower_and_seed_trait_summary

        # fix names using tnrs

          tnrsed_raw_trait_sp_names <- readRDS("data/tnrsed_raw_names.RDS")

          flower_and_seed_trait_summary %>%
            left_join(tnrsed_raw_trait_sp_names %>%
                        select(Name_submitted,Accepted_species,Name_score,
                               Genus_score,Specific_epithet_score) %>%
                         unique(),
                       by = c("SpeciesName" = "Name_submitted")) -> flower_and_seed_trait_summary

          flower_and_seed_species <- unique(c(wcvp_flower$taxon_name,wcvp_seed$taxon_name))

          flower_and_seed_trait_summary %>%
            filter(Accepted_species != "") %>%
            filter(Name_score > 0.53) %>%
            filter(Accepted_species %in% flower_and_seed_species) -> flower_and_seed_trait_summary

            flower_and_seed_trait_summary %>%
              rename(AccSpeciesName = Accepted_species) -> flower_and_seed_trait_summary


        #split into flower summary and seed summary

          flower_and_seed_trait_summary %>%
            ungroup()%>%
          filter(grepl(pattern = "seed", ignore.case = TRUE,x = TraitName)) %>% #filter to only traits pertaining to seeds
            select(AccSpeciesName,TraitName,n) %>%
            group_by(AccSpeciesName, TraitName)%>%
            summarize(n = sum(n)) -> seed_trait_summary


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


          n_seed_species <- length(unique(wcvp_seed$taxon_name)) #336 712


          seed_trait_summary %>%
            group_by(AccSpeciesName, TraitName) %>%
            count() %>%
            ungroup() %>%
            group_by(TraitName) %>%
            count() %>%
            collect() %>%
            mutate(pct_coverage_clean = (n/n_seed_species)*100) %>%
            rename(n_clean = n) -> seed_trait_list


        # seed_trait_coverage <-
        #   get_trait_coverage(wcvp = wcvp_seed,
        #                      trait_summary = seed_trait_summary)
        #
        #   saveRDS(object = seed_trait_coverage,
        #           file = "data/cleaning_raw_names/seed_trait_coverage.rds")
        #
        #
        #   seed_traits_focal <-
        #     seed_trait_coverage %>%
        #     filter(trait %in% seed_trait_list$TraitName[which(seed_trait_list$pct_coverage_clean >= 1)])
        #
        #
        #
        #   # focal dataset as per the others
        #
        #   seed_traits_focal <-
        #     seed_trait_coverage %>%
        #     filter(trait %in% seed_trait_list$TraitName[which(seed_trait_list$pct_coverage_clean >= 1)])
        #
        #
        #   saveRDS(object = seed_traits_focal,
        #           file = "data/cleaning_raw_names/focal_seed_trait_coverage.rds")
        #
        #


######################################################################


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


          n_flower_species <- length(unique(wcvp_flower$taxon_name)) #335 537


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


          # flower_trait_coverage <-
          #   get_trait_coverage(wcvp = wcvp_flower,
          #                      trait_summary = flower_trait_summary)

          # saveRDS(object = flower_trait_coverage,
          #         file = "data/cleaning_raw_names/flower_trait_coverage.rds")

          # focal dataset as per the others

          flower_traits_focal <-
            flower_trait_coverage %>%
            filter(trait %in% flower_trait_list$TraitName[which(flower_trait_list$pct_coverage_clean >= 1)])

          length(unique(flower_traits_focal$trait))#8


          # saveRDS(object = flower_traits_focal,
          #         file = "data/cleaning_raw_names/focal_flower_trait_coverage.rds")


##########################################

  #Trait coverage without ferns (for comparison with Rudbeck data)


        source("R/get_trait_coverage.R")
        wcvp_no_ferns <- readRDS("manual_downloads/WCVP/wcvp_cleaned_no_ferns.RDS")
        trait_summary_for_main_analysis<-  readRDS("data/cleaning_raw_names/trait_summary_for_main_analysis.RDS")


        #     trait_coverage_no_ferns <-
        #       get_trait_coverage(wcvp = wcvp_no_ferns,
        #                          trait_summary = trait_summary_for_main_analysis %>%
        #                            rename(AccSpeciesName = Accepted_species))
        #
        # saveRDS(object = trait_coverage_no_ferns,
        #         file = "data/cleaning_raw_names/focal_trait_coverage_no_ferns.rds")





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


############################################

# australia completeness

    #remotes::install_github("traitecoevo/austraits", dependencies = TRUE, upgrade = "ask",build_vignettes = TRUE)
    library(austraits)
    vignette("austraits")

    #version 3.0.2
    aus <- austraits::load_austraits(version = austraits::get_version_latest(),
                                     path = "data/austraits/")


    oz_traits <- aus$traits
    oz_trait_defs <- aus$definitions$traits

###############################




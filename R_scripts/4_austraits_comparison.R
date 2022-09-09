# Comparing Australia w and w/o Austraits


# Load packages

  library(sf)
  library(tidyverse)
  library(austraits)
  library(TNRS)
  library(ggpubr)
  source("R/get_trait_coverage.R")


# Load in data from 1_...

  wcvp <- readRDS("manual_downloads/WCVP/wcvp_cleaned.RDS")
  trait_summary_for_main_analysis <- readRDS("data/trait_summary_for_main_analysis.RDS")

# Load TDWG spatial data
  tdwg <- read_sf("manual_downloads/TDWG/old_lv3/level3.shp")

# Filter tdwg to only Australia

  tdwg %>%
    dplyr::filter(REGION_NAM == "Australia") -> tdwg

# Filter WCVP to Australia

  wcvp %>%
    dplyr::filter(area_code_l3 %in% tdwg$LEVEL_3_CO) -> wcvp

# Filter trait summary for Australian traits
  trait_summary_for_main_analysis %>%
    dplyr::filter(AccSpeciesName %in% wcvp$taxon_name)

# Get AusTraits

  austraits <- austraits::load_austraits(path = "data/austraits/",
                                         version = get_version_latest())


  austraits <- austraits$traits

# Load translation table (thanks Matthias)

  austry <- xlsx::read.xlsx("manual_downloads/austraits/AusTraits-TRY matches.xlsx",1)

# Convert Austraits names to TRY

  austry %>%
    dplyr::select(trait_name,TRY.name)%>%
    right_join(y = austraits,
               by = "trait_name") %>%
    rename(AccSpeciesName = taxon_name,
           TraitName = TRY.name) %>%
    select(AccSpeciesName, TraitName) %>%
    group_by(AccSpeciesName,TraitName)%>%
    summarise(n = n()) -> austraits

  austraits %>%
    dplyr::filter(TraitName %in% unique(trait_summary_for_main_analysis$TraitName)) -> austraits

# check names (note that "good" and "bad" are used for convenience, not a judgement of actual quality)

  austraits %>%
    dplyr::filter(AccSpeciesName %in% wcvp$taxon_name) -> good_names

  austraits %>%
    dplyr::filter(!AccSpeciesName %in% wcvp$taxon_name) -> bad_names

# attempt to resolve bad (== non-matching) names using TNRS

  tnrsed_bad <- TNRS(taxonomic_names = unique(bad_names$AccSpeciesName),
                       sources = "wcvp")

  bad_names %>%
    left_join(y = tnrsed_bad,
              by = c("AccSpeciesName" = "Name_submitted")) -> bad_names

  bad_names %>%
    ungroup %>%
    select(Accepted_species, TraitName, n) %>%
    rename(AccSpeciesName = Accepted_species) -> bad_names

  bad_names %>%
    dplyr::filter(!AccSpeciesName %in% wcvp$taxon_name) -> bad_names

  austraits <- bind_rows(good_names,
                         bad_names)

  rm(austry,bad_names,good_names,tnrsed_bad)

# Get coverage with AusTraits

  all_data <- bind_rows(austraits, trait_summary_for_main_analysis)

  all_data %>%
    group_by(AccSpeciesName,TraitName)%>%
    summarise(n = sum(n)) -> all_data

    trait_coverage_w_austraits <-
      get_trait_coverage(wcvp = wcvp,
                         trait_summary = all_data)

    saveRDS(object = trait_coverage_w_austraits,
            file = "data/focal_trait_coverage_australia_w_austraits.rds")



# Get coverage without AusTraits from the main analysis

  # Load overall trait completeness
  general_traits_one_percent_threshold <- read_rds("data/focal_trait_coverage.rds")

  ###########################################################
  # Mean coverage map

  tdwg_try <-
    general_traits_one_percent_threshold %>%
    group_by(area) %>%
    summarise(mean_coverage_focal = mean (completeness)) %>%
    inner_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))


  tdwg_austry <-
    trait_coverage_w_austraits %>%
    group_by(area) %>%
    summarise(mean_coverage_focal = mean (completeness)) %>%
    inner_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))



plot_try <-
  tdwg_general %>%
    st_transform(crs = st_crs(6933))%>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \nCompleteness\n(%)",
                        limits=c(0,30))+
    theme_minimal()+
  xlim(11000000,15000000)+
  geom_sf_text(mapping = aes(label = round(mean_coverage_focal*100,
                                           digits = 1)))+
  xlab(NULL)+ylab(NULL)+
  ggtitle("TRY")


plot_austry <-
  tdwg_austry %>%
    st_transform(crs = st_crs(6933))%>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \nCompleteness\n(%)",
                        limits=c(0,30))+
    xlim(11000000,15000000)+
    theme_minimal()+
    geom_sf_text(mapping = aes(label = round(mean_coverage_focal*100,
                                             digits = 1)),inherit.aes = TRUE)+
  xlab(NULL)+ylab(NULL)+
  ggtitle("TRY + AusTraits")






ggarrange(plot_try,
          plot_austry,
          ncol = 2,
          common.legend = TRUE)







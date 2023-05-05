# Comparing Australia w and w/o Austraits


# Load packages

  remotes::install_github("traitecoevo/austraits",
                          dependencies = TRUE,
                          upgrade = "ask",
                          build_vignettes = TRUE)
  library(sf)
  library(tidyverse)
  library(austraits)
  library(TNRS)
  library(ggpubr)
  source("R/get_trait_coverage.R")


# Load in data from 1_...

  wcvp <- readRDS("manual_downloads/WCVP/wcvp_cleaned.RDS")
  trait_summary_for_main_analysis <- readRDS("data/cleaning_raw_names/trait_summary_for_main_analysis.RDS")

# Load TDWG spatial data
  tdwg <- read_sf("manual_downloads/TDWG/old_lv3/level3.shp")

# Filter tdwg to only Australia

  tdwg %>%
    dplyr::filter(REGION_NAM == "Australia")%>%
    dplyr::filter(!grepl(pattern = "Is.",x = LEVEL_NAME))-> tdwg

# Filter WCVP to Australia

  wcvp %>%
    dplyr::filter(area_code_l3 %in% tdwg$LEVEL_3_CO) -> wcvp

# Filter trait summary for Australian traits
  trait_summary_for_main_analysis %>%
    rename(AccSpeciesName = Accepted_species)%>%
    dplyr::filter(AccSpeciesName %in% wcvp$taxon_name) -> trait_summary_for_main_analysis

# Get AusTraits

  austraits <- austraits::load_austraits(path = "data/austraits/",
                                         version = "3.0.2")

  austraits <- austraits$traits
  nrow(austraits) #997 808

# TNRS raw names

  library(TNRS)

  austraits %>%
    dplyr::select(original_name)%>%
    unique()%>%
    mutate(ID = 1:n())-> ausnames

    # TNRS::TNRS(taxonomic_names = ausnames[c("ID","original_name")],
    #            sources = "wcvp") -> tnrsed_raw_ausnames
    #
    # saveRDS(object = tnrsed_raw_ausnames,file = "data/cleaning_raw_names/tnrsed_raw_ausnames.rds")
    tnrsed_raw_ausnames <- readRDS("data/cleaning_raw_names/tnrsed_raw_ausnames.rds")


  austraits %>%
    left_join(tnrsed_raw_ausnames %>%
                dplyr::select(Name_submitted,Accepted_species,Name_score,Genus_score,Specific_epithet_score)%>%
                unique(),
              by = c("original_name" = "Name_submitted")) -> austraits

#toss names that don't match the wcvp

  austraits %>%
    dplyr::filter(Name_score > 0.53) %>%
    dplyr::filter(Accepted_species %in% wcvp$taxon_name) -> austraits #883 215 records



#remove the old names so there are no mistakes

  austraits %>%
    dplyr::select(-taxon_name,-original_name)%>%
    rename(taxon_name = Accepted_species) %>%
    dplyr::select(-Name_score, -Genus_score, -Specific_epithet_score) -> austraits


# Load translation table (thanks Matthias)

  austry <- xlsx::read.xlsx("manual_downloads/austraits/AusTraits-TRY matches.xlsx",1)

# Convert Austraits names to TRY

  austry %>%
    dplyr::select(trait_name,TRY.name)%>%
    right_join(y = austraits,
               by = "trait_name")%>%
    rename(AccSpeciesName = taxon_name,
           TraitName = TRY.name)%>%
    select(AccSpeciesName, TraitName) %>%
    group_by(AccSpeciesName,TraitName)%>%
    summarise(n = n()) -> austraits

  austraits %>%
    dplyr::filter(TraitName %in% unique(trait_summary_for_main_analysis$TraitName)) -> austraits

  rm(austry)


#coverage improvement with AusTraits

  nrow(trait_summary_for_main_analysis) # 70,103 species x trait combinations in TRY
  length(unique(trait_summary_for_main_analysis$AccSpeciesName))#9,809 species in TRY
  length(unique(trait_summary_for_main_analysis$TraitName)) #53

  nrow(austraits) # 147,228
  length(unique(austraits$AccSpeciesName)) #18,619
  length(unique(austraits$TraitName)) #34


# Get coverage with AusTraits

  all_data <- bind_rows(austraits, trait_summary_for_main_analysis)

  all_data %>%
    group_by(AccSpeciesName,TraitName) %>%
    summarise(n = sum(n)) -> all_data

    nrow(all_data) #187312
    length(unique(all_data$AccSpeciesName)) #18 866
    length(unique(all_data$TraitName)) #53

    length(unique(all_data$AccSpeciesName))/length(unique(wcvp$taxon_name))#89%

    # trait_coverage_w_austraits <-
    #   get_trait_coverage(wcvp = wcvp,
    #                      trait_summary = all_data)

    # saveRDS(object = trait_coverage_w_austraits,
    #         file = "data/cleaning_raw_names/focal_trait_coverage_australia_w_austraits.rds")

    trait_coverage_w_austraits <- readRDS("data/cleaning_raw_names/focal_trait_coverage_australia_w_austraits.rds")


    length(unique(all_data$AccSpeciesName))/length(unique(wcvp$taxon_name)) #89.0% of species are represented




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

  #crs_wintri <- "+proj=wintri +datum=WGS84 +no_defs +over"
  crs_wintri <- "+proj=wintri +datum=WGS84 +no_defs +over +proj=laea +lon_0=133.59375 +lat_0=-26.0513517 +datum=WGS84 +units=m +no_defs"


  grat_wintri <-
    st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
    st_transform_proj(crs = crs_wintri)



plot_try <-
  tdwg_try %>%
    #st_transform(crs = st_crs(6933))%>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \nCompleteness\n(%)",
                        limits=c(0,30))+
    #theme_minimal()+
  coord_sf(datum = NULL) +
  theme_map()+
  xlim(-2060834,1980902)+
    ylim(-2007507,1739014)+
  geom_sf_text(mapping = aes(label = round(mean_coverage_focal*100,
                                           digits = 1)))+
  xlab(NULL)+
    ylab(NULL)+
  ggtitle("TRY")+
    coord_sf(datum = NULL) +
    theme_map()



plot_austry <-
  tdwg_austry %>%
    #st_transform(crs = st_crs(6933))%>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \nCompleteness\n(%)",
                        limits=c(0,31))+
    #xlim(11000000,15000000)+
    #theme_minimal()+
    coord_sf(datum = NULL) +
    theme_map()+
    xlim(-2060834,1980902)+
    ylim(-2007507,1739014)+
    geom_sf_text(mapping = aes(label = round(mean_coverage_focal*100,
                                               digits = 1)),inherit.aes = TRUE)+
    xlab(NULL)+ylab(NULL)+
  ggtitle("TRY + AusTraits")+
  coord_sf(datum = NULL) +
  theme_map()






ggarrange(plot_try,
          plot_austry,
          ncol = 2,
          common.legend = TRUE)-> austry_plot

ggsave(plot = austry_plot,
       filename = "plots/austraits_plus_try.jpg",width = 10,height = 4.5, bg = "white")


ggsave(plot = austry_plot,
       filename = "plots/fig3_high_res.pdf",width = 10,height = 4.5, bg = "white",dpi = 600)





#' @author Brian Maitner
#' @description Scripts to make maps for the shortfalls paper

# Load libraries
library(tidyverse)
library(lme4)
library(sf)
library(tricolore)


###########################################################

# Main analysis

  tdwg <- read_sf("manual_downloads/TDWG/old_lv3/level3.shp")

  # Load overall trait completeness
  general_traits_one_percent_threshold <- read_rds("data/focal_trait_coverage.rds")

###########################################################
# Mean coverage map

  tdwg <-
  general_traits_one_percent_threshold %>%
    group_by(area) %>%
    summarise(mean_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))

tdwg %>% st_transform(crs = st_crs(6933))%>%
ggplot()+
  geom_sf(mapping = aes(fill = mean_coverage*100))+
  scale_fill_viridis_b(name = "Mean \ncoverage \n(%)")+
  theme_minimal()


tdwg %>% st_transform(crs = st_crs(6933))%>%
  ggplot()+
  geom_sf(mapping = aes(fill = mean_coverage*100))+
  scale_fill_gradient(low = "white",
                      high = "magenta",
                      name = "Mean \ncoverage\n(%)",
                      limits=c(0,100))+
  theme_minimal()


#######################################################

#Trivariate map

  rudbeck_data <- read.table("data\\variables_2022.csv",sep = ",",
                             header = TRUE)


  tdwg <-
    rudbeck_data %>%
    dplyr::select(LEVEL_3_CO,GEN_DUP,BIEN_OCCUR) %>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"))

  tdwg %>%
    mutate(TRAIT_COVERAGE = mean_coverage*100)->tdwg


#GEN_DUP,BIEN_OCCUR

  library(tricolore)
  tric <- Tricolore(tdwg, p1 = 'GEN_DUP', p2 = 'BIEN_OCCUR', p3 = 'mean_coverage')


  tdwg$rgb <- tric$rgb
  library(ggtern)

  tdwg %>%
    st_transform(crs = st_crs(6933))%>%
  ggplot() +
    # ...draw a polygon for each region...
    geom_sf(aes(fill = rgb, geometry = geometry), size = 0.1) +
    # ...and color each region according to the color code in the variable `rgb`
    scale_fill_identity() +
    annotation_custom(
      ggplotGrob(tric$key +
                   theme(plot.background = element_rect(fill = NA, color = NA)) +
                   labs(L = 'Genetic', T = 'Distribution', R = 'Traits')),
      xmin = -19367530, xmax = -5367530,
      ymin = -7342230 , ymax = 300000
    )+theme_minimal()

###########################

# Georeferenced map

  # Load overall trait completeness
  geo_traits_one_percent_threshold <- read_rds("data/focal_georeferenced_trait_coverage.rds")


  geo_traits_one_percent_threshold


  tdwg <-
    geo_traits_one_percent_threshold %>%
    group_by(area) %>%
    summarise(mean_geo_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))


  tdwg %>% st_transform(crs = st_crs(6933))%>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_geo_coverage*100))+
    scale_fill_viridis_b(name = "Mean \ngeoreferenced\ncoverage \n(%)")+
    theme_minimal()

  tdwg %>% st_transform(crs = st_crs(6933))%>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_geo_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \ngeoreferenced\ncoverage \n(%)",
                        limits=c(0,100))+
    theme_minimal()

  tdwg %>% st_transform(crs = st_crs(6933))%>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_geo_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \ngeoreferenced\ncoverage \n(%)")+
    theme_minimal()


######################################

  #Wood map

  wood_traits_focal_one_percent_threshold <- readRDS("data/focal_wood_trait_coverage.rds")


  wood_tdwg <-
    wood_traits_focal_one_percent_threshold %>%
    group_by(area) %>%
    summarise(mean_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))

  wood_tdwg %>% st_transform(crs = st_crs(6933))%>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \ncoverage\n(%)",
                        limits=c(0,100))+
    theme_minimal()+ggtitle("Wood Traits")

#############################################


  # Endemism map




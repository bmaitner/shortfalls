#' @author Brian Maitner
#' @description Scripts to make maps for the shortfalls paper

# Load libraries
library(tidyverse)
library(sf)
library(lwgeom)
library("ggplot2")
library(cowplot)

###########################################################

# Main analysis

  tdwg <- read_sf("manual_downloads/TDWG/old_lv3/level3.shp")

  # Load overall trait completeness
  general_traits_one_percent_threshold <- read_rds("data/cleaning_raw_names/focal_trait_coverage.rds")

  write.csv(x = general_traits_one_percent_threshold,
            file = "tables/focal_trait_coverage_per_country.csv",
            row.names = FALSE)

###########################################################
  #Summary stats of trait completeness
  min(general_traits_one_percent_threshold$completeness)
  countries_with_no_traits  <- general_traits_one_percent_threshold[which(general_traits_one_percent_threshold$completeness==0),]
  countries_with_all_traits  <- general_traits_one_percent_threshold[which(general_traits_one_percent_threshold$completeness==1),]


  countries_with_no_traits %>%
    inner_join(tdwg %>%
                 dplyr::select(LEVEL_3_CO,LEVEL_NAME)%>%
                 st_drop_geometry(),
               by = c("area" = "LEVEL_3_CO"))-> countries_with_no_traits

    unique(countries_with_no_traits$LEVEL_NAME) #all islands and antarctica


    countries_with_all_traits %>%
      inner_join(tdwg %>%
                   dplyr::select(LEVEL_3_CO,LEVEL_NAME)%>%
                   st_drop_geometry(),
                 by = c("area" = "LEVEL_3_CO"))-> countries_with_all_traits

    unique(countries_with_all_traits$LEVEL_NAME) #all islands and antarctica

    #mean and median completeness within botanical countries
    mean(general_traits_one_percent_threshold$completeness)#19.4
    median(general_traits_one_percent_threshold$completeness)#12.7




###########################################################
# Mean coverage map (unfortunately, the specified projection is a pain)

  tdwg <-
  general_traits_one_percent_threshold %>%
    group_by(area) %>%
    summarise(mean_coverage_focal = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))


  crs_wintri <- "+proj=wintri +datum=WGS84 +no_defs +over"
  world_wintri <- st_transform_proj(tdwg, crs = crs_wintri)

  grat_wintri <-
    st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
    st_transform_proj(crs = crs_wintri)


  ggplot(world_wintri) +
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean Trait\nCompleteness\n(%)",
                        limits=c(0,100))+
    #geom_sf(size = 0.5/.pt) +
    coord_sf(datum = NULL) +
    theme_map()-> overall_trait_completeness


  # summary stats for mean coverage
    min(tdwg$mean_coverage_focal)*100 #2.84%
    tdwg$LEVEL_NAME[which.min(tdwg$mean_coverage_focal)] #New Guinea

    max(tdwg$mean_coverage_focal)*100 #58.67
    tdwg$LEVEL_NAME[which.max(tdwg$mean_coverage_focal)] #"Føroyar" aka Faroe islands

    mean(tdwg$mean_coverage_focal)*100 #19.49
    median(tdwg$mean_coverage_focal)*100 #17.31

  # #############
  #
  # tdwg %>%
  #   st_transform(crs = st_crs(6933))%>%
  #   ggplot()+
  #   geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
  #   scale_fill_gradient(low = "white",
  #                       high = "magenta",
  #                       name = "Mean Trait\nCompleteness\n(%)",
  #                       limits=c(0,100))+
  #   theme_minimal()-> overall_trait_completeness

ggsave(plot = overall_trait_completeness,
       filename = "plots/focal_completeness.svg",width = 10,height = 10)

ggsave(plot = overall_trait_completeness,
       filename = "plots/focal_completeness.jpg",width = 10,height = 4.5, bg = "white")

ggsave(plot = overall_trait_completeness,
       filename = "plots/fig1_high_res.pdf",width = 10,height = 4.5, bg = "white",dpi = 600)


#######################################################

#Load Rudbeck data

  rudbeck_data <- read.table("data\\variables_2022.csv",sep = ",",
                             header = TRUE)

  focal_trait_coverage_no_ferns <- readRDS("data/cleaning_raw_names/focal_trait_coverage_no_ferns.rds")


  tdwg <-
    rudbeck_data %>%
    dplyr::select(LEVEL_3_CO,GEN_DUP,BIEN_OCCUR) %>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"))

  tdwg %>%
    mutate(TRAIT_COVERAGE = mean_coverage_focal*100)->tdwg

  tdwg <-
    focal_trait_coverage_no_ferns %>%
    group_by(area) %>%
    summarise(mean_coverage_focal_no_ferns = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))


  tdwg %>%
    mutate(TRAIT_COMPLETENESS_NO_FERNS = mean_coverage_focal_no_ferns*100) -> tdwg


##########################


  tdwg %>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Trait \nCompleteness\n(%)",
                        limits=c(0,100))+
    theme_minimal()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    coord_sf(datum = NULL) +
    theme_map()





  # 3 shortfalls 0 to 100

trait_completeness <-
  tdwg %>%
    #st_transform(crs = st_crs(6933)) %>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Trait \nCompleteness\n(%)",
                        limits=c(0,100))+
    #theme_minimal()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    coord_sf(datum = NULL) +
    theme_map()


  dist_completeness <-
    tdwg %>%
    #st_transform(crs = st_crs(6933)) %>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = BIEN_OCCUR*100))+
    scale_fill_gradient(low = "white",
                        high = "#74ee15",
                        name = "Distributional \nCompleteness\n(%)",
                        limits=c(0,100))+
    #theme_minimal()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    coord_sf(datum = NULL) +
    theme_map()




  phylo_completeness <-
    tdwg %>%
    #st_transform(crs = st_crs(6933)) %>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = GEN_DUP * 100))+
    scale_fill_gradient(low = "white",
                        high = "#00D1D0",
                        name = "Phylogenetic \nCompleteness\n(%)",
                        limits=c(0,100))+
    #theme_minimal()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    coord_sf(datum = NULL) +
    theme_map()


# 3 shortfalls free

  trait_completeness_free <-
    tdwg %>%
    #st_transform(crs = st_crs(6933)) %>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Trait \nCompleteness\n(%)")+
    #theme_minimal()+
    coord_sf(datum = NULL) +
    theme_map()



  dist_completeness_free <-
    tdwg %>%
    #st_transform(crs = st_crs(6933)) %>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = BIEN_OCCUR*100))+
    scale_fill_gradient(low = "white",
                        high = "#74ee15",
                        name = "Distributional \nCompleteness\n(%)")+
    #theme_minimal()+
    coord_sf(datum = NULL) +
    theme_map()




  phylo_completeness_free <-
    tdwg %>%
    #st_transform(crs = st_crs(6933)) %>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = GEN_DUP * 100))+
    scale_fill_gradient(low = "white",
                        high = "#00D1D0",
                        name = "Phylogenetic \nCompleteness\n(%)")+
    #theme_minimal()+
    coord_sf(datum = NULL) +
    theme_map()



library(ggpubr)


  ggarrange(trait_completeness,
            phylo_completeness,
            dist_completeness,ncol = 1)

  ggarrange(trait_completeness_free,
            phylo_completeness_free,
            dist_completeness_free,
            ncol = 1)

  ggarrange(phylo_completeness_free,
            dist_completeness_free,
            ncol = 1)

  trait_completeness_free

  # dist, phylo completeness + scatterplots with trait v phylo, trait vs dist.  (maybe colored by a ramp from one to the other?)

  tdwg %>%
    mutate(GENETIC_COMPLETENESS = GEN_DUP * 100,
           DISRIBUTIONAL_COMPLETENESS = BIEN_OCCUR *100)->tdwg


  library(ggpmisc)

  trait_v_gene <-
  ggplot(data = tdwg,mapping = aes(x= TRAIT_COMPLETENESS_NO_FERNS,
                                   y = GENETIC_COMPLETENESS,
                                   ))+
                                   #color = GENETIC_COMPLETENESS -TRAIT_COMPLETENESS_NO_FERNS))+
    geom_point()+
    #xlim(c(0,100))+
    #ylim(c(0,100))+
    # scale_color_gradient(low = "magenta",
    #                      high = "#00D1D0",
    #                      limits=c(-100,100),name= "Bias")+
    geom_abline(slope = 1,intercept = 0,lty=2)+
    #stat_smooth(method = "lm")+
    xlab("Trait Completeness (%)")+
    ylab("Phylogenetic Completeness (%)")+
    stat_correlation(method = "pearson",aes(label = paste(after_stat(r.label),
                                                          after_stat(p.value.label),
                                                          after_stat(n.label),
                                                          sep = "*\"; \"*")),
                     small.p = TRUE,
                     small.r = TRUE)+
    theme_bw()

  trait_v_dist <-

  ggplot(data = tdwg,mapping = aes(x= TRAIT_COMPLETENESS_NO_FERNS,
                                   y = DISRIBUTIONAL_COMPLETENESS
                                   ))+
                                   # ,color = DISRIBUTIONAL_COMPLETENESS -TRAIT_COMPLETENESS_NO_FERNS))+

    geom_point()+
    #xlim(c(0,100))+
    #ylim(c(0,100))+
    # scale_color_gradient(low = "magenta",
    #                      high = "#74ee15",
    #                      limits=c(-100,100),
    #                      name= "Bias")+
    geom_abline(slope = 1,intercept = 0,lty=2)+
    #stat_smooth(method = "lm",show.legend = TRUE)+
    xlab("Trait Completeness (%)")+
    ylab("Distributional Completeness (%)")+
#    stat_poly_line() +
 #   stat_poly_eq(aes(label = paste(after_stat(eq.label),
  #                                 after_stat(rr.label), sep = "*\", \"*")))+

    stat_correlation(method = "pearson",aes(label = paste(after_stat(r.label),
                                                          after_stat(p.value.label),
                                                          after_stat(n.label),
                                                          sep = "*\"; \"*")),
                     small.p = TRUE,
                     small.r = TRUE)+
    theme_bw()

  library(ggplot2)

  phylo_and_dist <-
  ggarrange(phylo_completeness_free,trait_v_gene,
            dist_completeness_free,trait_v_dist,
            ncol = 2,
            nrow = 2,
            widths = c(2,1),
            labels = "AUTO")

  ggsave(plot = phylo_and_dist, filename = "plots/phylo_and_dist_completeness.svg",width = 10,height = 4.5)
  ggsave(plot = phylo_and_dist,filename = "plots/phylo_and_dist_completeness_tall.svg",width = 10,height = 4.5)
  ggsave(plot = phylo_and_dist,filename = "plots/phylo_and_dist_completeness.jpg",width = 10,height = 4.5, bg = "white")
  ggsave(plot = phylo_and_dist,filename = "plots/phylo_and_dist_completeness_tall.jpg",width = 10,height = 10, bg = "white")
  ggsave(plot = phylo_and_dist,filename = "plots/phylo_and_dist_completeness_med.jpg",width = 10,height = 7, bg = "white")

  ggsave(plot = phylo_and_dist,filename = "plots/fig2_high_res.pdf",width = 10,height = 7, bg = "white",dpi = 600)



  # Distrbutional data
    max(tdwg$DISRIBUTIONAL_COMPLETENESS) #100%
    min(tdwg$DISRIBUTIONAL_COMPLETENESS) #0%
    mean(tdwg$DISRIBUTIONAL_COMPLETENESS) #47.12%
    median(tdwg$DISRIBUTIONAL_COMPLETENESS) #51.66%


  # Phylogenetic data
    max(tdwg$GENETIC_COMPLETENESS) #11.28%
    min(tdwg$GENETIC_COMPLETENESS) #1.34%
    mean(tdwg$GENETIC_COMPLETENESS) #5.67%
    median(tdwg$GENETIC_COMPLETENESS) #5.40%

  # Traits data
    max(tdwg$TRAIT_COVERAGE) #58.68%
    min(tdwg$TRAIT_COVERAGE) #2.84%
    mean(tdwg$TRAIT_COVERAGE) #19.49%
    median(tdwg$TRAIT_COVERAGE) #17.31%


  #Where do we have more trait than distributional data?

    tdwg %>%
      mutate(dist_v_trait = DISRIBUTIONAL_COMPLETENESS - TRAIT_COMPLETENESS_NO_FERNS)%>%
      ggplot(mapping = aes(fill = dist_v_trait))+
      geom_sf()+
      scale_fill_gradient2(low = "magenta",mid="white",high = "#74ee15") #looks like Russia and parts of northern africa may be trait biased

    tdwg %>%
      mutate(dist_v_trait = DISRIBUTIONAL_COMPLETENESS - TRAIT_COMPLETENESS_NO_FERNS)%>%
      st_drop_geometry()->tdwg_info


###########################

# Georeferenced map

  # Load overall trait completeness
  geo_traits_one_percent_threshold <- read_rds("data/cleaning_raw_names/focal_georeferenced_trait_coverage.rds")


  geo_traits_one_percent_threshold


  tdwg <-
    geo_traits_one_percent_threshold %>%
    group_by(area) %>%
    summarise(mean_geo_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))


  georef_plot <-
  tdwg %>%
    #st_transform(crs = st_crs(6933))%>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = mean_geo_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \nGeoreferenced\nTrait\nCompleteness \n(%)")+
    #theme_minimal()
      coord_sf(datum = NULL) +
      theme_map()


  geo_v_overall <-
    ggplot(data = tdwg,
           mapping = aes(x= mean_coverage_focal*100,
                                     y = mean_geo_coverage*100
    ))+
    # ,color = DISRIBUTIONAL_COMPLETENESS -TRAIT_COMPLETENESS_NO_FERNS))+

    geom_point()+
    geom_abline(slope = 1,intercept = 0,lty=2)+
    xlab("Focal Trait Completeness (%)")+
    ylab("Georeferenced Trait Completeness (%)")+
    stat_correlation(method = "pearson",aes(label = paste(after_stat(r.label),
                                                          after_stat(p.value.label),
                                                          after_stat(n.label),
                                                          sep = "*\"; \"*")),
                     small.p = TRUE,
                     small.r = TRUE)+
    theme_bw()

  #combine map and plot
  geo_plot_and_map <-
  ggarrange(georef_plot, geo_v_overall,
            ncol = 2,
            widths = c(2,1),
            labels = "AUTO")

  ggsave(plot = geo_plot_and_map,
         filename = "plots/geo_completeness_and_plot.jpg",width = 10,height = 3.6)
  ggsave(plot = geo_plot_and_map,
         filename = "plots/geo_completeness_and_plot.svg",width = 10,height = 3.6)



  # merge with predictor variables
  socio_vars <- read.csv("manual_downloads/Darwinian_shortfalls/socioeco_var.csv")

  # combine trait completeness and predictor variables
  general_traits_one_percent_threshold <-
    merge(x = general_traits_one_percent_threshold,
          y = socio_vars,
          by.x= "area",
          by.y = "LEVEL_3_CO")

  tdwg %>%
    full_join(y = socio_vars) -> tdwg


######################################

  #Wood map

  wood_traits_focal_one_percent_threshold <- readRDS("data/cleaning_raw_names/focal_wood_trait_coverage.rds")


  tdwg <-
    wood_traits_focal_one_percent_threshold %>%
    group_by(area) %>%
    summarise(wood_mean_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))



tdwg %>%
  #st_transform(crs = st_crs(6933))%>%
  st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = wood_mean_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \ncoverage\n(%)",
                        limits=c(0,100))+
    #theme_minimal()+
    ggtitle("Wood Traits")+
    coord_sf(datum = NULL) +
    theme_map()


#############################################

  # flower map

  flower_traits_focal_one_percent_threshold <- readRDS("data/cleaning_raw_names/focal_flower_trait_coverage.rds")


  tdwg <-
    flower_traits_focal_one_percent_threshold %>%
    group_by(area) %>%
    summarise(flower_mean_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))


  tdwg %>%
    #st_transform(crs = st_crs(6933))%>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = flower_mean_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \ncoverage\n(%)",
                        limits=c(0,100))+
    #theme_minimal()+
    ggtitle("Flower Traits")+
    coord_sf(datum = NULL) +
    theme_map()


  #############################################


  # seed map

  seed_traits_focal_one_percent_threshold <- readRDS("data/cleaning_raw_names/focal_seed_trait_coverage.rds")


  tdwg <-
    seed_traits_focal_one_percent_threshold %>%
    group_by(area) %>%
    summarise(seed_mean_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))

  tdwg %>%
    #st_transform(crs = st_crs(6933))%>%
    st_transform_proj(crs = crs_wintri) %>%
    ggplot()+
    geom_sf(data = grat_wintri, color = "gray30", size = 0.25/.pt,alpha=.5) +
    geom_sf(mapping = aes(fill = seed_mean_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \ncoverage\n(%)",
                        limits=c(0,100))+
    #theme_minimal()+
    ggtitle("Seed Traits")+
    coord_sf(datum = NULL) +
    theme_map()

  ##############################################




  library(ggpubr)

  trait_subsets_fixed <-
  ggarrange(tdwg %>%
              # st_transform(crs = st_crs(6933))%>%
              st_transform_proj(crs = crs_wintri) %>%
              ggplot()+
              geom_sf(data = grat_wintri,
                      color = "gray30",
                      size = 0.25/.pt,
                      alpha=.5) +
              geom_sf(mapping = aes(fill = wood_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)",
                                  limits=c(0,100))+
              #theme_minimal()+
              ggtitle("A. Wood Traits")+
              coord_sf(datum = NULL) +
              theme_map()

            ,
            tdwg %>%
              #st_transform(crs = st_crs(6933))%>%
              st_transform_proj(crs = crs_wintri) %>%
              ggplot()+
              geom_sf(data = grat_wintri,
                      color = "gray30",
                      size = 0.25/.pt,
                      alpha=.5) +
              geom_sf(mapping = aes(fill = flower_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)",
                                  limits=c(0,100))+
              #theme_minimal()+
              ggtitle("B. Flower Traits")+
              coord_sf(datum = NULL) +
              theme_map()
            ,

            tdwg %>%
              #st_transform(crs = st_crs(6933))%>%
              st_transform_proj(crs = crs_wintri) %>%
              ggplot()+
              geom_sf(data = grat_wintri,
                      color = "gray30",
                      size = 0.25/.pt,
                      alpha=.5) +
              geom_sf(mapping = aes(fill = seed_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)",
                                  limits=c(0,100))+
              #theme_minimal()+
              ggtitle("C. Seed Traits")+
              coord_sf(datum = NULL) +
              theme_map(),
            ncol = 1,
            common.legend = TRUE,
            legend = "right")


  trait_subsets_free <-
  ggarrange(tdwg %>%
              #st_transform(crs = st_crs(6933))%>%
              st_transform_proj(crs = crs_wintri) %>%
              ggplot()+
              geom_sf(data = grat_wintri,
                      color = "gray30",
                      size = 0.25/.pt,
                      alpha=.5) +
              geom_sf(mapping = aes(fill = wood_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)")+
              #theme_minimal()+
              ggtitle("A. Wood Traits")+
              coord_sf(datum = NULL) +
              theme_map()
            ,
            tdwg %>%
              #st_transform(crs = st_crs(6933))%>%
              st_transform_proj(crs = crs_wintri) %>%
              ggplot()+
              geom_sf(data = grat_wintri,
                      color = "gray30",
                      size = 0.25/.pt,
                      alpha=.5) +
              geom_sf(mapping = aes(fill = flower_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)")+
              #theme_minimal()+
              ggtitle("B. Flower Traits")+
              coord_sf(datum = NULL) +
              theme_map()
            ,

            tdwg %>%
              #st_transform(crs = st_crs(6933))%>%
              st_transform_proj(crs = crs_wintri) %>%
              ggplot()+
              geom_sf(data = grat_wintri,
                      color = "gray30",
                      size = 0.25/.pt,
                      alpha=.5) +
              geom_sf(mapping = aes(fill = seed_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)")+
              #theme_minimal()+
              ggtitle("C. Seed Traits")+
              coord_sf(datum = NULL) +
              theme_map(),
            ncol = 1,
            legend = "right")

    ggsave(plot = trait_subsets_fixed,
           filename = "plots/trait_subsets.jpg",width = 10,height = 4.5*3)
    ggsave(plot = trait_subsets_fixed,
           filename = "plots/trait_subsets.svg",width = 10,height = 4.5*3)


  ##############################################

  #correlations

  tdwg_info <- st_drop_geometry(tdwg)

  tri <- tdwg_info %>%
    dplyr::select(GEN_DUP,BIEN_OCCUR,mean_coverage_focal)


  library(corrplot)

  data_for_cor<-
  tdwg_info %>%
    dplyr::select(mean_coverage_focal,
           mean_geo_coverage,
           wood_mean_coverage,
           flower_mean_coverage,
           seed_mean_coverage,
           GEN_DUP,
           BIEN_OCCUR)



  correlations <- cor(data_for_cor,
                      method = "pearson")

  colnames(correlations)<- c("Focal Trait Completeness",
                             "Georeferenced Trait Completeness",
                             "Wood Trait Completeness",
                             "Flower Trait Completeness",
                             "Seed Trait Completeness",
                             "Phylogenetic Completeness",
                             "Distributional Completeness")

  rownames(correlations)<- c("Focal Trait Completeness",
                             "Georeferenced Trait Completeness",
                             "Wood Trait Completeness",
                             "Flower Trait Completeness",
                             "Seed Trait Completeness",
                             "Phylogenetic Completeness",
                             "Distributional Completeness")


  correlation_pvals <- cor.mtest(data_for_cor,
                                 conf.level = 0.95,
                                 method="pearson",alternative="two.sided",exact = TRUE)


  library(Cairo)

  Cairo(file="plots/data_type_correlations.jpg",
        type="jpg",
        units="in",
        width=10,
        height=10,
        dpi=300)


  corrplot( correlations ,
            addCoef.col="black",
            order = "FPC",
            tl.col = "black",
            tl.cex = c(rep(.75,6),1),
            p.mat = correlation_pvals$p)


  ## When the device is off, file writing is completed.
  dev.off()




#####################################################################


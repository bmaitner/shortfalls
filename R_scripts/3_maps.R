#' @author Brian Maitner
#' @description Scripts to make maps for the shortfalls paper

# Load libraries
library(tidyverse)
library(lme4)
library(sf)
library(tricolore)
library(expss)

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
    summarise(mean_coverage_focal = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))

tdwg %>% st_transform(crs = st_crs(6933))%>%
ggplot()+
  geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
  scale_fill_viridis_b(name = "Mean \ncoverage \n(%)")+
  theme_minimal()


tdwg %>% st_transform(crs = st_crs(6933))%>%
  ggplot()+
  geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
  scale_fill_gradient(low = "white",
                      high = "magenta",
                      name = "Mean \nCompleteness\n(%)",
                      limits=c(0,100))+
  theme_minimal()

#######################################################

#Trivariate map

  rudbeck_data <- read.table("data\\variables_2022.csv",sep = ",",
                             header = TRUE)

  focal_trait_coverage_no_ferns <- readRDS("data/focal_trait_coverage_no_ferns.rds")


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
    mutate(TRAIT_COMPLETENESS_NO_FERNS = mean_coverage_focal_no_ferns*100)->tdwg


#GEN_DUP,BIEN_OCCUR

  library(tricolore)
  tric <- Tricolore(tdwg, p1 = 'GEN_DUP', p2 = 'BIEN_OCCUR', p3 = 'mean_coverage_focal')


  tric$key$data

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





##########################

# # two sets of choropleth maps
#
#
# library(biscale)
# library(cowplot)
#
#
# bi_class<-  bi_class(.data = data.frame(x=1:100,y=100:1),x = x,y = y,
#            style = "equal",
#            dim = 4,
#            keep_factors = TRUE,
#            dig_lab = 2)
#
#
# tdwg %>%
#   ggplot() +
#   geom_sf( mapping = aes(fill = bi_class), color = "black", size = 0.1, show.legend = FALSE) +
#   bi_scale_fill(pal = "GrPink2", dim = 4) +
#   labs(
#     title = "Race and Income in St. Louis, MO",
#     subtitle = "Gray Pink (GrPink) Palette"
#   )+
#   bi_theme()->test_map
#
#
# legend <- bi_legend(pal = "GrPink2",
#                     dim = 4,
#                     xlab = "Traits",
#                     ylab = "Occurrences",
#                     size = 8)
#
# finalPlot <- ggdraw() +
#   draw_plot(test_map, 0, 0, 1, 1) +
#   draw_plot(legend, 0.2, .65, 0.2, 0.2)
#
# finalPlot
# ?bi_class
#
#


##########################

  # 3 shortfalls 0 to 100

trait_completeness <-
  tdwg %>%
    st_transform(crs = st_crs(6933)) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Trait \nCompleteness\n(%)",
                        limits=c(0,100))+
    theme_minimal()


  dist_completeness <-
    tdwg %>%
    st_transform(crs = st_crs(6933)) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = BIEN_OCCUR*100))+
    scale_fill_gradient(low = "white",
                        high = "#74ee15",
                        name = "Distributional \nCompleteness\n(%)",
                        limits=c(0,100))+
    theme_minimal()



  phylo_completeness <-
    tdwg %>%
    st_transform(crs = st_crs(6933)) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = GEN_DUP * 100))+
    scale_fill_gradient(low = "white",
                        high = "#00D1D0",
                        name = "Phylogenetic \nCompleteness\n(%)",
                        limits=c(0,100))+
    theme_minimal()

# 3 shortfalls free

  trait_completeness_free <-
    tdwg %>%
    st_transform(crs = st_crs(6933)) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = mean_coverage_focal*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Trait \nCompleteness\n(%)")+
    theme_minimal()


  dist_completeness_free <-
    tdwg %>%
    st_transform(crs = st_crs(6933)) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = BIEN_OCCUR*100))+
    scale_fill_gradient(low = "white",
                        high = "#74ee15",
                        name = "Distributional \nCompleteness\n(%)")+
    theme_minimal()



  phylo_completeness_free <-
    tdwg %>%
    st_transform(crs = st_crs(6933)) %>%
    ggplot()+
    geom_sf(mapping = aes(fill = GEN_DUP * 100))+
    scale_fill_gradient(low = "white",
                        high = "#00D1D0",
                        name = "Phylogenetic \nCompleteness\n(%)")+
    theme_minimal()


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

  tdwg$TRAIT_COMPLETENESS_NO_FERNS

  tdwg %>%
    mutate(GENETIC_COMPLETENESS = GEN_DUP * 100,
           DISRIBUTIONAL_COMPLETENESS = BIEN_OCCUR *100)->tdwg

  tdwg %>%
    mutate(gen_v_trait = GENETIC_COMPLETENESS/TRAIT_COMPLETENESS_NO_FERNS,
           dist_v_trait = DISRIBUTIONAL_COMPLETENESS / TRAIT_COMPLETENESS_NO_FERNS)->tdwg

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
                     small.r = TRUE)

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
                     small.r = TRUE)

  library(ggplot2)

  ggarrange(phylo_completeness_free,trait_v_gene,
            dist_completeness_free,trait_v_dist,
            ncol = 2,
            nrow = 2,
            widths = c(2,1),
            labels = "AUTO")



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

  wood_traits_focal_one_percent_threshold <- readRDS("data/focal_wood_trait_coverage.rds")


  tdwg <-
    wood_traits_focal_one_percent_threshold %>%
    group_by(area) %>%
    summarise(wood_mean_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))

tdwg %>% st_transform(crs = st_crs(6933))%>%
    ggplot()+
    geom_sf(mapping = aes(fill = wood_mean_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \ncoverage\n(%)",
                        limits=c(0,100))+
    theme_minimal()+ggtitle("Wood Traits")

#############################################

  # flower map

  flower_traits_focal_one_percent_threshold <- readRDS("data/focal_flower_trait_coverage.rds")


  tdwg <-
    flower_traits_focal_one_percent_threshold %>%
    group_by(area) %>%
    summarise(flower_mean_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))

  tdwg %>% st_transform(crs = st_crs(6933))%>%
    ggplot()+
    geom_sf(mapping = aes(fill = flower_mean_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \ncoverage\n(%)",
                        limits=c(0,100))+
    theme_minimal()+ggtitle("Flower Traits")

  #############################################


  # seed map

  seed_traits_focal_one_percent_threshold <- readRDS("data/focal_seed_trait_coverage.rds")


  tdwg <-
    seed_traits_focal_one_percent_threshold %>%
    group_by(area) %>%
    summarise(seed_mean_coverage = mean (completeness))%>%
    right_join(x = tdwg,
               y = .,
               by = c("LEVEL_3_CO"="area"))

  tdwg %>% st_transform(crs = st_crs(6933))%>%
    ggplot()+
    geom_sf(mapping = aes(fill = seed_mean_coverage*100))+
    scale_fill_gradient(low = "white",
                        high = "magenta",
                        name = "Mean \ncoverage\n(%)",
                        limits=c(0,100))+
    theme_minimal()+ggtitle("Seed Traits")

  ##############################################




  library(ggpubr)


  ggarrange(tdwg %>% st_transform(crs = st_crs(6933))%>%
              ggplot()+
              geom_sf(mapping = aes(fill = wood_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)",
                                  limits=c(0,100))+
              theme_minimal()+ggtitle("A. Wood Traits")
            ,
            tdwg %>% st_transform(crs = st_crs(6933))%>%
              ggplot()+
              geom_sf(mapping = aes(fill = flower_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)",
                                  limits=c(0,100))+
              theme_minimal()+ggtitle("B. Flower Traits")
            ,

            tdwg %>% st_transform(crs = st_crs(6933))%>%
              ggplot()+
              geom_sf(mapping = aes(fill = seed_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)",
                                  limits=c(0,100))+
              theme_minimal()+ggtitle("C. Seed Traits"),
            ncol = 1,
            common.legend = TRUE,
            legend = "right")


  ggarrange(tdwg %>% st_transform(crs = st_crs(6933))%>%
              ggplot()+
              geom_sf(mapping = aes(fill = wood_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)")+
              theme_minimal()+ggtitle("A. Wood Traits")
            ,
            tdwg %>% st_transform(crs = st_crs(6933))%>%
              ggplot()+
              geom_sf(mapping = aes(fill = flower_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)")+
              theme_minimal()+ggtitle("B. Flower Traits")
            ,

            tdwg %>% st_transform(crs = st_crs(6933))%>%
              ggplot()+
              geom_sf(mapping = aes(fill = seed_mean_coverage*100))+
              scale_fill_gradient(low = "white",
                                  high = "magenta",
                                  name = "Mean \nCompleteness\n(%)")+
              theme_minimal()+ggtitle("C. Seed Traits"),
            ncol = 1,
            legend = "right")









  ##############################################


  # Endemism map


  ##############################################

  #correlations

  tdwg_info <- st_drop_geometry(tdwg)

  tri<-tdwg_info%>%
    dplyr::select(GEN_DUP,BIEN_OCCUR,mean_coverage_focal)


  library(corrplot)
  corrplot()
?corrplot


  data_for_cor<-
  tdwg_info %>%
    dplyr::select(mean_coverage_focal,
           mean_geo_coverage,
           wood_mean_coverage,
           flower_mean_coverage,
           seed_mean_coverage,
           GEN_DUP,
           BIEN_OCCUR)



  correlations <- cor(data_for_cor,method = "pearson")
  colnames(correlations)<- c("General Traits","Georeferenced Traits",
                             "Wood Traits","Flower Traits",
                             "Seed Traits","Genes","Occurrences")

  rownames(correlations)<- c("General Traits","Georeferenced Traits",
                             "Wood Traits","Flower Traits",
                             "Seed Traits","Genes","Occurrences")


  correlation_pvals <- cor.mtest(data_for_cor,
                                 conf.level = 0.95,
                                 method="pearson",alternative="two.sided",exact = TRUE)


  corrplot(  cor(data_for_cor),
             addCoef.col="black")


  corrplot( correlations ,
             addCoef.col="black",
             order = "FPC",
             tl.col = "black",
            p.mat = correlation_pvals$p,
            lowCI.mat = correlation_pvals$lowCI,
            uppCI.mat = correlation_pvals$uppCI,plotCI = "circle")


  corrplot( correlations ,
            addCoef.col="black",
            order = "FPC",
            tl.col = "black",
            p.mat = correlation_pvals$p)

#####################################################################



library(sf)


  hist(preds$mean_diff_from_pred)

plot(preds$exp_mean_completeness,
     preds$mean_completeness)


    preds %>%
      st_transform(crs = st_crs(6933))%>%
      ggplot()+
      geom_sf(mapping = aes(fill = (mean_completeness - exp_mean_completeness)*100))+
      scale_fill_gradient(low = "white",
                          high = "magenta",
                          name = "Mean \ncoverage\n(%)")+
      theme_minimal()+ggtitle("Expect % - Obs. %")






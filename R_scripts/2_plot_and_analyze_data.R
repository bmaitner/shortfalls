#' @author Brian Maitner
#' @description Scripts to wrangle data for the trait shortfalls manuscript

# Load libraries
library(tidyverse)
library(lme4)

# Get TRY data

library(sf)
#tdwg <- read_sf("manual_downloads/TDWG/level3.shp")
tdwg <- read_sf("manual_downloads/TDWG/old_lv3/level3.shp")

#Load in the trait data
  #traits <- arrow::open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

#examine the data structure
  # traits$schema
  # traits$metadata





########################################################

# Trait completeness model: 1% threshold and broadly applicable traits
  #via a single model

  # Load overall trait completeness
    # general_traits_one_percent_threshold <- read_rds("data/focal_trait_coverage.rds")
    general_traits_one_percent_threshold <- read_rds("data/cleaning_raw_names/focal_trait_coverage.rds")

  # merge with predictor variables
    socio_vars <- read.csv("manual_downloads/Darwinian_shortfalls/socioeco_var.csv")
    #socio_vars[2:10] <- scale(socio_vars[2:10])

  # combine trait completeness and predictor variables
    general_traits_one_percent_threshold <-
    merge(x = general_traits_one_percent_threshold,
          y = socio_vars,
          by.x= "area",
          by.y = "LEVEL_3_CO")

  # add species richness
    wcvp <- readRDS("manual_downloads/WCVP/wcvp_cleaned.RDS")

    wcvp %>%
        group_by(area_code_l3) %>%
        summarise(richness = n()) %>%
      merge(x = general_traits_one_percent_threshold,
            y = .,
            by.x = "area",
            by.y = "area_code_l3",
            all.x = TRUE) ->
      general_traits_one_percent_threshold

    # Add mean range size

      socio_vars %>%
        dplyr::select(LEVEL_3_CO, AREA_SQKM) %>%
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
        summarise(mean_species_range = mean(na.omit(range_size))) %>%
        merge(x = general_traits_one_percent_threshold,
              y = .,
              by.x="area",
              by.y = "area_code_l3") -> general_traits_one_percent_threshold

      # Add endemism

        wcvp %>%
          group_by(taxon_name) %>%
          summarize(n_countries = n()) %>%
          merge(x = wcvp,
                y = .) %>%
          group_by(area_code_l3) %>%
          filter(n_countries == 1) %>%
          summarise(endemics = n()) %>%
          merge(x = general_traits_one_percent_threshold,
              y = .,
              by.x = "area",
              by.y = "area_code_l3",
              all.x = TRUE) %>%
          mutate(endemics = replace_na(endemics,0)) %>%
          mutate(endemism = endemics/richness) %>%
          dplyr::select(-endemics) -> general_traits_one_percent_threshold

  #rescale predictors

    general_traits_one_percent_threshold[4:ncol(general_traits_one_percent_threshold)] <- scale(general_traits_one_percent_threshold[4:ncol(general_traits_one_percent_threshold)])
    general_traits_one_percent_threshold <- na.omit(general_traits_one_percent_threshold)

  #Simple model

    all_variables <-
      glmer(data = general_traits_one_percent_threshold,
           formula = completeness ~ AREA_SQKM + GDP_SUM + GDP_CAPITA + ROAD_DENSITY + POP_COUNT + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP
           + richness + mean_species_range + endemism +
             (1|trait),
           family = "binomial",
           control = glmerControl(optimizer="bobyqa",
                                  optCtrl=list(maxfun=2e5)),
           na.action = "na.fail")

    plot(y = fitted.values(all_variables),
         x = general_traits_one_percent_threshold$completeness)

    plot(y = resid(all_variables),
         x = general_traits_one_percent_threshold$completeness)

    library(MuMIn)
    summary(all_variables)
    MuMIn::r.squaredGLMM(all_variables)
    r2glmm::r2beta(model = all_variables)

    #check residuals

    resids <- residuals(all_variables)


    # dredged_all_variables <-
    # dredge(global.model = all_variables)
    # dredged_all_variables


    #Moran's I
    tdwg %>%
      inner_join(y = general_traits_one_percent_threshold,
                 by = c("LEVEL_3_CO" = "area")) -> tdwg_combined

    library(spdep)
    library(tidyverse)
    library(gridExtra)
    library(NLMR)
    library(DHARMa)


    #Trying to do them all at once (if this fails, iterate and record)

      nb <- poly2nb(pl = st_make_valid(tdwg_combined),
                    queen = TRUE)

      # saveRDS(object = nb,file = "data/polygon_nb.RDS")

      lw <- nb2listw(nb, style="W", zero.policy=TRUE)

      comp.lag <- lag.listw(lw, tdwg_combined$completeness)

      plot(comp.lag ~ tdwg_combined$completeness, pch=16, asp=1)
      M1 <- lm(comp.lag ~ tdwg_combined$completeness)
      abline(M1, col="blue")

      I <- moran(tdwg_combined$completeness, lw, length(nb), Szero(lw))[1]

      moran_test_combined <- moran.test(tdwg_combined$completeness,
                 lw,
                 alternative = "greater")


      #check whether the autocorrelation is present in the residuals


      moran_test_residuals <- moran.test(residuals(all_variables),
                                        lw,
                                        alternative = "two.sided")


      # Moran I statistic standard deviate = 450.37, p-value < 2.2e-16
      # alternative hypothesis: two.sided
      # sample estimates:
      #   Moran I statistic       Expectation          Variance
      # 4.340564e-01     -5.127416e-05      9.290920e-07


      MC_test_combined <- moran.mc(residuals(all_variables),
                    lw,
                    nsim = 999,
                    alternative = "two.sided")

      # Monte-Carlo simulation of Moran I
      #
      # data:  residuals(all_variables)
      # weights: lw
      # number of simulations + 1: 1000
      #
      # statistic = 0.43406, observed rank = 1000, p-value < 2.2e-16
      # alternative hypothesis: two.sided


  moran_test_combined
  MC_test_combined
  MC_test_combined_two_sided

    plot(MC_test_combined)

    # tdwg_dist <- st_distance(x = st_make_valid(tdwg_combined))
    # saveRDS(object = tdwg_dist,file = "data/tdwg_dist.RDS")
    tdwg_dist <- readRDS(file = "data/tdwg_dist.RDS")

    sims <- simulateResiduals(all_variables)
    testOutliers(simulationOutput = sims,type = "bootstrap")
    testSpatialAutocorrelation(simulationOutput = sims,
                               distMat = tdwg_dist)

library(sf)
library(tidyverse)

  tdwg_cents <- st_centroid(st_make_valid(tdwg))
  tdwg_cents <- cbind(tdwg_cents,st_coordinates(tdwg_cents))

  general_traits_one_percent_threshold %>%
    inner_join(y = tdwg_cents,by = c("area"= "LEVEL_3_CO")) -> general_traits_w_coords


################################

  #spamm to account for spatial autocorrelation

  wcvp %>%
    group_by(area_code_l3) %>%
    summarise(richness_untf = n()) %>%
    merge(x = general_traits_w_coords,
          y = .,
          by.x = "area",
          by.y = "area_code_l3",
          all.x = TRUE) ->
    general_traits_w_coords



  #trying spamm
    # https://datascienceplus.com/spatial-regression-in-r-part-1-spamm-vs-glmmtmb/
    # https://kimura.univ-montp2.fr/~rousset/spaMM/spaMMintro.pdf
  library(spaMM)

  m_spamm <-
  general_traits_w_coords %>%
    mutate(species_w_data = completeness * richness_untf) %>%
    mutate(species_wo_data = (1-completeness) * richness_untf) %>%
  fitme(cbind(species_w_data,species_wo_data) ~ AREA_SQKM + GDP_SUM + GDP_CAPITA + ROAD_DENSITY + POP_COUNT + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP
                   + richness + mean_species_range + endemism +
                     (1|trait) + Matern(1 | X + Y),
                   data = .,
                   family = "binomial") # for spamm, binomial models need to include successes and failures (e.g. samples and no samples)


  m_spamm_null <-
    general_traits_w_coords %>%
    mutate(species_w_data = completeness * richness_untf) %>%
    mutate(species_wo_data = (1-completeness) * richness_untf) %>%
    fitme(cbind(species_w_data,species_wo_data) ~  (1|trait) + Matern(1 | X + Y),
          data = .,
          family = "binomial") # for spamm, binomial models need to include successes and failures (e.g. samples and no samples)


  #Test whether model is significant improvement over null
    spaMM::LRT(object = m_spamm,
               object2 = m_spamm_null) # p << 0.001

  #Look at model results
    summary(m_spamm)


    # ------------ Fixed effects (beta) ------------
    #   Estimate Cond. SE  t-value
    # (Intercept)        -1.827439 0.167811 -10.8898
    # AREA_SQKM          -0.009068 0.015193  -0.5969
    # GDP_SUM            -0.006973 0.008457  -0.8246
    # GDP_CAPITA         -0.010274 0.008403  -1.2226
    # ROAD_DENSITY        0.007441 0.007724   0.9634
    # POP_COUNT           0.004552 0.010248   0.4442
    # POP_DENSITY         0.019438 0.010070   1.9303
    # SECURITY           -0.014420 0.010232  -1.4093
    # RESEARCH_EXP        0.060051 0.014366   4.1799
    # EDUCATION_EXP      -0.005070 0.007721  -0.6567
    # richness           -0.089983 0.014143  -6.3624
    # mean_species_range  0.408990 0.017557  23.2944
    # endemism           -0.125873 0.012587 -10.0005

    sims <- DHARMa::simulateResiduals(fittedModel = m_spamm,integerResponse = TRUE,n = 1000)
    ?simulateResiduals
    plot(sims)
    DHARMa::testOutliers(simulationOutput = sims,
                 type = "bootstrap")

  #Check confidence intervals
    confint(object = m_spamm,
                   parm = c("(Intercept)","AREA_SQKM","GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
                            "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
                            "EDUCATION_EXP",
                            "richness","mean_species_range","endemism"))

    #CIs

    # lower (Intercept) upper (Intercept)
    # -2.201892         -1.463093
    # lower AREA_SQKM upper AREA_SQKM
    # -0.03921223      0.02116688
    # lower GDP_SUM upper GDP_SUM
    # -0.023619696   0.009742648
    # lower GDP_CAPITA upper GDP_CAPITA
    # -0.026964859      0.006244498
    # lower ROAD_DENSITY upper ROAD_DENSITY
    # -0.007744537        0.022702515
    # lower POP_COUNT upper POP_COUNT
    # -0.01559710      0.02469701
    # lower POP_DENSITY upper POP_DENSITY
    # -0.0003513056      0.0393199362
    # lower SECURITY upper SECURITY
    # -0.034531197    0.005745719
    # lower RESEARCH_EXP upper RESEARCH_EXP
    # 0.03165934         0.08863481
    # lower EDUCATION_EXP upper EDUCATION_EXP
    # -0.02027907          0.01011657
    # lower richness upper richness
    # -0.11802884    -0.06189039
    # lower mean_species_range upper mean_species_range
    # 0.3725413                0.4456265
    # lower endemism upper endemism
    # -0.1514902     -0.1002516

  #Also check AICs
    AIC(m_spamm)
    AIC(m_spamm_null)


  plot(m_spamm) #documentation suggests this may be incorrect for binomials and recommends use of dharma package

  # Partial Dependence Plots
    plot_effects(m_spamm,"mean_species_range")
    plot_effects(m_spamm,"AREA_SQKM")


#################################

  #Taxonomic Biases

  family_trait_coverage <- readRDS("data/focal_trait_coverage_family.rds")
  trait_list <- readRDS("data/trait_list_w_coverage.RDS")

  #need overall completeness for comparison
  source("R/test_family_coverage.R")
  family_coverage_tests <- test_family_coverage(family_trait_coverage = family_trait_coverage,
                                                trait_list = trait_list)

  family_coverage_tests %>%
    mutate(significant = case_when(p_value <= 0.05 ~ 1,
                                   p_value > 0.05 ~ 0)) %>%
    mutate(sig_greater = case_when(significant == 1 & exp_completeness > completeness ~ 1,
                                   significant == 1 & exp_completeness < completeness ~ 0,
                                   significant == 0 ~ 0)) %>%
    mutate(sig_lesser = case_when(significant == 1 & exp_completeness < completeness ~ 1,
                                   significant == 1 & exp_completeness > completeness ~ 0,
                                   significant == 0 ~ 0)) -> family_coverage_tests


  length(which(family_coverage_tests$significant == 1))/nrow(family_coverage_tests) # 36.3% differ from expectation

  family_coverage_tests %>%
    group_by(family) %>%
    summarize( n_different = sum(significant),
               n_greater = sum(sig_greater),
               n_lesser = sum(sig_lesser),
               n_spp = unique(n_spp_in_fam)) %>%
    arrange(-n_different,n_greater) -> family_coverage_deviations

  family_coverage_tests %>%
    group_by(trait) %>%
    summarize( n_different = sum(significant),
               n_greater = sum(sig_greater),
               n_lesser = sum(sig_lesser)) %>%
    arrange(-n_different,n_greater) -> trait_coverage_deviations


#############################

  #Intraspecific


    # Trait completeness model: 1% threshold and broadly applicable traits
    #via a single model

    # Load overall trait completeness
    geo_traits_one_percent_threshold <- read_rds("data/focal_georeferenced_trait_coverage.rds")


    # merge with predictor variables
    socio_vars <- read.csv("manual_downloads/Darwinian_shortfalls/socioeco_var.csv")
    #socio_vars[2:10] <- scale(socio_vars[2:10])

    # combine trait completeness and predictor variables
    geo_traits_one_percent_threshold <-
      merge(x = geo_traits_one_percent_threshold,
            y = socio_vars,
            by.x= "area",
            by.y = "LEVEL_3_CO")

    # add species richness
    wcvp <- readRDS("manual_downloads/WCVP/wcvp_cleaned.RDS")

    wcvp %>%
      group_by(area_code_l3) %>%
      summarise(richness = n()) %>%
      merge(x = geo_traits_one_percent_threshold,
            y = .,
            by.x = "area",
            by.y = "area_code_l3",
            all.x = TRUE) ->
      geo_traits_one_percent_threshold

    # Add mean range size

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
      summarise(mean_species_range = mean(na.omit(range_size))) %>%
      merge(x = geo_traits_one_percent_threshold,
            y = .,
            by.x="area",
            by.y = "area_code_l3") -> geo_traits_one_percent_threshold

    # Add endemism

    wcvp %>%
      group_by(taxon_name) %>%
      summarize(n_countries = n()) %>%
      merge(x = wcvp,
            y = .) %>%
      group_by(area_code_l3) %>%
      filter(n_countries == 1) %>%
      summarise(endemics = n()) %>%
      merge(x = geo_traits_one_percent_threshold,
            y = .,
            by.x = "area",
            by.y = "area_code_l3",
            all.x = TRUE) %>%
      mutate(endemics = replace_na(endemics,0)) %>%
      mutate(endemism = endemics/richness) %>%
      select(-endemics) -> geo_traits_one_percent_threshold


    #rescale predictors
    geo_traits_one_percent_threshold[4:ncol(geo_traits_one_percent_threshold)] <- scale(geo_traits_one_percent_threshold[4:ncol(geo_traits_one_percent_threshold)])
    geo_traits_one_percent_threshold <- na.omit(geo_traits_one_percent_threshold) #shouldn't be any NAs, but just to be sure


    #Get coordinates needed for spatial
      geo_traits_one_percent_threshold %>%
        inner_join(y = tdwg_cents, by = c("area"= "LEVEL_3_CO")) -> geo_traits_w_coords

    #Get un-tf-ed richness
      wcvp %>%
        group_by(area_code_l3) %>%
        summarise(richness_untf = n()) %>%
        merge(x = geo_traits_w_coords,
              y = .,
              by.x = "area",
              by.y = "area_code_l3",
              all.x = TRUE) ->
        geo_traits_w_coords

    #Generate a list of island vs mainland polygons


      shape_sf <- st_read("manual_downloads/Darwinian_shortfalls/level3.shp")

      intersects <- st_intersects(shape_sf)
      intmat <- as.matrix(intersects)
      shape_sf$is_island <- 0
      shape_sf$is_island[which(rowSums(intmat)<=1)]<-1


        shape_sf %>%
          select(area = LEVEL_3_CO, is_island)%>%
          st_drop_geometry() %>%
          inner_join(geo_traits_w_coords) -> geo_traits_w_coords


      #trying spamm
      library(spaMM)

      geo_m_spamm <-
        geo_traits_w_coords %>%
        mutate(species_w_data = completeness * richness_untf) %>%
        mutate(species_wo_data = (1-completeness) * richness_untf) %>%
        fitme(cbind(species_w_data,species_wo_data) ~ AREA_SQKM + GDP_SUM + GDP_CAPITA + ROAD_DENSITY + POP_COUNT + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP
              + richness + mean_species_range + endemism +
                (1|trait) + Matern(1 | X + Y),
              data = .,
              family = "binomial") # for spamm, binomial models need to include successes and failures (e.g. samples and no samples)


      geo_m_spamm_null <-
        geo_traits_w_coords %>%
        mutate(species_w_data = completeness * richness_untf) %>%
        mutate(species_wo_data = (1-completeness) * richness_untf) %>%
        fitme(cbind(species_w_data,species_wo_data) ~  (1|trait) + Matern(1 | X + Y),
              data = .,
              family = "binomial") # for spamm, binomial models need to include successes and failures (e.g. samples and no samples)


      #Test whether model is significant improvement over null
      spaMM::LRT(object = geo_m_spamm,
                 object2 = geo_m_spamm_null) # p ~ 0

      geo_sims <- simulateResiduals(fittedModel = geo_m_spamm,integerResponse = TRUE,n = 1000)
      plot(geo_sims)

      #Look at model results
      summary(geo_m_spamm)

      #Check confidence intervals
        # confint(object = geo_m_spamm,
        #         parm = c("(Intercept)","AREA_SQKM","GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
        #                  "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
        #                  "EDUCATION_EXP",
        #                  "richness","mean_species_range","endemism"))
        #



      #Check confidence intervals using bs: convergence errors are thrown when using profiling

      #This is super slow, so I'm doing them individually
        # confint(object = geo_m_spamm,
        #         parm = c("(Intercept)","AREA_SQKM","GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
        #                  "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
        #                  "EDUCATION_EXP",
        #                  "richness","mean_species_range","endemism"),
        #         boot_args = list(nb_cores=2, nsim=199, seed=123))


      # confint(object = geo_m_spamm,
      #         parm = c("(Intercept)"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      #(-8.548, -8.055 )

      # confint(object = geo_m_spamm,
      #         parm = c("AREA_SQKM"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      #( 0.8464,  0.9900 )

      # confint(object = geo_m_spamm,
      #         parm = c("GDP_SUM"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # ( 0.3803,  0.5201 )

      # confint(object = geo_m_spamm,
      #         parm = c("GDP_CAPITA"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # ( 0.2308,  0.6863 )

      # confint(object = geo_m_spamm,
      #         parm = c("ROAD_DENSITY"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # (-0.3849, -0.0415 )

      # confint(object = geo_m_spamm,
      #         parm = c("POP_COUNT"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # ( 0.0161,  0.0723 )

      # confint(object = geo_m_spamm,
      #         parm = c("POP_DENSITY"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # (-0.6493, -0.1608 )

      # confint(object = geo_m_spamm,
      #         parm = c("SECURITY"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # ( 0.1544,  0.3886 )

      # confint(object = geo_m_spamm,
      #         parm = c("RESEARCH_EXP"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # ( 1.217,  1.623 )

      # confint(object = geo_m_spamm,
      #         parm = c("EDUCATION_EXP"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # ( 0.0675,  0.2363 )

      # confint(object = geo_m_spamm,
      #         parm = c("richness"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # ( 1.0157,  1.3320 )

      # confint(object = geo_m_spamm,
      #         parm = c("mean_species_range"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # (-0.3868,  0.3891 )

      # confint(object = geo_m_spamm,
      #         parm = c("endemism"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      # (-0.2434, -0.0449 )



      # ------------ Fixed effects (beta) ------------
      #   Estimate Cond. SE  t-value
      # (Intercept)        -8.18861   0.1795 -45.6066
      # AREA_SQKM           0.88467   0.1638   5.4002
      # GDP_SUM             0.39369   0.1681   2.3420
      # GDP_CAPITA          0.25658   0.2620   0.9794
      # ROAD_DENSITY       -0.26328   0.2883  -0.9132
      # POP_COUNT           0.02963   0.1736   0.1706
      # POP_DENSITY        -0.25743   0.2659  -0.9682
      # SECURITY            0.24383   0.2056   1.1862
      # RESEARCH_EXP        1.30454   0.2289   5.6984
      # EDUCATION_EXP       0.15835   0.1914   0.8273
      # richness            0.97937   0.2087   4.6927
      # mean_species_range -0.17278   0.2121  -0.8144
      # endemism           -0.07242   0.1919  -0.3773


##############################


    #Check richness
    general_traits_one_percent_threshold %>%
      select(area,richness)%>%
      unique() %>%
      inner_join(x = tdwg,y=.,by = c('LEVEL_3_CO'="area")) -> tdwg_w_richness


    ggplot()+
      geom_sf(data = tdwg_w_richness,mapping = aes(fill = richness))

    wcvp_w_unplaced %>%
      group_by(area_code_l3) %>%
      summarise(richness_w_unplaced = n()) %>%
      inner_join(x = tdwg_w_richness,y=.,by = c('LEVEL_3_CO'="area_code_l3")) -> tdwg_w_richness

    wcvp %>%
      group_by(area_code_l3) %>%
      summarise(richness_untf = n()) %>%
      inner_join(x = tdwg_w_richness,y=.,by = c('LEVEL_3_CO'="area_code_l3")) -> tdwg_w_richness
##############################################

############################

#venn diagram

  #need to know which species have any trait, dist, gene data


#Load data products from Rudbeck
  gen.data <- data.table::fread("manual_downloads/Darwinian_shortfalls/genbank_entries_w_duplicates_2022.csv", header = T, quote = "", sep = NULL) # As produced in the genbank_download.R script
  #gen.data <- fread("manual_downloads/Darwinian_shortfalls/genbank_entries_2022.csv", header = T, quote = "", sep = NULL) # As produced in the genbank_download.R script
  load("manual_downloads/Darwinian_shortfalls/BIEN_in_WCSP_regions_sep2021.RData") # R workspace with objects "spec.list" and "res". Spec.list is a list of all L3 regions with the species recorded there in BIEN, translated to WCSP ID's. Res is a df with the lengths of the elements of spec.list.
  name.id <- data.table::fread("data\\ID_and_Names_2022.csv", header = T) # As produced in the wcvp_subset_2021.R script
  rm(fin)


  # BIEN data from Rudbeck et al.

    bien.list <-
    lapply(X = names(spec.list),
           FUN = function(x){

             if(length(spec.list[[x]])==0){

               return(data.frame(country = x,
                          species = NA))

             }else{

               return(
                 data.frame(country = x,
                          species = spec.list[[x]])
                 )

             }

           }) %>%
      do.call(what = rbind) %>%
      inner_join(y = name.id,
                 by = c("species"="plant_name_id"))


  # Genetic data from Rudbeck
    gen.data <-
    gen.data %>%
      mutate(ID = 1:nrow(gen.data))

  # Pull names of species with trait data
    #tnrsed_trait_sp_names <- readRDS("data/tnrsed_names.RDS")
    trait_summary <- readRDS(file = "data/cleaning_raw_names/trait_summary_overall.RDS")

  # Pull the wcvp (without ferns for consistency with Rudbeck et al)

    wcvp_no_ferns <- readRDS("manual_downloads/WCVP/wcvp_cleaned_no_ferns.RDS") %>%
      dplyr::select(family, taxon_name) %>% unique()


  wcvp_no_ferns %>%
    mutate(trait_data = case_when(taxon_name %in% trait_summary$Accepted_species ~ 1,
                                  !taxon_name %in% trait_summary$Accepted_species ~ 0))%>%
    mutate(distribution_data = case_when(taxon_name %in% bien.list$taxon_name ~ 1,
                                  !taxon_name %in% bien.list$taxon_name ~ 0))%>%
    mutate(genetic_data = case_when(taxon_name %in% gen.data$species ~ 1,
                                  !taxon_name %in% gen.data$species ~ 0))-> data_availability


    trait_data_vec <- data_availability$taxon_name[which(data_availability$trait_data==1)]
    dist_data_vec <- data_availability$taxon_name[which(data_availability$distribution_data==1)]
    genetic_data_vec <- data_availability$taxon_name[which(data_availability$genetic_data==1)]

  #Euler
  library(eulerr)

data_availability%>%
  dplyr::select(c("trait_data","distribution_data","genetic_data"))%>%
  mutate(species = 1)%>%
  euler()->eu

eu %>%
  plot(labels = c("Trait","Distribution","Phylogeny","Seed Plants"),
       quantities =TRUE,
       fills=list(fill = c( "#FF80F7","#74ee15","#00D1D0",NA), alpha = 0.5))



cor.test(x = data_availability$trait_data,y = data_availability$distribution_data)
cor.test(x = data_availability$trait_data,y = data_availability$genetic_data)

library(Cairo)

Cairo(file="plots/euler_plot.jpg",
      type="jpg",
      units="in",
      width=10,
      height=10,
      dpi=300)


## Now render the plot
eu %>%
  plot(labels = c("Trait","Distribution","Phylogeny","Seed Plants"),
       quantities =TRUE,
       fills=list(fill = c( "#FF80F7","#74ee15","#00D1D0",NA), alpha = 0.5))


## When the device is off, file writing is completed.
dev.off()



############################


#Trait and variable correlations

library(corrplot)

# Load overall trait completeness
general_traits_one_percent_threshold <- read_rds("data/cleaning_raw_names/focal_trait_coverage.rds")

# merge with predictor variables
socio_vars <- read.csv("manual_downloads/Darwinian_shortfalls/socioeco_var.csv")
#socio_vars[2:10] <- scale(socio_vars[2:10])

# combine trait completeness and predictor variables
general_traits_one_percent_threshold <-
  merge(x = general_traits_one_percent_threshold,
        y = socio_vars,
        by.x= "area",
        by.y = "LEVEL_3_CO")

# add species richness
wcvp <- readRDS("manual_downloads/WCVP/wcvp_cleaned.RDS")

wcvp %>%
  group_by(area_code_l3) %>%
  summarise(richness = n()) %>%
  merge(x = general_traits_one_percent_threshold,
        y = .,
        by.x = "area",
        by.y = "area_code_l3",
        all.x = TRUE) ->
  general_traits_one_percent_threshold

# Add mean range size

socio_vars %>%
  dplyr::select(LEVEL_3_CO, AREA_SQKM) %>%
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
  summarise(mean_species_range = mean(na.omit(range_size))) %>%
  merge(x = general_traits_one_percent_threshold,
        y = .,
        by.x="area",
        by.y = "area_code_l3") -> general_traits_one_percent_threshold

# Add endemism

wcvp %>%
  group_by(taxon_name) %>%
  summarize(n_countries = n()) %>%
  merge(x = wcvp,
        y = .) %>%
  group_by(area_code_l3) %>%
  filter(n_countries == 1) %>%
  summarise(endemics = n()) %>%
  merge(x = general_traits_one_percent_threshold,
        y = .,
        by.x = "area",
        by.y = "area_code_l3",
        all.x = TRUE) %>%
  mutate(endemics = replace_na(endemics,0)) %>%
  mutate(endemism = endemics/richness) %>%
  dplyr::select(-endemics) -> general_traits_one_percent_threshold

general_traits_wide <-
general_traits_one_percent_threshold %>%
  pivot_wider(names_from = "trait", values_from = "completeness")

general_traits_wide <-
general_traits_wide %>%
  mutate("Mean Completeness" = rowMeans(dplyr::select(general_traits_wide,colnames(general_traits_wide)[14:66])))

correlations <- cor(general_traits_wide[14:67],method = "pearson")
correlations_long_names <- correlations
cnc <- colnames(correlations)
cnc_short <- sapply(X = cnc,FUN = function(x){substr(x = x,start = 1,stop = 22)})

  cnc_short <-
  lapply(X = 1:length(cnc),FUN = function(x){

    if(str_length(cnc_short)[x]!=str_length(cnc)[x]){
      cnc_short[x] <- paste(cnc_short[x],"...",sep = "")

    }else{
      cnc_short[x]
    }

  })%>%
    unlist()%>%
    as.character()


colnames(correlations) <- cnc_short
rownames(correlations) <- cnc_short
correlation_pvals <- cor.mtest(general_traits_wide[14:67],
                               conf.level = 0.95,
                               method="pearson",alternative="two.sided",exact = TRUE)



#correlations of individual traits with overall mean


corrplot( correlations ,
          addCoef.col="black",
          order = "FPC",
          tl.col = "black",
          p.mat = correlation_pvals$p,
          lowCI.mat = correlation_pvals$lowCI,
          uppCI.mat = correlation_pvals$uppCI,plotCI = "circle")

  #check that the names are correct
    names_df <- data.frame(colnames(correlations),colnames(correlation_pvals$p))
library(Cairo)

Cairo(file="plots/focal_trait_correlations.jpg",
      type="jpg",
      units="in",
      width=10,
      height=10,
      dpi=300)


## Now render the plot
corrplot( correlations ,
          order = "FPC",
          tl.col = "black",
          #tl.cex = 0.75,
          tl.cex = c(.9,rep(.75,53)),
          p.mat = correlation_pvals$p,
          pch.cex = .75)

## When the device is off, file writing is completed.
dev.off()


  #Get full names

corrplot( correlations_long_names ,
          order = "FPC",
          tl.col = "black",
          #tl.cex = 0.75,
          tl.cex = c(.9,rep(.75,53)),
          p.mat = correlation_pvals$p,
          pch.cex = .75)->corrplot_long_names


colnames(corrplot_long_names$corr)%>%
  writeLines(con = "plots/full_focal_correlation_names.txt",sep = ", ")


#correlations of predictor variables with overall mean

pred_corrs <- general_traits_wide[c(2:13,67)]
pred_corrs <- scale(pred_corrs)

pred_correlations <- cor(pred_corrs,method = "pearson")


pred_pvals <- cor.mtest(pred_corrs,
                               conf.level = 0.95,
                               method="pearson",alternative="two.sided",exact = TRUE)



#correlations of predictors (weak as expected)


corrplot( pred_correlations,
          addCoef.col="black",
          order = "FPC",
          tl.col = "black",
          p.mat = pred_pvals$p)




####################################



# Does trait coverage of all traits show the same pattern, or just the focal traits?

#Note that these coverage estimates are wrong since they compare all traits to the full set of species
#so, e.g. bark trait completeness doesn't accound for whether the species have bark

  total_trait_coverage <- readRDS("data/cleaning_raw_names/total_trait_coverage.RDS")

  total_coverage_wide <-
  total_trait_coverage %>%
    pivot_wider(names_from = trait,
                values_from = completeness)



  general_traits_one_percent_threshold <- read_rds("data/cleaning_raw_names/focal_trait_coverage.rds")

  general_traits_one_percent_threshold %>%
    group_by(area)%>%
    summarise(mean_completeness_focal = mean(completeness))%>%
    inner_join(total_coverage_wide)%>%
    dplyr::select(-area)%>%
    cor() -> correlations

  cor_focal <- correlations[,"mean_completeness_focal",drop=FALSE]%>%
    as.data.frame()%>%
    rownames_to_column(var = "TraitName")


  hist(cor_focal$mean_completeness_focal, breaks = 100)

  cor_focal %>%
    filter(mean_completeness_focal < 0) #370 traits with correlations less than zero


  trait_summary <- readRDS(file = "data/cleaning_raw_names/trait_summary_overall.RDS")


  species_per_trait <-
  trait_summary %>%
    group_by(TraitName)%>%
    summarise(n_species = n())

  cor_focal %>%
    inner_join(y = species_per_trait) -> cor_focal


  cor_focal %>%
    mutate(focal_trait = if_else(condition = TraitName %in% general_traits_one_percent_threshold$trait,
                                 true = "Yes",false = "No")) -> cor_focal


  cor_focal %>%
  ggplot(mapping = aes(x = n_species,
                       y = mean_completeness_focal,
                       color = focal_trait))+
    geom_point()+
    scale_x_log10(labels = scales::label_number())+
    geom_hline(yintercept = 0,lty=2)+
    xlab("Number of Species with Data")+
    ylab("Correlation with Mean Focal Trait Completeness")+
    labs(color = "Focal Trait")+
    theme_bw() -> cor_focal_plot

  ggsave(plot = cor_focal_plot,
         filename = "plots/all_traits_vs_focal_traits.jpg",width = 8,height = 5, bg = "white")



  # what fraction of traits show positive correlations with out mean trait coverage?

    length(which(cor_focal$mean_completeness_focal<0))/nrow(cor_focal) #negative
    length(which(cor_focal$mean_completeness_focal>0))/nrow(cor_focal) #positive


############################


    #Completeness vs climate change/uncertainty

    source("R/get_climate_data.R")

    climate_data <-
      get_climate_data(out_directory = "data/climate/",
                       out_parquet_directory = "data/climate_averages/",
                       tdwg = tdwg)

    gc()

    mean_completeness <-
      general_traits_one_percent_threshold %>%
      group_by(area)%>%
      summarise(mean_completeness = mean(completeness))

    unique(climate_data$changes_by_years$years)


    climate_data$changes_by_years%>%
      inner_join(y = mean_completeness,
                 by = c("country" = "area"))-> climate_by_year

    climate_data$changes_overall%>%
      inner_join(y = mean_completeness,
                 by = c("country" = "area"))-> climate_overall

    library(ggpmisc)
    library(ggpubr)

    bio_lookup <- data.frame(bio = 1:19,
                             variable =
                               c("Annual Mean Temperature",
                                 "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
                                 "Isothermality (BIO2/BIO7) (×100)",
                                 "Temperature Seasonality (standard deviation ×100)",
                                 "Max Temperature of Warmest Month",
                                 "Min Temperature of Coldest Month",
                                 "Temperature Annual Range (BIO5-BIO6)",
                                 "Mean Temperature of Wettest Quarter",
                                 "Mean Temperature of Driest Quarter",
                                 "Mean Temperature of Warmest Quarter",
                                 "Mean Temperature of Coldest Quarter",
                                 "Annual Precipitation",
                                 "Precipitation of Wettest Month",
                                 "Precipitation of Driest Month",
                                 "Precipitation Seasonality (Coefficient of Variation)",
                                 "Precipitation of Wettest Quarter",
                                 "Precipitation of Driest Quarter",
                                 "Precipitation of Warmest Quarter",
                                 "Precipitation of Coldest Quarter"),
                             variable_short =
                               c("Annual Mean Temp.",
                                 "Mean Diurnal Range)",
                                 "Isothermality",
                                 "Temperature Seasonality",
                                 "Max Temp. Warmest Month",
                                 "Min Temp. Coldest Month",
                                 "Temp. Annual Range",
                                 "Mean Temp. Wettest Quarter",
                                 "Mean Temp. Driest Quarter",
                                 "Mean Temp. Warmest Quarter",
                                 "Mean Temp. Coldest Quarter",
                                 "Annual Precip.",
                                 "Precip. Wettest Month",
                                 "Precip. Driest Month",
                                 "Precip. Seasonality",
                                 "Precip. Wettest Quarter",
                                 "Precip. Driest Quarter",
                                 "Precip. Warmest Quarter",
                                 "Precip. Coldest Quarter"))



    climate_overall %>%
      inner_join(y = bio_lookup,
                 by = c("bio"="bio"))%>%
      ggplot(mapping = aes(x = mean_difference, y = mean_completeness))+
      facet_wrap(~variable_short, scales = "free")+
      stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
      stat_poly_line() +
      geom_point()+
      ylab("Mean Trait Completeness")+
      xlab("Mean Difference from Present")

    climate_overall %>%
      inner_join(y = bio_lookup,
                 by = c("bio"="bio")) %>%
      ggplot(mapping = aes(x = sd_difference, y = mean_completeness))+
      facet_wrap(~variable_short, scales = "free")+
      stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
      stat_poly_line() +
      geom_point()+
      ylab("Mean Trait Completness")+
      xlab("SD Difference from Present")


    climate_overall %>%
      dplyr::select(bio,country,mean_completeness)->clim_comp
    climate_overall %>%
      dplyr::select(bio,country,mean_difference)->clim_mean
    climate_overall %>%
      dplyr::select(bio,country,sd_difference)->clim_sd

    clim_mean %>%
      inner_join(y = bio_lookup)%>%
      ungroup()%>%
      dplyr::select(-variable,-bio)%>%
      pivot_wider(names_from = variable_short,
                  values_from = mean_difference, names_prefix = "Diff. ")->clim_mean

    clim_sd %>%
      inner_join(y = bio_lookup)%>%
      ungroup()%>%
      dplyr::select(-variable,-bio)%>%
      pivot_wider(names_from = variable_short,
                  values_from = sd_difference, names_prefix = "SD. Diff. ")->clim_sd

    clim_comp %>%
      ungroup()%>%
      dplyr::select(country,mean_completeness)%>%unique()->clim_comp

    library(corrplot)


    clim_comp%>%
      inner_join(clim_mean)%>%
      inner_join(clim_sd)%>%
      dplyr::select(-country)%>%
      rename("Mean Completeness"= mean_completeness)%>%
      corrplot::cor.mtest()->clim_mean_p



    clim_comp%>%
      inner_join(clim_mean)%>%
      inner_join(clim_sd)%>%
      dplyr::select(-country)%>%
      rename("Mean Completeness" = mean_completeness)%>%
      na.omit()%>%
      cor()%>%
      corrplot( is.corr = TRUE,
                tl.col = "black",
                p.mat = clim_mean_p$p,
                pch.cex = .75,
                tl.cex = c(rep(.75,14),0.9,rep(.75,24)),
                order = "FPC")



    library(Cairo)

    Cairo(file="plots/correlations_w_climate.jpg",
          type="jpg",
          units="in",
          width=10,
          height=10,
          dpi=300)


    ## Now render the plot
    clim_comp%>%
      inner_join(clim_mean)%>%
      inner_join(clim_sd)%>%
      dplyr::select(-country)%>%
      rename("Mean Completeness" = mean_completeness)%>%
      na.omit()%>%
      cor()%>%
      corrplot( is.corr = TRUE,
                tl.col = "black",
                p.mat = clim_mean_p$p,
                pch.cex = .75,
                tl.cex = c(rep(.75,14),1,rep(.75,24)),
                order = "FPC")

    ## When the device is off, file writing is completed.
    dev.off()


    clim_comp%>%
      inner_join(clim_mean)%>%
      inner_join(clim_sd)%>%
      dplyr::select(-country)%>%
      rename("Mean Completeness" = mean_completeness)%>%
      na.omit()%>%
      cor()->test

    test <-
      data.frame(mean_completeness= test[,1],
                 p = round(clim_mean_p$p[,1],digits = 3))


    ###########################

    # Completeness vs human impact

    lbii <- raster::raster("manual_downloads/lbii.asc") #https://data.nhm.ac.uk/dataset/global-map-of-the-biodiversity-intactness-index-from-newbold-et-al-2016-science
    lbii
    crs(lbii) <- 4326
    plot(lbii) #higher values = more intact

    lbii_ea <-
      projectRaster(from = lbii,
                    res = 10000,
                    crs = 6933,
                    method = "bilinear")

    plot(lbii_ea)

    tdwg$biodiv_intactness <-
      raster::extract(x = lbii_ea,
                      y = tdwg) %>%
      lapply(X = ., FUN = function(x){mean(x, na.rm=TRUE)}) %>%
      unlist()



    human_footprint <- raster::raster("manual_downloads/hfp_global_geo_grid/hf_v2geo/") #https://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-footprint-geographic/data-download#openModal
    plot(human_footprint) #higher values = less intact

    human_footprint_ea <-
      projectRaster(from = human_footprint,
                    res = 10000,
                    crs = 6933,
                    method = "bilinear")

    tdwg$human_footprint <-
      raster::extract(x = human_footprint_ea,
                      y = tdwg) %>%
      lapply(X = ., FUN = function(x){mean(x, na.rm=TRUE)}) %>%
      unlist()


    general_traits_one_percent_threshold %>%
      group_by(area) %>%
      summarize(mean_trait_completeness = mean(completeness))%>%
      inner_join(y = tdwg,
                 by = c("area"="LEVEL_3_CO"))->tdwg_w_human

    plot(human_footprint_ea)
    cor.test(x = tdwg_w_human$mean_trait_completeness,
             y = tdwg_w_human$human_footprint)# Pearson = 0.18, p = 0.001

    plot(lbii_ea)
    cor.test(x = tdwg_w_human$mean_trait_completeness,
             y = tdwg_w_human$biodiv_intactness)# Pearson = - 0.04, p = 0.524








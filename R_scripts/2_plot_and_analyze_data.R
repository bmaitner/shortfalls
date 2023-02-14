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
traits$schema
traits$metadata





########################################################

# Trait completeness model: 1% threshold and broadly applicable traits
  #via a single model

  # Load overall trait completeness
    general_traits_one_percent_threshold <- read_rds("data/focal_trait_coverage.rds")


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


      # Moran I statistic standard deviate = 445.45, p-value < 2.2e-16
      # alternative hypothesis: two.sided
      # sample estimates:
      #   Moran I statistic       Expectation          Variance
      # 4.135883e-01     -4.940956e-05      8.622745e-07



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
      # statistic = 0.41359, observed rank = 1000, p-value < 2.2e-16
      # alternative hypothesis: two.sided



  moran_test_combined
  MC_test_combined
  MC_test_combined_two_sided

    plot(MC_test_combined)

    tdwg_dist <- st_distance(x = st_make_valid(tdwg_combined))
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


  m_spamm_w_area <-
    general_traits_w_coords %>%
    mutate(species_w_data = completeness * richness_untf) %>%
    mutate(species_wo_data = (1-completeness) * richness_untf) %>%
    fitme(cbind(species_w_data,species_wo_data) ~ AREA_SQKM + GDP_SUM + GDP_CAPITA + ROAD_DENSITY + POP_COUNT + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP
          + richness + mean_species_range + endemism +
            (1|trait) + Matern(1 | X + Y) + (1|area),
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

    summary(m_spamm_w_area)
  #Look at model results
    summary(m_spamm)


    # ------------ Fixed effects (beta) ------------
    #   Estimate Cond. SE  t-value
    # (Intercept)        -1.855277 0.154177 -12.0334
    # AREA_SQKM          -0.008144 0.014850  -0.5484
    # GDP_SUM            -0.006493 0.008456  -0.7678
    # GDP_CAPITA         -0.010972 0.008538  -1.2851
    # ROAD_DENSITY        0.007242 0.007795   0.9290
    # POP_COUNT           0.003991 0.010130   0.3939
    # POP_DENSITY         0.021014 0.010137   2.0730
    # SECURITY           -0.012117 0.010281  -1.1786
    # RESEARCH_EXP        0.060333 0.014348   4.2050
    # EDUCATION_EXP      -0.006242 0.007773  -0.8030
    # richness           -0.089006 0.014061  -6.3302
    # mean_species_range  0.399891 0.017347  23.0521
    # endemism           -0.128450 0.012497 -10.2784

    sims <- simulateResiduals(fittedModel = m_spamm,integerResponse = TRUE,n = 1000)
    ?simulateResiduals
    plot(sims)
    testOutliers(simulationOutput = sims,
                 type = "bootstrap")



  #Check confidence intervals
    confint(object = m_spamm,
                   parm = c("(Intercept)","AREA_SQKM","GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
                            "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
                            "EDUCATION_EXP",
                            "richness","mean_species_range","endemism"))

    #CIs

    # lower (Intercept) upper (Intercept)
    # -2.195835         -1.525245
    # lower AREA_SQKM upper AREA_SQKM
    # -0.03758786      0.02137629
    # lower GDP_SUM upper GDP_SUM
    # -0.02313497    0.01021238
    # lower GDP_CAPITA upper GDP_CAPITA
    # -0.027943928      0.005807283
    # lower ROAD_DENSITY upper ROAD_DENSITY
    # -0.008084352        0.022632342
    # lower POP_COUNT upper POP_COUNT
    # -0.01593169      0.02390306
    # lower POP_DENSITY upper POP_DENSITY
    # 0.00109606        0.04104828
    # lower SECURITY upper SECURITY
    # -0.032323488    0.008141573
    # lower RESEARCH_EXP upper RESEARCH_EXP
    # 0.03192910         0.08893986
    # lower EDUCATION_EXP upper EDUCATION_EXP
    # -0.021540199         0.009049627
    # lower richness upper richness
    # -0.11699809    -0.06098231
    # lower mean_species_range upper mean_species_range
    # 0.3639214                0.4361106
    # lower endemism upper endemism
    # -0.1540563     -0.1027810


  #Also check AICs
    AIC(m_spamm)
    AIC(m_spamm_null)


  plot(m_spamm) #documentation suggests this may be incorrect for binomials and recommends use of dharma package

  # Partial Dependence Plots
    plot_effects(m_spamm,"mean_species_range")
    plot_effects(m_spamm,"AREA_SQKM")


################################
  # Mean completeness: Trade-off: using mean completeness increases our R2, but decreases our power (only mean range size remains significant)
  #

    mean_completeness_glm <-
    general_traits_one_percent_threshold %>%
      group_by(area) %>%
      summarize(mean_trait_completeness = mean(completeness)) %>%
      merge(x = .,
            y = general_traits_one_percent_threshold,
            by="area")%>%
      select(-trait) %>%
      select(-completeness) %>%
      unique() %>%
      glm(formula = mean_trait_completeness ~ AREA_SQKM + GDP_SUM + GDP_CAPITA + ROAD_DENSITY + POP_COUNT + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP
            + richness + mean_species_range + endemism ,
            family = "binomial")

    summary(mean_completeness_glm)
    r2glmm::r2beta(model = mean_completeness_glm)

    performance::check_model()

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

        stop("Brian work on this bit. need to add ")

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

      # confint(object = geo_m_spamm,
      #         parm = c("AREA_SQKM"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("GDP_SUM"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("GDP_CAPITA"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("ROAD_DENSITY"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("POP_COUNT"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("POP_DENSITY"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))
      #
      # confint(object = geo_m_spamm,
      #         parm = c("SECURITY"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("RESEARCH_EXP"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("EDUCATION_EXP"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("richness"),
      #         boot_args = list(nb_cores=4, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("mean_species_range"),
      #         boot_args = list(nb_cores=4, nsim=100, seed=123))

      # confint(object = geo_m_spamm,
      #         parm = c("endemism"),
      #         boot_args = list(nb_cores=4, nsim=100, seed=123))



      # ------------ Fixed effects (beta) ------------
      #   Estimate Cond. SE  t-value
      # (Intercept)        -8.22194   0.1800 -45.6787
      # AREA_SQKM           0.90282   0.1657   5.4489
      # GDP_SUM             0.39607   0.1699   2.3313
      # GDP_CAPITA          0.22684   0.2655   0.8545
      # ROAD_DENSITY       -0.25824   0.2908  -0.8881
      # POP_COUNT           0.02202   0.1755   0.1255
      # POP_DENSITY        -0.23971   0.2692  -0.8906
      # SECURITY            0.26630   0.2084   1.2779
      # RESEARCH_EXP        1.31657   0.2315   5.6860
      # EDUCATION_EXP       0.15842   0.1939   0.8172
      # richness            0.96996   0.2109   4.5993
      # mean_species_range -0.22202   0.2152  -1.0317
      # endemism           -0.08735   0.1939  -0.4504

      # lower (Intercept) upper (Intercept)
      # -8.221943         -8.221941
      # lower AREA_SQKM upper AREA_SQKM
      # 0.902815        0.902817

  #################################

      #Geo data omitting islands

      geo_m_spamm_no_islands <-
        geo_traits_w_coords %>%
        filter(is_island == 0) %>%
        mutate(species_w_data = completeness * richness_untf) %>%
        mutate(species_wo_data = (1-completeness) * richness_untf) %>%
        fitme(cbind(species_w_data,species_wo_data) ~ AREA_SQKM + GDP_SUM + GDP_CAPITA + ROAD_DENSITY + POP_COUNT + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP
              + richness + mean_species_range + endemism +
                (1|trait) + Matern(1 | X + Y),
              data = .,
              family = "binomial") # for spamm, binomial models need to include successes and failures (e.g. samples and no samples)


      geo_m_spamm_null_no_islands <-
        geo_traits_w_coords %>%
        filter(is_island == 0) %>%
        mutate(species_w_data = completeness * richness_untf) %>%
        mutate(species_wo_data = (1-completeness) * richness_untf) %>%
        fitme(cbind(species_w_data,species_wo_data) ~  (1|trait) + Matern(1 | X + Y),
              data = .,
              family = "binomial") # for spamm, binomial models need to include successes and failures (e.g. samples and no samples)


      #Test whether model is significant improvement over null
      spaMM::LRT(object = geo_m_spamm_no_islands,
                 object2 = geo_m_spamm_null_no_islands) # p ~ 0

      # ------------ Fixed effects (beta) ------------
      #   Estimate Cond. SE   t-value
      # (Intercept)        -7.29017   0.2354 -30.96531
      # AREA_SQKM           0.79138   0.2531   3.12614
      # GDP_SUM             0.52173   0.2466   2.11553
      # GDP_CAPITA         -0.83286   0.4039  -2.06226
      # ROAD_DENSITY        0.03441   0.5469   0.06291
      # POP_COUNT          -0.08070   0.1693  -0.47664
      # POP_DENSITY        -0.22415   0.4507  -0.49731
      # SECURITY           -0.44147   0.2799  -1.57717
      # RESEARCH_EXP        1.93047   0.2938   6.57026
      # EDUCATION_EXP       0.09629   0.2096   0.45944
      # richness            0.53163   0.2360   2.25302
      # mean_species_range  0.02315   0.2339   0.09897
      # endemism            0.03499   0.3186   0.10981

      # confint(object = geo_m_spamm_no_islands,
      #         parm = c("(Intercept)"),
      #         boot_args = list(nb_cores=2, nsim=100, seed=123))

      # Intervals :
      #   Level      Normal              Basic              Percentile
      # 95%   (-7.498, -7.356 )   (-7.500, -7.352 )   (-7.228, -7.080 )

      # confint(object = geo_m_spamm_no_islands,
      #         parm = c("AREA_SQKM"),
      #         boot_args = list(nb_cores=4, nsim=100, seed=123))

      # Intervals :
      #   Level      Normal              Basic              Percentile
      # 95%   ( 0.8168,  0.9213 )   ( 0.8174,  0.9189 )   ( 0.6638,  0.7653 )

      # confint(object = geo_m_spamm_no_islands,
      #         parm = c("GDP_SUM"),
      #         boot_args = list(nb_cores=4, nsim=100, seed=123))
      #
      # Intervals :
      #   Level      Normal              Basic              Percentile
      # 95%   ( 0.4911,  0.5592 )   ( 0.4932,  0.5641 )   ( 0.4794,  0.5503 )
#
#       confint(object = geo_m_spamm_no_islands,
#               parm = c("GDP_CAPITA"),
#               boot_args = list(nb_cores=4, nsim=100, seed=123))

      # Intervals :
      #   Level      Normal              Basic              Percentile
      # 95%   (-0.9833, -0.7833 )   (-1.0008, -0.7961 )   (-0.8696, -0.6649 )

      # confint(object = geo_m_spamm_no_islands,
      #         parm = c("ROAD_DENSITY"),
      #         boot_args = list(nb_cores=4, nsim=100, seed=123))

      # Intervals :
      #   Level      Normal              Basic              Percentile
      # 95%   (-0.1258,  0.1551 )   (-0.1182,  0.1540 )   (-0.0851,  0.1870 )

      # confint(object = geo_m_spamm_no_islands,
      #         parm = c("POP_COUNT"),
      #         boot_args = list(nb_cores=4, nsim=100, seed=123))

      # Intervals :
      #   Level      Normal              Basic              Percentile
      # 95%   (-0.1216, -0.0637 )   (-0.1192, -0.0623 )   (-0.0991, -0.0422 )

      confint(object = geo_m_spamm_no_islands,
              parm = c("POP_DENSITY"),
              boot_args = list(nb_cores=4, nsim=100, seed=123))

      confint(object = geo_m_spamm_no_islands,
              parm = c("SECURITY"),
              boot_args = list(nb_cores=4, nsim=100, seed=123))

      confint(object = geo_m_spamm_no_islands,
              parm = c("RESEARCH_EXP"),
              boot_args = list(nb_cores=4, nsim=100, seed=123))

      confint(object = geo_m_spamm_no_islands,
              parm = c("EDUCATION_EXP"),
              boot_args = list(nb_cores=4, nsim=100, seed=123))

      confint(object = geo_m_spamm_no_islands,
              parm = c("richness"),
              boot_args = list(nb_cores=4, nsim=100, seed=123))

      confint(object = geo_m_spamm_no_islands,
              parm = c("mean_species_range"),
              boot_args = list(nb_cores=4, nsim=100, seed=123))

      confint(object = geo_m_spamm_no_islands,
              parm = c("endemism"),
              boot_args = list(nb_cores=4, nsim=100, seed=123))


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

    # should be ~ 119,405 wcvp species


  # Pull names of species with trait data
    tnrsed_trait_sp_names <- readRDS("data/tnrsed_names.RDS")


  # Pull the wcvp (without ferns for consistency with Rudbeck et al)

    wcvp_no_ferns <- readRDS("manual_downloads/WCVP/wcvp_cleaned_no_ferns.RDS") %>%
      dplyr::select(family, taxon_name) %>% unique()


  wcvp_no_ferns %>%
    mutate(trait_data = case_when(taxon_name %in% tnrsed_trait_sp_names$Accepted_species ~ 1,
                                  !taxon_name %in% tnrsed_trait_sp_names$Accepted_species ~ 0))%>%
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

############################


#Trait and variable correlations

library(corrplot)

# Load overall trait completeness
general_traits_one_percent_threshold <- read_rds("data/focal_trait_coverage.rds")

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


general_traits_wide %>%
  mutate(mean_completeness = rowMeans(dplyr::select(general_traits_wide,colnames(general_traits_wide)[14:67]))) -> general_traits_wide

correlations <- cor(general_traits_wide[14:68],method = "pearson")

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


colnames(correlations)<-cnc_short
rownames(correlations)<- cnc_short

correlation_pvals <- cor.mtest(general_traits_wide[14:68],
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


corrplot( correlations ,
          order = "FPC",
          tl.col = "black",
          tl.cex = 0.75,
          p.mat = correlation_pvals$p,
          pch.cex = .75)


#correlations of predictor variables with overall mean

pred_corrs <- general_traits_wide[c(2:13,68)]
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

  total_trait_coverage <- readRDS("data/total_trait_coverage.RDS")

  total_coverage_wide <-
  total_trait_coverage %>%
    pivot_wider(names_from = trait,
                values_from = completeness)



  general_traits_one_percent_threshold <- read_rds("data/focal_trait_coverage.rds")

  general_traits_one_percent_threshold %>%
    group_by(area)%>%
    summarise(mean_completeness_focal = mean(completeness))%>%
    inner_join(total_coverage_wide)%>%
    select(-area)%>%
    cor() -> correlations

  cor_focal <- correlations[,"mean_completeness_focal",drop=FALSE]%>%
    as.data.frame()%>%
    rownames_to_column(var = "TraitName")


  hist(cor_focal$mean_completeness_focal, breaks = 100)

  cor_focal %>%
    filter(mean_completeness_focal < 0) #370 traits with correlations less than zero


  trait_summary <- readRDS(file = "data/trait_summary_overall.RDS")


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
    xlab("Number of species with data")+
    ylab("Correlation with mean focal trait completeness")+
    labs(color = "Focal trait")


  # what fraction of traits show positive correlations with out mean trait coverage?

    length(which(cor_focal$mean_completeness_focal<0))/nrow(cor_focal)
    length(which(cor_focal$mean_completeness_focal>0))/nrow(cor_focal)

######################################

# get_citations <- function(directory, out_file = NULL){
#
#   packages <-
#   list.files(path = directory,
#              recursive = TRUE,
#              pattern = ".R$",
#              full.names = TRUE) %>%
#   questionr::qscan(load = FALSE) %>%
#     unlist() %>%
#     as.vector()
#
#
#   out <- sapply(packages,
#          FUN = function(x){
#
#           cite <- tryCatch(expr = toBibtex(citation(x)),
#                    error = function(e){})
#
#           if(!is.null(cite) & is.null(out_file)){
#             write(x = cite,
#                   file = out_file,
#                   append = TRUE)
#             }
#
#
#          })
#
#   return(out)
#
#
# }


############################









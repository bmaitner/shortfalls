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

      moran_test_two_sided <- moran.test(tdwg_combined$completeness,
                                        lw,
                                        alternative = "two.sided")


      MC_test_combined <- moran.mc(tdwg_combined$completeness,
                    lw,
                    nsim = 999,
                    alternative = "greater")

      MC_test_combined_two_sided <- moran.mc(tdwg_combined$completeness,
                                   lw,
                                   nsim = 999,
                                   alternative = "two.sided")

  moran_test_combined
  MC_test_combined
  MC_test_combined_two_sided

  ?moran.mc
  plot(MC_test_combined)
  plot(MC_test_combined_two_sided)

  #check whether the autocorrelation is present in the residuals

    tdwg_dist <- st_distance(x = st_make_valid(tdwg_combined))
    sims <- simulateResiduals(all_variables)
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
    # (Intercept)        -1.971818 0.155355 -12.6923
    # AREA_SQKM          -0.002878 0.015211  -0.1892
    # GDP_SUM            -0.006243 0.008596  -0.7263
    # GDP_CAPITA         -0.011658 0.008598  -1.3559
    # ROAD_DENSITY        0.007931 0.007888   1.0055
    # POP_COUNT           0.002478 0.010361   0.2392
    # POP_DENSITY         0.022487 0.010265   2.1907
    # SECURITY           -0.011157 0.010419  -1.0709
    # RESEARCH_EXP        0.059530 0.014590   4.0801
    # EDUCATION_EXP      -0.006837 0.007870  -0.8687
    # richness           -0.089564 0.014333  -6.2489
    # mean_species_range  0.385158 0.017740  21.7113
    # endemism           -0.125133 0.012753  -9.8120


  #Check confidence intervals
    confint(object = m_spamm,
                   parm = c("AREA_SQKM","GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
                            "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
                            "EDUCATION_EXP",
                            "richness","mean_species_range","endemism"))

    confint(object = m_spamm,
            parm = c("(Intercept)"))


    # lower (Intercept) upper (Intercept)
    # -2.313630         -1.635251
    # lower AREA_SQKM upper AREA_SQKM
    # -0.03304453      0.02733951
    # lower GDP_SUM upper GDP_SUM
    # -0.02317210    0.01076032
    # lower GDP_CAPITA upper GDP_CAPITA
    # -0.028761533      0.005237228
    # lower ROAD_DENSITY upper ROAD_DENSITY
    # -0.007579177        0.023494870
    # lower POP_COUNT upper POP_COUNT
    # -0.01790140      0.02284565
    # lower POP_DENSITY upper POP_DENSITY
    # 0.002317483       0.042809915
    # lower SECURITY upper SECURITY
    # -0.031635492    0.009369981
    # lower RESEARCH_EXP upper RESEARCH_EXP
    # 0.03068766         0.08857076
    # lower EDUCATION_EXP upper EDUCATION_EXP
    # -0.022325566         0.008648484
    # lower richness upper richness
    # -0.11823147    -0.06085545
    # lower mean_species_range upper mean_species_range
    # 0.3484147                0.4221169
    # lower endemism upper endemism
    # -0.15153308    -0.09865601

  #Also check AICs
    AIC(m_spamm)
    AIC(m_spamm_null)


  plot(m_spamm) #documentation suggests this may be incorrect for binomials and recommends use of dharma package

  #Partial Dependence Plots
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




#################################

    # Alternatively, may be able to fit each separatemodels for each trait using e.g. lme4::lmList()
      #this appears to be less powerful, since each model is fitted with little data.  I think its better to stick with random effects of trait

    library(lme4)

    lm_list_results <-
    lme4::lmList(formula = completeness ~ AREA_SQKM + GDP_SUM + GDP_CAPITA +
                                ROAD_DENSITY + POP_COUNT + POP_DENSITY +
                                SECURITY + RESEARCH_EXP + EDUCATION_EXP +
                                richness + mean_species_range + endemism |trait,
                 family = "binomial",
                 data = general_traits_one_percent_threshold)

    summary(lm_list_results)



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

      #Look at model results
      summary(geo_m_spamm)


      #Check confidence intervals
        # confint(object = geo_m_spamm,
        #         parm = c("(Intercept)","AREA_SQKM","GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
        #                  "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
        #                  "EDUCATION_EXP",
        #                  "richness","mean_species_range","endemism"))
        #

        confint(object = geo_m_spamm,
                parm = c("GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
                         "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
                         "EDUCATION_EXP",
                         "richness","mean_species_range","endemism"))





      #Check confidence intervals using bs: convergence errors are thrown when using profiling

      #This is super slow, so I'm doing them individually
        # confint(object = geo_m_spamm,
        #         parm = c("(Intercept)","AREA_SQKM","GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
        #                  "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
        #                  "EDUCATION_EXP",
        #                  "richness","mean_species_range","endemism"),
        #         boot_args = list(nb_cores=2, nsim=199, seed=123))


      confint(object = geo_m_spamm,
              parm = c("(Intercept)"),
              boot_args = list(nb_cores=2, nsim=100, seed=123))

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
          tl.cex = .75,
          order = "FPC")


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
         y = tdwg_w_human$biodiv_intactness)# Pearson = - 0.04, p = 0.497


############################

#venn diagram

  #need to know which species have any trait, dist, gene data


#Load data products from Rudbeck
  gen.data <- fread("manual_downloads/Darwinian_shortfalls/genbank_entries_w_duplicates_2022.csv", header = T, quote = "", sep = NULL) # As produced in the genbank_download.R script
  #gen.data <- fread("manual_downloads/Darwinian_shortfalls/genbank_entries_2022.csv", header = T, quote = "", sep = NULL) # As produced in the genbank_download.R script
  load("manual_downloads/Darwinian_shortfalls/BIEN_in_WCSP_regions_sep2021.RData") # R workspace with objects "spec.list" and "res". Spec.list is a list of all L3 regions with the species recorded there in BIEN, translated to WCSP ID's. Res is a df with the lengths of the elements of spec.list.
  name.id <- fread("data\\ID_and_Names_2022.csv", header = T) # As produced in the wcvp_subset_2021.R script
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




  library(VennDiagram)
    #https://r-graph-gallery.com/14-venn-diagramm.html

  vd<-
  venn.diagram(x = list(trait_data_vec[1:10],dist_data_vec[1:10],genetic_data_vec[1:10]),
               category.names = c("traits","distrbution","genes"),
               filename = "#venn1.svg",
               output=TRUE)


  set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
  set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
  set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

  set1 <- paste(trait_data_vec , sep="")
  set2 <- paste(dist_data_vec , sep="")
  set3 <- paste(genetic_data_vec , sep="")


  # Diagramt
  venn.diagram(
    x = list(set1, set2, set3),
    category.names = c("trait" , "dist" , "genetic"),
    filename = '#venn_diagramm.png',
    output=TRUE
  )

  venn.diagram(
    x = list(set1, set2, set3),
    category.names = c("trait" , "dist" , "genetic"),
    filename = '#venn_diagramm.png',
    output=TRUE
  )


  #Euler
  library(eulerr)

data_availability%>%
  dplyr::select(c("trait_data","distribution_data","genetic_data"))%>%
  mutate(species = 1)%>%
  euler()->eu

eu %>%
  plot(labels = c("Traits","Locations","Genes","Seed Plants"),
       quantities =TRUE,
       fills=list(fill = c("#FF80F7", "#00D1D0","#CFB000",NA), alpha = 0.5))



cor.test(x = data_availability$trait_data,y = data_availability$distribution_data)
cor.test(x = data_availability$trait_data,y = data_availability$genetic_data)

############################


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



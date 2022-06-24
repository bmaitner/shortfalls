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
          select(-endemics) -> general_traits_one_percent_threshold


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


#################################

  #Trivariate chloropleth maps

    general_traits_one_percent_threshold %>%
      group_by(area) %>%
      summarize(mean_trait_completeness = mean(completeness)) %>%
      merge(x = .,
            y = general_traits_one_percent_threshold,
            by="area")%>%
      select(-trait) %>%
      select(-completeness) %>%
      unique() -> mean_completeness




    load("manual_downloads/Darwinian_shortfalls/BIEN_in_WCSP_regions_sep2021.RData")


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
    geo_traits_one_percent_threshold <- na.omit(geo_traits_one_percent_threshold)


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
      confint(object = geo_m_spamm,
              parm = c("(Intercept)","AREA_SQKM","GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
                       "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
                       "EDUCATION_EXP",
                       "richness","mean_species_range","endemism"))


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



############################

  #Make plot of proportion of bryophytes

####################################

    DependentVar<-rnorm(60,mean = 9.5,sd=.5)
    Temperature<-rep(c(23.5,18,26.1,24.7,20.8,20),each=10)
    SpatData<-data.frame(x = rep(runif(6,0,100),each=10), y = rep(runif(6,0,100),each=10))
    SpatData$Temperature<-Temperature
    SpatData$DependentVar<-DependentVar
    library(stringr)

    SpatData$coords <- paste(SpatData$x,", ",SpatData$y)
    coords <- c(unique(SpatData$coords))
    x_unique <- c(str_extract(coords, "^.+(?=,)"))
    y_unique <- c(str_extract(coords, "(?<=, ).+$"))

    SpatLM <- lm(DependentVar~Temperature,data = SpatData)

    sims<-simulateResiduals(SpatLM)
    simsrecalc<-recalculateResiduals(sims,group = SpatData$coords)
    testSpatialAutocorrelation(simsrecalc, x = x_unique, y = y_unique)







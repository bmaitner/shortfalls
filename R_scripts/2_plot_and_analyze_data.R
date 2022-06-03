#' @author Brian Maitner
#' @description Scripts to wrangle data for the trait shortfalls manuscript

# Load libraries
library(tidyverse)
library(lme4)

# Get TRY data


library(sf)
#tdwg <- read_sf("manual_downloads/TDWG/level3.shp")
tdwg <- read_sf("manual_downloads/TDWG/old_lv3/level3.shp")
plot(tdwg[1])


#Load in the trait data
  traits <- arrow::open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

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

      MC_test_combined <- moran.mc(tdwg_combined$completeness,
                    lw,
                    nsim = 999,
                    alternative = "greater")

      MC_test_combined_v2 <- moran.mc(tdwg_combined$completeness,
                                   lw,
                                   nsim = 9999,
                                   alternative = "greater")



  moran_test_combined
  MC_test_combined

  ?moran.mc
  plot(MC_test_combined)

  #check whether the autocorrelation is present in the residuals

    tdwg_dist <- st_distance(x = st_make_valid(tdwg_combined))
    sims <- simulateResiduals(all_variables)
    testSpatialAutocorrelation(simulationOutput = sims,
                               distMat = tdwg_dist)


#################################

    #Trying glmmfields
    library(glmmfields) #doesn't seem to work with random effects

library(sf)
library(tidyverse)

  tdwg_cents <- st_centroid(st_make_valid(tdwg))
  tdwg_cents <- cbind(tdwg_cents,st_coordinates(tdwg_cents))

  general_traits_one_percent_threshold %>%
    inner_join(y = tdwg_cents,by = c("area"= "LEVEL_3_CO")) -> general_traits_w_coords



# m_spatial <- glmmfields(completeness ~ AREA_SQKM + GDP_SUM + GDP_CAPITA + ROAD_DENSITY + POP_COUNT + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP
#                         + richness + mean_species_range + endemism +
#                           (1|trait),
#                         data = general_traits_w_coords,
#                         family = binomial(),
#                         lat = "y", lon = "x",
#                         nknots = 12,
#                         iter = 500,
#                         chains = 1,
#                         prior_intercept = student_t(3, 0, 10),
#                         prior_beta = student_t(3, 0, 3),
#                         prior_sigma = half_t(3, 0, 3),
#                         prior_gp_theta = half_t(3, 0, 10),
#                         prior_gp_sigma = half_t(3, 0, 3),
#                         seed = 123 # passed to rstan::sampling()
# )

######################################################

#glmmtmb with coordinate distance

library(DHARMa)

all_variables <-
  glmer(data = general_traits_one_percent_threshold,
        formula = completeness ~ AREA_SQKM + GDP_SUM + GDP_CAPITA + ROAD_DENSITY + POP_COUNT + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP
        + richness + mean_species_range + endemism +
          (1|trait),
        family = "binomial",
        control = glmerControl(optimizer="bobyqa",
                               optCtrl=list(maxfun=2e5)),
        na.action = "na.fail")

  sims <- simulateResiduals(all_variables)
  plot(sims)
  sims2 <- recalculateResiduals(simulationOutput = sims,
                       group = general_traits_one_percent_threshold$trait)

  sa_test <-
  testSpatialAutocorrelation(simulationOutput = sims2,
                             x = unique(general_traits_w_coords$X),
                             y = unique(general_traits_w_coords$Y),
                             plot = FALSE)

  length(resid(all_variables))/length(unique(general_traits_one_percent_threshold$trait))


  library(lme4)

  testData = createData(sampleSize = 100, overdispersion = 0.5, family = poisson())
  fittedModel <- glmer(observedResponse ~ Environment1 + (1|group),
                       family = "poisson", data = testData)

  simulationOutput <- simulateResiduals(fittedModel = fittedModel)

  # standard plot
  plot(simulationOutput)

  # one of the possible test, for other options see ?testResiduals / vignette
  testDispersion(simulationOutput)

  # the calculated residuals can be accessed via
  residuals(simulationOutput)

  # transform residuals to other pdf, see ?residuals.DHARMa for details
  residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))

  # get residuals that are outside the simulation envelope
  outliers(simulationOutput)

  # calculating aggregated residuals per group
  simulationOutput2 = recalculateResiduals(simulationOutput, group = testData$group)
  plot(simulationOutput2, quantreg = FALSE)

  # calculating residuals only for subset of the data
  simulationOutput3 = recalculateResiduals(simulationOutput, sel = testData$group == 1 )
  plot(simulationOutput3, quantreg = FALSE)


################################

  wcvp %>%
    group_by(area_code_l3) %>%
    summarise(richness_untf = n()) %>%
    merge(x = general_traits_w_coords,
          y = .,
          by.x = "area",
          by.y = "area_code_l3",
          all.x = TRUE) ->
    general_traits_w_coords



  all_variables <-
    glmer(data = general_traits_one_percent_threshold,
          formula = completeness ~ AREA_SQKM + GDP_SUM + GDP_CAPITA + ROAD_DENSITY + POP_COUNT + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP
          + richness + mean_species_range + endemism +
            (1|trait),
          family = "binomial",
          control = glmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=2e5)),
          na.action = "na.fail")


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


  #Check confidence intervals
    confint(object = m_spamm,
                   parm = c("AREA_SQKM","GDP_SUM","GDP_CAPITA","ROAD_DENSITY",
                            "POP_COUNT","POP_DENSITY","SECURITY","RESEARCH_EXP",
                            "EDUCATION_EXP",
                            "richness","mean_species_range","endemism"))

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







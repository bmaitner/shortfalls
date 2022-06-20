library(parzer)
library(sf)



get_countries <- function(useful_md, tdwg){

  #Come up with a single, consensus lat/lon

    #First, take decimal lat/long where possible

    useful_md %>%
      mutate(assigned_lat = case_when(!is.na(`Latitude (decimal degrees)`) ~ `Latitude (decimal degrees)`),
             assigned_lon = case_when(!is.na(`Longitude (decimal degrees)`) ~ `Longitude (decimal degrees)`))%>%
    filter(!is.na(assigned_lon) & !is.na(assigned_lat))%>%
    mutate(assigned_lat = as.numeric(assigned_lat),
           assigned_lon = as.numeric(assigned_lon))-> useful_md_working

    #check that all longitude grid observations are a subset of longitude
      if(any(which(!is.na(useful_md$`Longitude grid`) & is.na(useful_md$Longitude)))){stop("write more code")}

    #assign any convenient lat/lon values (ignore the dms format for now)
      next_batch <-
      useful_md%>%
        filter(!ObservationID %in% useful_md_working$ObservationID)

      next_batch %>%
        mutate(assigned_lat = case_when(!is.na(as.numeric(Latitude)) ~ as.numeric(Latitude))) %>%
        mutate(assigned_lon = case_when(!is.na(as.numeric(Longitude)) ~ as.numeric(Longitude))) %>%
        filter(!is.na(assigned_lon) & !is.na(assigned_lat))-> next_batch

      useful_md_working %>%
        bind_rows(next_batch) -> useful_md_working

    #handle the harder stuff using a parser

      next_batch <-
        useful_md%>%
        filter(!ObservationID %in% useful_md_working$ObservationID)

      next_batch$parse_lat <-
      parzer::parse_lat(lat = next_batch$Latitude)

      next_batch$parse_lon <-
        parzer::parse_lon(lon = next_batch$Longitude)

      next_batch %>%
        mutate(assigned_lat = case_when(!is.na(parse_lat) ~ as.numeric(parse_lat))) %>%
        mutate(assigned_lon = case_when(!is.na(parse_lon) ~ as.numeric(parse_lon))) %>%
        filter(!is.na(assigned_lon) & !is.na(assigned_lat))-> next_batch


      useful_md_working %>%
        bind_rows(next_batch %>%select(colnames(useful_md_working))) -> useful_md_working

      #give up on anything without lat/long

      next_batch <-
        useful_md%>%
        filter(!ObservationID %in% useful_md_working$ObservationID)

      next_batch %>%
        filter(is.na(Latitude) | is.na(Longitude))%>%
        mutate(assigned_lat = NA,
               assigned_lon = NA) -> next_batch


      useful_md_working %>%
        bind_rows(next_batch %>%select(colnames(useful_md_working))) -> useful_md_working


      #whats left: anomalously formatted lat/lons

      next_batch <-
        useful_md %>%
        filter(!ObservationID %in% useful_md_working$ObservationID)


      message(round(nrow(next_batch)/nrow(useful_md)*100,digits = 2), " percent of MD being tossed due to unparsable lat/long")
      message("==", round(nrow(next_batch)/length(which(!is.na(useful_md$Latitude) & !is.na(useful_md$Longitude)))*100,digits = 2), " percent of MD with coordinates being tossed due to unparsable lat/long")



    # Now, assign botanical country to places with coordinates

      useful_md_working %>%
        select(assigned_lat, assigned_lon) %>%
        filter(!is.na(assigned_lat) & !is.na(assigned_lon)) %>%
        unique()%>%
        st_as_sf(coords = c("assigned_lon","assigned_lat"),
                 crs = st_crs("WGS84"),
                 remove = FALSE) -> point_data


      tdwg_countries <- st_intersection(x = point_data,
                                        y = st_make_valid(tdwg))

      tdwg_countries %>%
        st_set_geometry(NULL) -> tdwg_countries

        useful_md_working %>%
          left_join(y = tdwg_countries) -> useful_md_working

        rm(tdwg_countries, point_data)
        gc()


    # use GNRS to clean both the the tdwg and TRY country/state names and try to match them up

        library(GNRS)

        useful_md_working %>%
          select(`Location Country`,`Location: state`) %>%
          unique() -> try_countries

      GNRSed_TRY_countries <-
        GNRS_super_simple(country = try_countries$`Location Country`,
                          state_province = try_countries$`Location: state`)

      GNRSed_TRY_countries %>%
        filter(!is.na(country)) %>%
        filter(country != "") %>%
        select(country_verbatim, state_province_verbatim, country,state_province) -> GNRSed_TRY_countries


      #First, we'll try to match TDWG level names
        gnrsed_bot_countries <- GNRS_super_simple(country = tdwg$LEVEL_NAME)

        gnrsed_bot_countries %>%
          filter(country != "") -> gnrsed_bot_countries


        gnrsed_bot_countries <-
        gnrsed_bot_countries %>%
          filter(match_score_country > .95 | #keep anything with a high country match score
                  grepl(pattern = "Island",x = country)) %>% #keep islands (the abbreviation hurts the score)
          filter(match_score_state_province == "") #this filter excludes NSW from showing up as Wales

      #Next, try to get any that correspond to states
        gnrsed_bot_countries_v2 <- GNRS_super_simple(country = tdwg$REGION_NAM,state_province = tdwg$LEVEL_NAME)

        gnrsed_bot_countries_v2 <-
        gnrsed_bot_countries_v2 %>%
          filter(!user_id %in% gnrsed_bot_countries$user_id) %>% #toss anything already matched
          filter(country != "" & state_province != "") %>%  #toss anything that wasn't resolved
          filter(!country_verbatim %in% c("Southern Africa")) %>% #erroneously resolved to South Africa
          filter(!country %in% c("Mexico")) #toss mexico (erroneously resolves states)


        gnrsed_bot_countries %>%
          bind_rows(gnrsed_bot_countries_v2) -> gnrsed_bot_countries_total

        rm(gnrsed_bot_countries,gnrsed_bot_countries_v2)

      #Now merge the TRY countries to the TDWG countries to create a name key
        GNRSed_TRY_countries %>%
          rename(TRY_country_name = country_verbatim,
                 TRY_state_name = state_province_verbatim) %>%
          merge(y = gnrsed_bot_countries_total %>% rename(TDWG_country_name = country_verbatim, TDWG_state_name = state_province_verbatim  ),
                by = c("country","state_province") ) -> name_key

        rm(GNRSed_TRY_countries, gnrsed_bot_countries_total)

      #Merge the useful fields (TDWG_country_name) of the name key to the useful_md_working

          name_key %>%
            select(TRY_country_name,
                   TRY_state_name,
                   TDWG_country_name,
                   TDWG_state_name) %>%
            right_join(y = useful_md %>% mutate(`Location: state` = replace_na(`Location: state`,"")),
                       by = c("TRY_country_name" = "Location Country",
                            "TRY_state_name"  = "Location: state"),
                       keep=TRUE) -> useful_md_plus

      # Merge the useful fields of the extracted bits to the full dataset
          useful_md_working %>%
            rename("TDWG_country_imputed" = "LEVEL_NAME") %>%
            select("ObservationID", "TDWG_country_imputed")%>%
            right_join(y = useful_md_plus,
                       by = "ObservationID") -> useful_md_working
          rm(useful_md_plus)
      # Where we can infer country from multiple sources, how many are accurate?

        message("\n ",round(length(which(useful_md_working$TDWG_country_imputed != useful_md_working$TDWG_country_name))/
          length(which(!is.na(useful_md_working$TDWG_country_imputed) & !is.na(useful_md_working$TDWG_country_name))),digits = 2)*100,
          " percent of countries inferrable from both points and declared country do not match, tossing these")

        useful_md_working %>%
          filter(is.na(TDWG_country_imputed) | is.na(TDWG_country_name)|
                   TDWG_country_imputed == TDWG_country_name) -> useful_md_working

        #Toss anything where we couldn't infer a country

        useful_md_working %>%
          filter(!is.na(TDWG_country_imputed) | !is.na(TDWG_country_name)) -> useful_md_working


        #Now, assign a single country to the md
          useful_md_working %>%
            select(ObservationID, TDWG_country_imputed, TDWG_country_name) -> useful_md_working

          useful_md_working %>%
            mutate(LEVEL_NAME = case_when(!is.na(TDWG_country_imputed) ~ TDWG_country_imputed,
                                          !is.na(TDWG_country_name) & is.na(TDWG_country_imputed) ~ TDWG_country_name)) %>%
            select(ObservationID, LEVEL_NAME) -> useful_md_working


        # Append the assigned country to the md and return
          useful_md %>%
            left_join(y = useful_md_working,
                      by = "ObservationID")-> useful_md

        return(useful_md)


}

























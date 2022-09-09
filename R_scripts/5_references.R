
traits <- arrow::open_dataset(sources = "manual_downloads/TRY/TRY_parquet/")

traits %>%
  dplyr::select(Reference)%>%
  collect()%>%
  unique() -> try_refs

austraits <- austraits::load_austraits(path = "data/austraits/",
                                       version = "3.0.2")

aus_refs <- austraits$sources


citation("parzer")

citation("GNRS")

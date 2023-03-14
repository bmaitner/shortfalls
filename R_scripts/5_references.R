
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

#Make sure we get all R functions

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

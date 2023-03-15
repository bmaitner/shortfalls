library(tidyverse)


#Make sure we get all R packages

  get_citations <- function(directory, out_file = NULL){

    packages <-
    list.files(path = directory,
               recursive = TRUE,
               pattern = ".R$",
               full.names = TRUE) %>%
    questionr::qscan(load = FALSE) %>%
      unlist() %>%
      as.vector()%>%
      unique()


    out <- sapply(packages,
           FUN = function(x){

            cite <- tryCatch(expr = toBibtex(citation(x)),
                     error = function(e){})

            if(!is.null(cite) & !is.null(out_file)){
              write(x = cite,
                    file = out_file,
                    append = TRUE)
              }


           })


    return(out)


  }


get_citations(directory = getwd(),
              out_file = "data/r_packages_used.bib")


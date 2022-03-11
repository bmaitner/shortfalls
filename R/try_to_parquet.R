
library(arrow)
library(data.table)


#' @author Brian Maitner
#' @description Quick and dirty function to convert TRY's .txt files to a more useful parquet format
#' @param file A file to be converted
#' @param output_directory Where you want all the files deposited
#' @param batch_size The number of lines to load in at once
try_to_parquet <- function(file,
                           output_directory = "manual_downloads/TRY/TRY_parquet/",
                           batch_size = 80000){

  #Setup variables

    i <- 0
    error_found <- FALSE

  #Iterate through batches

    while(!error_found){

        data <- fread(file = file,
            nrows = batch_size,
            skip = i)

        if(i == 0){

          col_names <- colnames(data)
        }else{

          colnames(data) <- col_names

        }



      tryCatch(
        expr = write_parquet(x = data,
                             sink = file.path(output_directory,paste(basename(file),".",as.integer(i),".gz.parquet",sep = "")),
                             compression = "gzip"),
        error = function(e){
          error_found <- TRUE
          message(paste("Finished converting TRY file",file,"to parquet"))
          return(invisible(NULL))

        }


      )


      #check whether you're done
      if(nrow(data) < batch_size){

        message(paste("Finished converting TRY file",file,"to parquet"))
        return(invisible(NULL))

      }


      i <- i+nrow(data)

      print(i)

      rm(data)


    } #while loop




}# end fx



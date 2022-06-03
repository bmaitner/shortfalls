library(raster)
library(rasterVis)
library(animation)
library(latticeExtra)
library(stars)
library(ggplot2)
library(gganimate)
ani.options(convert = "C:/Program Files/ImageMagick-7.1.0-Q16-HDRI/convert.exe")

shortfall_animations <- function(countries, format=c("gif","mp4")){


  saveGIF({
    for(i in 5:55){

      #countries[,5]

      plot_i <-
      ggplot(data = countries[,i])+
        geom_sf(aes(fill =  100*  get(colnames(countries)[i])))+
        scale_fill_viridis_c(option = "plasma")+
        labs(fill = "%")+
        ggtitle(  colnames(countries)[i]  )

      plot(plot_i)
      ani.pause()
      #save plot

      #paste("Unicron/biome_",i,sep = "")


    }#i biomes loop




  },

  if("mp4" %in% format){
    movie.name = paste("shortfalls",".mp4",sep = "")
  }else{movie.name = paste("shortfalls",".gif",sep = "") }


  )


}

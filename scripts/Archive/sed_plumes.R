gc()

library(tidyverse)
library(dplyr)
library(rgdal)
library(raster)
library(gdistance)
library(reshape2)
library(sf)
library(seegSDM)

source('C:/Data/50Reefs/Felipe/Codes/function_plume_diffusion_model_2.R')

setwd("C:/Data/50Reefs/Felipe")

#READ OCEANS RASTER
ocean_raster<-raster("./Sed_plumes/ocean_raster_1km.tif")

#read other plumes
plu<-raster("./rivers_coast/POUR_POINTS/CAITIE/ero_dep_tiff.tif")

#READ POUR POINTS SHAPEFILE
pour_points<-st_read("./rivers_coast/POUR_POINTS/CAITIE/All_nutrients_at_pourpoints.shp")   #50R_pp 

pour_points<-spTransform(as(pour_points,"Spatial"), crs(plu))
ocean_raster <- projectRaster(ocean_raster, crs = plu)

#pour_points<-pour_points[c(1,2)]

#two options to apply the function
re<-list()#to save your rasters
pb <- txtProgressBar(min = 0, max = length(pour_points), style = 3)
for (i in 1:nrow(pour_points)){
  setTxtProgressBar(pb, i) 
  re[[i]]<-  sed_flow(  
    pp_sf=pour_points[i,],#list of pour points (each element of the list should be a point)
    initial.value = pour_points[i,"N"][[1]],#value with sediment export value
    ocean_raster,
    buffer_value=60000,#buffer from pour point (in meters)-change value to make nutrients go farther 
    resolution_raster=1,#in km
    searchdistance=9000,#search distance to snap points to the coast,
    decay_percent=0.005,#5 % of sediment decay
    remaining_left=0.0005#stop when 99.5% of the sediment has spread from the coast
  )
}
gc()
#add all rasters and save

newextent<-extent(Reduce(extend,re))

re = lapply(re, function(r){extend(r,newextent)})

#set all rasters to same extent
for (i in 1:length(re)){
  r<-re[[i]]
  r<-extend(r,newextent)
  r[is.na(r[])] <- 0 
  re[[i]]<-r
  print(i)
}

save(re,file = "./rivers_coast/POUR_POINTS/CAITIE/N_plumes_fiji.R")

#merge all rasters
sed_plume_merge<-Reduce("+", re)

gc()

plot(sed_plume_merge)

#save rasters

writeRaster(sed_plume_merge,"./rivers_coast/POUR_POINTS/CAITIE/N_plumes_fiji.tif",overwrite = T)

gc()

library(tidyverse)
library(dplyr)
library(rgdal)
library(raster)
library(gdistance)
library(reshape2)
library(sf)
#library(seegSDM)#discontinued

source('function_plume_diffusion_model_2.R')

#READ OCEANS RASTER
ocean_raster<-raster("ocean_raster_1km.tif")

#read other plumes
plu<-raster("ero_dep_tiff.tif")

#READ POUR POINTS SHAPEFILE
pour_points<-st_read("All_nutrients_at_pourpoints.shp")   #50R_pp 

pour_points<-spTransform(as(pour_points,"Spatial"), crs(plu))
ocean_raster <- projectRaster(ocean_raster, crs = plu)

pour_points<-st_as_sf(pour_points)

#pour_points<-pour_points[c(1,2)]

#This code works well for a relatively small area.
#If you want to work at a global scale, you need to divide them
#by groups using something like this:

#divide pour points by groups depending on location
p_extent=st_as_sfc(st_bbox(pour_points))
p_extent_sub = st_make_grid(p_extent,n = c(40,20)) %>%
  st_as_sf() %>% 
  #rename(geom=x) %>%
  mutate(group=1:nrow(.))
#tm_shape(p_extent_sub)+tm_polygons()
pour_points = st_join(st_as_sf(pour_points), p_extent_sub, join = st_intersects) %>% 
  mutate(ID=1:nrow(.)) 

#convert to spatial object
pour_points<-as(pour_points,"Spatial")
#create list by group
#point_list<-split(pour_points,as.character(pour_points@data$group))

point_list<-split(pour_points,as.character(pour_points$group))


#The you can apply the function to each group in the
#list

#I used a loop because I wanted to keep track of some
#pour points. It is fairly quickly but you can also paralellize.

pb <- txtProgressBar(min = 0, max = length(point_list), style = 3)
sed_plume_merge<-list()
for (p in 100:length(point_list)) {
  setTxtProgressBar(pb, p)
  focus_points<-point_list[[p]]
tb <- txtProgressBar(min = 0, max = nrow(focus_points), style = 3)
re<-list()#to save your rasters
for (i in 1:nrow(focus_points)){#only the first 10 pour points
  setTxtProgressBar(tb, i) 
  re[[i]]<-  sed_flow(  
    pp_sf=focus_points[i,],#list of pour points (each element of the list should be a point)
    initial.value = focus_points[i,"N"][[1]],#value with sediment export value
    ocean_raster,
    buffer_value=60000,#buffer from pour point (in meters)-change value to make nutrients go farther 
    resolution_raster=1,#in km
    searchdistance=9000,#search distance to snap points to the coast,
    decay_percent=0.005,#5 % of sediment decay
    remaining_left=0.0005#stop when 99.5% of the sediment has spread from the coast
  )
  gc()
}

newextent<-extent(Reduce(extend,re))

re = lapply(re, function(r){extend(r,newextent)})

#set all rasters to same extent
for (i in 1:length(re)){
  r<-re[[i]]
  r<-extend(r,newextent)
  r[is.na(r[])] <- 0 
  re[[i]]<-r
  #print(i)
}
#merge all rasters
sed_plume_merge[[p]]<-Reduce("+", re)
gc()
print(p)
}
#add all rasters and save

sed<-plyr::compact(sed_plume_merge)

newextent<-extent(Reduce(extend,sed))

re = lapply(sed, function(r){extend(r,newextent)})

for (i in 1:length(re)){
  r<-re[[i]]
  r<-extend(r,newextent)
  r[is.na(r[])] <- 0 
  re[[i]]<-r
  print(i)
}


save(re,file = "N_plumes_fij_1.R")

#merge all rasters
re<-Reduce("+", re)

gc()

plot(re)

#save rasters

writeRaster(re,"N_plumes_fiji_1.tif",overwrite = T)

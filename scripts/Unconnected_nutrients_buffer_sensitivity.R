library(tidyverse)
library(sf)
#library(startR)
library(here)
library(raster)
library(gdalUtils)
library(tabularaster)
library(units)
library(rgdal)
library(lwgeom)
library(pbmcapply)



# STP watersheds

stp<-st_read(here("raw_data/Sewage_plants/Fiji Waste Water Catchment.shp")) %>% 
  mutate(STP = ifelse(TREATMENT_ == "NATABUA WWTP", "Natabua WWTP",
                      ifelse(TREATMENT_ == "Navakai PS", "Navakai WWTP", as.character(TREATMENT_)))) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(STP) %>% 
  summarise(area = sum(area, na.rm = T)) %>% 
  filter(!is.na(STP) == T) 


# River data

rivers_proj<-st_read(here("projected_data/Rivers/Fiji_rivers.shp"))
rivers_proj<-st_transform(rivers_proj, st_crs(stp))

# N and P pollution

N_cell<-raster(here("output_data/Nutrients/Fiji_N_per_cell.tif"))
P_cell<-raster(here("output_data/Nutrients/Fiji_P_per_cell.tif"))

# Buffer river

#For now I do 500 m each side - determine a more defensible value?
  
#  Could also buffer coastline?
  
#  Could also do something like Steph Borelle's paper on plastics:

#https://science.sciencemag.org/content/suppl/2020/09/16/369.6510.1515.DC1


buff<-seq(from = 100, to = 2000, by = 100) %>% as.list()

create_buff<-function(x){
  
buffer_rivers<-st_buffer(rivers_proj, dist = x) %>% 
  mutate(ID = 1)

diff_buffer<-st_difference(buffer_rivers, st_union(stp))

st_write(diff_buffer, paste0("/home/kuempel/Fiji_sewage/output_data/Rivers/buffers/River_buffer_stprm_",x,".shp"), append = F)

buffer_ras<-rasterize(diff_buffer, N_cell, field = "ID")

writeRaster(buffer_ras, paste0("/home/kuempel/Fiji_sewage/output_data/Rivers/buffers/buffer_rivers_",x,".tif"),overwrite = T)

}

pbmclapply(buff, create_buff, mc.cores = 5, mc.style = "ETA")


rm(list=ls())
gc()

# N and P pollution

N_cell<-raster(here("output_data/Nutrients/Fiji_N_per_cell.tif"))
P_cell<-raster(here("output_data/Nutrients/Fiji_P_per_cell.tif"))

# Nutrients within watersheds

#Watershed boundaries

#Level 12 basins

# Read in data and calculate area

fiji_basins12<-st_read(here("projected_data/Basins/Fiji_basins12.shp")) %>% 
  group_by(MAIN_BAS) %>% 
  summarise(area_km2 = sum(area_km2, na.rm = T))

fiji_basins12<-fiji_basins12 %>% 
  mutate(ID = seq(from = 1, to =nrow(fiji_basins12), 1))


# Change Buffer area to N and P pollution values


buffer_list<-list.files(path = here("output_data//Rivers/buffers"), pattern = "buffer_rivers_", full.names = T)

buffer_N_P<-function(x){
  
  buffer_ras<-raster(x)
  
  N_unconnected<-N_cell * buffer_ras
  P_unconnected<-P_cell * buffer_ras
  
  basin_N_unconnected<-raster::extract(N_unconnected, fiji_basins12, fun=sum, na.rm=TRUE, df=TRUE)
  basin_P_unconnected<-raster::extract(P_unconnected, fiji_basins12, fun=sum, na.rm=TRUE, df=TRUE)
  
  basins_nutri_df<-full_join(basin_N_unconnected, fiji_basins12, by = "ID") %>% 
    rename(N = layer) %>% 
    full_join(., basin_P_unconnected, by = "ID") %>% 
    rename(P = layer) %>% 
    dplyr::select(ID, N, P, MAIN_BAS)
  
  write.csv(basins_nutri_df, paste0("/home/kuempel/Fiji_sewage/output_data/Nutrients/buffers/N_P_unconnected_per_basin_",x,"_buff.csv"))
}

pbmclapply(buffer_list, buffer_N_P, mc.cores = 15, mc.style = "ETA")






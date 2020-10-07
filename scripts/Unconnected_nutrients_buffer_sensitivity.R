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
  filter(!is.na(STP) == T,
         !STP %in% c("ACS WWTP", "Wailada WWTP", "Naboro WWTP")) 


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

N_cell_septic<-raster(here("output_data/Nutrients/Septic_N.tif"))
P_cell_septic<-raster(here("output_data/Nutrients/Septic_P.tif"))
N_cell_direct<-raster(here("output_data/Nutrients/Direct_N.tif"))
P_cell_direct<-raster(here("output_data/Nutrients/Direct_P.tif"))

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

buffer_N_P<-function(x,
                     N_ras,
                     P_ras,
                     name){
  
  buffer_ras<-raster(x)
  num<-as.numeric(regmatches(x, gregexpr("[[:digit:]]+", x)))
  
  N<-N_ras * buffer_ras
  P<-P_ras * buffer_ras
  
  basin_N<-raster::extract(N, fiji_basins12, fun=sum, na.rm=TRUE, df=TRUE)
  basin_P<-raster::extract(P, fiji_basins12, fun=sum, na.rm=TRUE, df=TRUE)
  
  basins_nutri_df<-full_join(basin_N, fiji_basins12, by = "ID") %>% 
    rename(N = layer) %>% 
    full_join(., basin_P, by = "ID") %>% 
    rename(P = layer) %>% 
    dplyr::select(ID, N, P, MAIN_BAS)
  
  
  write.csv(basins_nutri_df, paste0("/home/kuempel/Fiji_sewage/output_data/Nutrients/buffers/N_P_", name,"_per_basin_",num,"_buff.csv"))

  }

pbmclapply(buffer_list, buffer_N_P, N_cell_septic, P_cell_septic, "septic", mc.cores = 1, mc.style = "ETA")

pbmclapply(buffer_list, buffer_N_P, N_cell_direct, P_cell_direct, "direct", mc.cores = 1, mc.style = "ETA")


# Add in no buffer value

pop_ras_proj<-raster(here("projected_data/Population/population_fji_proj.tif"))

stp<-st_read(here("raw_data/Sewage_plants/Fiji Waste Water Catchment.shp")) %>% 
  mutate(STP = ifelse(TREATMENT_ == "NATABUA WWTP", "Natabua WWTP",
                      ifelse(TREATMENT_ == "Navakai PS", "Navakai WWTP", as.character(TREATMENT_)))) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(STP) %>% 
  summarise(area = sum(area, na.rm = T)) %>% 
  filter(!is.na(STP) == T,
         !STP %in% c("ACS WWTP", "Naboro WWTP", "Wailada WWTP")) 

stp<-stp %>% 
  mutate(stp_pop = extract(pop_ras_proj, ., fun = sum, na.rm = T),
         ID = 1:nrow(stp))


#stp_ras<-rasterize(stp, pop_ras_proj, field = "ID")
#writeRaster(stp_ras, here("output_data/Sewage_plants/stp_raster_residential.tif"), overwrite = T)
stp_ras<-raster(here("output_data/Sewage_plants/stp_raster_residential.tif"))

stp_ras2<-stp_ras
stp_ras2[stp_ras2>0]<-0
stp_ras2[is.na(stp_ras2) == T]<-1

N_cell_septic_all<-N_cell_septic * stp_ras2
P_cell_septic_all<-P_cell_septic * stp_ras2
N_cell_direct_all<-N_cell_direct * stp_ras2
P_cell_direct_all<-P_cell_direct * stp_ras2



septic_basin_N<-raster::extract(N_cell_septic_all, fiji_basins12, fun=sum, na.rm=TRUE, df=TRUE)
septic_basin_P<-raster::extract(P_cell_septic_all, fiji_basins12, fun=sum, na.rm=TRUE, df=TRUE)

direct_basin_N<-raster::extract(N_cell_direct_all, fiji_basins12, fun=sum, na.rm=TRUE, df=TRUE)
direct_basin_P<-raster::extract(P_cell_direct_all, fiji_basins12, fun=sum, na.rm=TRUE, df=TRUE)

septic_basins_nutri_df<-full_join(septic_basin_N, fiji_basins12, by = "ID") %>% 
  rename(N = layer) %>% 
  full_join(., septic_basin_P, by = "ID") %>% 
  rename(P = layer) %>% 
  dplyr::select(ID, N, P, MAIN_BAS)

direct_basins_nutri_df<-full_join(direct_basin_N, fiji_basins12, by = "ID") %>% 
  rename(N = layer) %>% 
  full_join(., direct_basin_P, by = "ID") %>% 
  rename(P = layer) %>% 
  dplyr::select(ID, N, P, MAIN_BAS)

# Name with large number even though it is technically no buffer

write.csv(septic_basins_nutri_df, "/home/kuempel/Fiji_sewage/output_data/Nutrients/buffers/N_P_septic_per_basin_no_buff.csv")

write.csv(direct_basins_nutri_df, "/home/kuempel/Fiji_sewage/output_data/Nutrients/buffers/N_P_direct_per_basin_no_buff.csv")







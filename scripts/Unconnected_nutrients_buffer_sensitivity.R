library(tidyverse)
library(sf)
#library(startR)
library(here)
library(gdalUtils)
library(tabularaster)
library(units)
library(rgdal)
library(lwgeom)
library(pbmcapply)
library(fasterize)
#devtools::install_github("rspatial/raster")
library(raster)

t_crs<-"+proj=aea +lat_1=-9 +lat_2=-26 +lat_0=-17 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"



# STP watersheds

stp<-st_read(here("raw_data/Sewage_plants/Fiji Waste Water Catchment.shp")) %>% 
  st_transform(., crs = t_crs) %>% 
  mutate(STP = ifelse(TREATMENT_ == "NATABUA WWTP", "Natabua WWTP",
                      ifelse(TREATMENT_ == "Navakai PS", "Navakai WWTP", as.character(TREATMENT_)))) %>% 
  mutate(area = st_area(.)) %>% 
  group_by(STP) %>% 
  summarise(area = sum(area, na.rm = T)) %>% 
  filter(!is.na(STP) == T,
         !STP %in% c("ACS WWTP", "Wailada WWTP", "Naboro WWTP")) 


# River data

rivers_proj<-st_read(here("projected_data/Rivers/Fiji_rivers.shp"))
rivers_proj<-st_transform(rivers_proj, crs=t_crs)
st_crs(rivers_proj)<-t_crs

# N and P pollution

N_unconnected<-raster(here("output_data/Nutrients/N_unconnected.tif"))

# Buffer river

#  Could also buffer coastline?
  
#  Could also do something like Steph Borelle's paper on plastics:

#https://science.sciencemag.org/content/suppl/2020/09/16/369.6510.1515.DC1


buff<-seq(from = 10, to = 50, by = 10) %>% as.list()

create_buff<-function(x){
  
buffer_rivers<-st_buffer(rivers_proj, dist = x) %>% 
  mutate(ID = 1) %>% 
  st_cast(., "POLYGON")

buffer_ras<-fasterize(buffer_rivers, N_unconnected)

writeRaster(buffer_ras, paste0("/home/kuempel/Fiji_sewage/output_data/Rivers/buffers/buffer_rivers_",x,".tif"),overwrite = T)

}

pbmclapply(buff, create_buff, mc.cores = 5, mc.style = "ETA")


rm(list=ls())
gc()

t_crs<-"+proj=aea +lat_1=-9 +lat_2=-26 +lat_0=-17 +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


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
  mutate(ID = seq(from = 1, to =nrow(fiji_basins12), 1)) %>% 
  st_transform(., crs = t_crs)


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

pbmclapply(buffer_list, buffer_N_P, N_cell_septic, P_cell_septic, "septic", mc.cores = 2, mc.style = "ETA")

pbmclapply(buffer_list, buffer_N_P, N_cell_direct, P_cell_direct, "direct", mc.cores = 2, mc.style = "ETA")


# Add in no buffer value
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



N_P_buffers<-list.files(here("output_data/Nutrients/buffers"), pattern = "septic", full.names = T)

leach<-seq(0.1,1,by = 0.1)
df2<-c()
df3<-c()

for(i in 1:length(N_P_buffers)){
  name<-N_P_buffers[i]
  buff<-parse_number(name)
  c<-read.csv(N_P_buffers[i])
  n<-sum(c$N, na.rm = T)
  p<-sum(c$P, na.rm = T)
  
  for(k in 1:length(leach)){
    rm(df)
    print(leach[k])
    n_leach<-n*leach[k]
    p_leach<-p*leach[k]
    
    df<-data.frame(buffer = buff, leach_rate = leach[k], N = n_leach, P = p_leach)
    
    df2<-rbind(df2, df)
  }
  
  df3<-rbind(df3,df2)
  
}

df3<-df3 %>% 
  mutate(type = "septic")

N_P_buffers_direct<-list.files(here("output_data/Nutrients/buffers"), pattern = "direct", full.names = T)

df4<-c()
df5<-c()

for(i in 1:length(N_P_buffers_direct)){
  name<-N_P_buffers_direct[i]
  buff<-parse_number(name)
  c<-read.csv(N_P_buffers_direct[i])
  n<-sum(c$N, na.rm = T)
  p<-sum(c$P, na.rm = T)
  
  for(k in 1:length(leach)){
    rm(df)
    print(leach[k])
    n_leach<-n*leach[k]
    p_leach<-p*leach[k]
    
    df<-data.frame(buffer = buff, leach_rate = leach[k], N = n_leach, P = p_leach)
    
    df4<-rbind(df4, df)
  }
  
  df5<-rbind(df5,df4)
  
}

df5<- df5 %>% 
  mutate(type = "direct")

done<-rbind(df3, df5)

write.csv(done, here("output_data/Sensitivity/Leaching_buffer_sensitivity.csv"))


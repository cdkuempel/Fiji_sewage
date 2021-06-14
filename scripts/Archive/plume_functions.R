#pp_sf<-point_list[[1]]
#buffer_value<-80000
#resolution_raster<-c(1000,1000)
#searchdistance<-90000

N_plume<-function(pp_sf,#list of pour points
                     ocean_raster,
                     buffer_value,#buffer from pour point
                     resolution_raster,#in meters
                     searchdistance
){
  #set initial sediment value
  initial.value<-pp_sf@data$N
  #create buffer around point and clip ocean raster
  point<-st_as_sf(pp_sf,coords=c("longitude.y","latitude.y"))#first need to convert to sf object
  point_buff<-st_buffer(point, buffer_value)#
  r1<-crop(ocean_raster,extent(point_buff),snap = "near")
  
  #create raster to match point location to ocean raster
  #for a new version we can use the function nearest.raster.point, but I had already solved this issue 
  #rasterizing the point
  
  point_1000<-st_buffer(point, 500)#
  #create raster of pour point
  r2  <- raster(ncol=1, nrow=1, xmn=extent(point_1000)[1], xmx=extent(point_1000)[2], ymn=extent(point_1000)[3], ymx=extent(point_1000)[4], res = resolution_raster)
  crs(r2)<-crs(ocean_raster)
  r2[]<-initial.value
  
  #now sum to match rasters of different extent
  
  extend_all =
    function(rasters){
      extent(Reduce(extend,rasters))
    }
  
  sum_all =
    function(rasters, extent){
      re = lapply(rasters, function(r){extend(r,extent)})
      Reduce("+", re)
    }
  
  rasterOptions(tolerance = 1)
  r3 = sum_all(list(r1,r2), extend_all(list(r1,r2)))
  
  #identify if points are on the coast
  
  if(length(which(r3[]>0))< 1){
    snap_point<-nearestLand(coordinates(r2)[],r1,searchdistance)
    pp_sf@coords<-snap_point
    point<-st_as_sf(pp_sf,coords=c("longitude.y","latitude.y"))
    point_1000<-st_buffer(point, 500)#
    #create raster of pour point
    r2  <- raster(ncol=1, nrow=1, xmn=extent(point_1000)[1], xmx=extent(point_1000)[2], ymn=extent(point_1000)[3], ymx=extent(point_1000)[4], res = resolution_raster)
    crs(r2)<-crs(ocean_raster)
    r2[]<-initial.value
    r3 = sum_all(list(r1,r2), extend_all(list(r1,r2)))
  }
  
  #find exact pour point position
  pour_point <- which(r3[] >0)
  
  fromCoords<-pp_sf@coords
  
  raster_temp<-r3
  r<-raster_temp
  
  #create transition layer to be used in cost layer
  transition_layer<-transition(r1,transitionFunction=mean, directions=8)
  
  transition_layer<-geoCorrection(transition_layer,scl=FALSE)
  
  cost_layer<-accCost(transition_layer, fromCoords)
  #calculate distance frequencies
  
  
  initial.value<-pp_sf@data$N
  
  values<-rasterToPoints(cost_layer)[,3]
  frequency_distances<-as.data.frame(table(round(values,1)))
  frequency_distances<-subset(frequency_distances,frequency_distances[1] != Inf)
  frequency_distances[,1]<-as.character(frequency_distances[,1])
  frequency_distances[,1]<-as.numeric(frequency_distances[,1])
  new.values<-initial.value
  raster_temp[pour_point]<-initial.value*0.005  ##Decay parameter
  
  for(i in 2:nrow(frequency_distances)){
    cellstoreplace<-which(round(cost_layer[],1)==frequency_distances[i,1])
    percell = initial.value * 0.005  ## decay parameter
    sumcells = percell * length(cellstoreplace)
    remain = initial.value - sumcells
    initial.value = remain
    raster_temp[cellstoreplace]<-percell
    if(remain < pp_sf@data$N*0.0005) {
      break
    }
  }
  raster_temp[which(raster_temp[]==0)]<-NA
  if(length(unique(raster_temp[]))>1){
    raster_temp<-trim(raster_temp)
  }
  return(raster_temp)
}



P_plume<-function(pp_sf,#list of pour points
                     ocean_raster,
                     buffer_value,#buffer from pour point
                     resolution_raster,#in meters
                     searchdistance #max searching distance to snap the points to the coast
                     
){
  #set initial sediment value
  initial.value<-pp_sf@data$P
  #create buffer around point and clip ocean raster
  point<-st_as_sf(pp_sf,coords=c("longitude.y","latitude.y"))#first need to convert to sf object
  point_buff<-st_buffer(point, buffer_value)#
  r1<-crop(ocean_raster,extent(point_buff),snap = "near")
  
  #create raster to match point location to ocean raster
  #for a new version we can use the function nearest.raster.point, but I had already solved this issue 
  #rasterizing the point
  
  point_1000<-st_buffer(point, 500)#
  #create raster of pour point
  r2  <- raster(ncol=1, nrow=1, xmn=extent(point_1000)[1], xmx=extent(point_1000)[2], ymn=extent(point_1000)[3], ymx=extent(point_1000)[4])
  crs(r2)<-crs(ocean_raster)
  r2[]<-initial.value
  
  #now sum to match rasters of different extent
  
  extend_all =
    function(rasters){
      extent(Reduce(extend,rasters))
    }
  
  sum_all =
    function(rasters, extent){
      re = lapply(rasters, function(r){extend(r,extent)})
      Reduce("+", re)
    }
  
  rasterOptions(tolerance = 1)
  r3 = sum_all(list(r1,r2), extend_all(list(r1,r2)))
  
  #identify if points are on the coast
  
  if(length(which(r3[]>0))< 1){
    snap_point<-nearestLand(coordinates(r2)[],r1,searchdistance)
    pp_sf@coords<-snap_point
    point<-st_as_sf(pp_sf,coords=c("longitude.y","latitude.y"))
    point_1000<-st_buffer(point, 500)#
    #create raster of pour point
    r2  <- raster(ncol=1, nrow=1, xmn=extent(point_1000)[1], xmx=extent(point_1000)[2], ymn=extent(point_1000)[3], ymx=extent(point_1000)[4])
    crs(r2)<-crs(ocean_raster)
    r2[]<-initial.value
    r3 = sum_all(list(r1,r2), extend_all(list(r1,r2)))
  }
  
  #find exact pour point position
  pour_point <- which(r3[] >0)
  
  fromCoords<-pp_sf@coords
  
  raster_temp<-r3
  r<-raster_temp
  
  #create transition layer to be used in cost layer
  transition_layer<-transition(r1,transitionFunction=mean, directions=8)
  
  transition_layer<-geoCorrection(transition_layer,scl=FALSE)
  
  cost_layer<-accCost(transition_layer, fromCoords)
  #calculate distance frequencies
  
  
  initial.value<-pp_sf@data$P
  
  values<-rasterToPoints(cost_layer)[,3]
  frequency_distances<-as.data.frame(table(round(values,1)))
  frequency_distances<-subset(frequency_distances,frequency_distances[1] != Inf)
  frequency_distances[,1]<-as.character(frequency_distances[,1])
  frequency_distances[,1]<-as.numeric(frequency_distances[,1])
  new.values<-initial.value
  raster_temp[pour_point]<-initial.value*0.005  ##Decay parameter
  
  for(i in 2:nrow(frequency_distances))
  {
    cellstoreplace<-which(round(cost_layer[],1)==frequency_distances[i,1])
    percell = initial.value * 0.005  ## decay parameter
    sumcells = percell * length(cellstoreplace)
    remain = initial.value - sumcells
    initial.value = remain
    raster_temp[cellstoreplace]<-percell
    if(remain < pp_sf@data$N*0.0005) {
      break
    }
  }
  raster_temp[which(raster_temp[]==0)]<-NA
  if(length(unique(raster_temp[]))>1){
    raster_temp<-trim(raster_temp)
  }
  return(raster_temp)
}




coral_func<-function(x){ # Y is plume data, X is coral list
  sub<-x[1]
  if(is.null(intersect(extent(all), extent(sub)))){ #57 is FALSE
    focus_corals<-x[1]
    #focus_corals<-st_cast(st_as_sf(focus_corals),"POLYGON")
    focus_corals$Fiji_connected_N_plumes<-0
    focus_corals$Fiji_connected_P_plumes<-0
    focus_corals$Fiji_septic_N_plumes<-0
    focus_corals$Fiji_septic_P_plumes<-0
    focus_corals$Fiji_direct_N_plumes<-0
    focus_corals$Fiji_direct_P_plumes<-0
    focus_corals$Fiji_connected_tourism_N_plumes<-0
    focus_corals$Fiji_connected_tourism_P_plumes<-0
  }else{
    focus_corals<-x[1]
    #focus_corals<-st_cast(st_as_sf(focus_corals),"POLYGON")
    #focus_corals$ID<-1:nrow(focus_corals)
    #  done<-c()
    # for(i in 1:nrow(focus_corals)){
    #    sub_coral<-focus_corals[i,]
    if(is.null(intersect(extent(all), extent(focus_corals)))){
      focus_corals$Fiji_connected_N_plumes<-0
      focus_corals$Fiji_connected_P_plumes<-0
      focus_corals$Fiji_septic_N_plumes<-0
      focus_corals$Fiji_septic_P_plumes<-0
      focus_corals$Fiji_direct_N_plumes<-0
      focus_corals$Fiji_direct_P_plumes<-0
      focus_corals$Fiji_connected_tourism_N_plumes<-0
      focus_corals$Fiji_connected_tourism_P_plumes<-0
    }else{
      raster.polygon<-crop(all,extent(focus_corals),snap = "near")
      nutrients<-cellStats(raster.polygon,"mean") %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column(., "Type") %>% 
        rename(Value = ".") %>% 
        pivot_wider(., names_from = Type, values_from = Value)
      
      focus_corals<-cbind(focus_corals, nutrients)
    }
    
    #done<-rbind(done, focus_corals)
    #focus_corals<-done
  }
  focus_corals
}

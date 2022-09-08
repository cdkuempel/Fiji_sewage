sed_flow<-function(#shapefile of pour points
                   pp_sf,
                   #the column from the attribute table of your pour
                   #point shapefile that has the 
                   #sediment load initial value
                   initial.value,
                   #an ocean raster where the coast has NA values
                   ocean_raster,
                   #buffer from  point to generate a sediment plume
                   #(basically the raster template)
                   buffer_value,
                   resolution_raster,#in km
                   #max searching distance to snap the points to the coast
                   #this is when a pour point does no allign with the coastline from the
                   #ocean raster
                   searchdistance, 
                   decay_percent,
                   remaining_left
                   
){
  #set initial sediment value
  #initial.value<-pp_sf@data$sed_per_area
  #create buffer around point and clip ocean raster
  #point<-st_as_sf(pp_sf)#first need to convert to sf object
  point_buff<-st_buffer(st_as_sf(pp_sf), buffer_value)#
  r1<-crop(ocean_raster,extent(point_buff),snap = "near")
  #plot(r1)
  #check whether is an sf or sp object and extract coordinates
  if(class(pp_sf)[1]=="sf"){
    pp_sf<-as(pp_sf,"Spatial")
  }
  
  #create raster to match point location to ocean raster
  #for a new version we can use the function nearest.raster.point, but I had already solved this issue 
  #rasterizing the point
  
  point_1000<-st_buffer(st_as_sf(pp_sf), 500)#
  #plot(point_1000, add = T)
  #create raster of pour point
  r2  <- raster(ncol=1, nrow=1, xmn=extent(point_1000)[1], 
                xmx=extent(point_1000)[2], ymn=extent(point_1000)[3], 
                ymx=extent(point_1000)[4])
  crs(r2)<-crs(ocean_raster)
  r2[]<-initial.value
  
  #now sum to match rasters of different extent
  r2e<-extend(r2,r1)
  r2e = resample(r2e, r1, "ngb")
  r3 = stack(r1,r2e)
  r3 <- calc(r3, sum)
  #identify if points are on the coast####
  #use function nearestLand from package seegSDM
  nearestLand <- function (points, raster, max_distance) {
    nearest <- function(lis, raster) {
      neighbours <- matrix(lis[[1]], ncol = 2)
      point <- lis[[2]]
      land <- !is.na(neighbours[, 2])
      if (!any(land)) {
        return(c(NA, NA))
      }
      else {
        coords <- xyFromCell(raster, neighbours[land, 1])
        if (nrow(coords) == 1) {
          return(coords[1, ])
        }
        dists <- sqrt((coords[, 1] - point[1])^2 + (coords[, 
                                                           2] - point[2])^2)
        return(coords[which.min(dists), ])
      }
    }
    neighbour_list <- raster::extract(raster, points, buffer = max_distance, 
                              cellnumbers = TRUE)
    neighbour_list <- lapply(1:nrow(points), function(i) {
      list(neighbours = neighbour_list[[i]], point = as.numeric(points[i, 
                                                                       ]))
    })
    return(t(sapply(neighbour_list, nearest, raster)))
  }
  
  #check if point is in the coast####
  
  if(length(which(r3[]>0))< 1){
    snap_point<-nearestLand(coordinates(r2)[],r1,searchdistance)
    pp_sf@coords<-snap_point
    #check whether is an sf or sp object and extract coordinates
    point_1000<-st_buffer(st_as_sf(pp_sf), 500)#
    #plot(point_1000,add = T)
    #create raster of pour point
    r2  <- raster(ncol=1, nrow=1, xmn=extent(point_1000)[1], xmx=extent(point_1000)[2], ymn=extent(point_1000)[3], ymx=extent(point_1000)[4])
    crs(r2)<-crs(ocean_raster)
    r2[]<-initial.value
    r2<-extend(r2,r1)
    r2e = resample(r2e, r1, "ngb")
    r3 = stack(r1,r2e)
    r3 <- calc(r3, sum)
  }
  #find exact pour point position
  pour_point <- which(r3[] >0)
  
  #raster template to save plume
  raster_temp<-r3
  
  #create transition layer to be used in cost layer
  transition_layer<-transition(r1,transitionFunction=mean, directions=8)
  
  #plot(raster(transition_layer))
  transition_layer<-geoCorrection(transition_layer,scl=FALSE)
  
  fromCoords<-pp_sf@coords
  #points(fromCoords)

  cost_layer<-accCost(transition_layer, fromCoords)
  #plot(cost_layer)
  
  #extract distance frequencies in the cost layer
  
  values<-rasterToPoints(cost_layer)[,3]
  frequency_distances<-as.data.frame(table(round(values,1)))
  frequency_distances<-subset(frequency_distances,frequency_distances[1] != Inf)
  frequency_distances[,1]<-as.character(frequency_distances[,1])
  frequency_distances[,1]<-as.numeric(frequency_distances[,1])
  remain<-initial.value
  raster_temp[pour_point]<-initial.value*decay_percent  ##Decay parameter
  
  for(i in 2:nrow(frequency_distances))
  {
    cellstoreplace<-which(round(cost_layer[],1)==frequency_distances[i,1])
    percell = remain * decay_percent  ## decay parameter
    sumcells = percell * length(cellstoreplace)
    remain = remain - sumcells
    #new.values = remain
    raster_temp[cellstoreplace]<-percell
    #plot(raster_temp)
    if(remain < initial.value*remaining_left) {
      break
    }
  }
  raster_temp[which(raster_temp[]==0)]<-NA
  if(length(unique(raster_temp[]))>1){
    raster_temp<-trim(raster_temp)
  }
  #plot(raster_temp)
  return(raster_temp)
}

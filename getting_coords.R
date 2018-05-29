library(rgdal)
library(ggmap)

getting_coords <- function(locations){
  locations_coords <- data.frame(name = character(0), lon = numeric(0L), lat = numeric(0L))
  
  for(i in 1:length(locations)){
    print(locations[i])
    geocode <- geocode(locations[i])
    if(is.na(geocode$lon) & is.na(geocode$lat)){
      while(is.na(geocode$lon) & is.na(geocode$lat)){
        #this aims to tackle the submission limit set by the api
        Sys.sleep(2)
        geocode <- geocode(locations[i])
      }
    }
    locations_coords_row <- data.frame(name = locations[i], lon = geocode$lon, lat = geocode$lat)
    locations_coords <- rbind(locations_coords, locations_coords_row)
  }
  
  coords <- locations_coords[,c("lon","lat")]
  
  #get projection used in ONS data
  #proj4string(nhsengland)
  #should be:
  #"+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36 +units=m +no_defs +ellps=airy +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894"
  
  #add in a baseline projection - this, I think, is the most basic projection
  locations_coords_spdf <- SpatialPointsDataFrame(coords = coords, data = locations_coords, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  #then transform the coords into the system used in the ONS data
  locations_coords_transformed <- spTransform(locations_coords_spdf, CRS(proj4string(nhsengland)))
  
  #convert to dataframe for use in ggplot
  locations_coords_transformed_df <- as.data.frame(locations_coords_transformed)
  
  return(locations_coords_transformed_df)
}
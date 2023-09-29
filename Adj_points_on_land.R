###############################################################################
#            Identifying / adjusting / removing locations on land
###############################################################################


rm(list=ls())

library(geosphere)
library(maptools)
library(sp)
library(raster)
library(rgdal)
library(ggplot2)
library(mapview)
library(rworldmap)
#library(rgdal)
data("wrld_simpl")


#------------------------------ LOAD AND FORMAT DATA --------------------------


crds_temp <- read.csv(file = "crds_all_for_UD.csv", header = TRUE, sep = ",", 
                      na = "NA")
crds_temp <- subset(crds_temp, !is.na(Lon_smooth))
crds_temp$ID_GLS <- crds_temp$ID
crds_temp$ID_GLS <- gsub("_1","",crds_temp$ID_GLS)
crds_temp$ID_GLS <- gsub("_2","",crds_temp$ID_GLS)
GLS_index <- read.csv("Abbotts_GLS_Index.csv", header = TRUE)
GLS_index$Animal_ID_ABBBS <- as.numeric(GLS_index$Animal_ID_ABBBS)

# Add animal ID and sex to the crds file
for(i in 1:nrow(crds_temp)){
  ID_temp <- crds_temp$ID_GLS[i]
  animal_ID <- GLS_index$Animal_ID_ABBBS[GLS_index$ID==ID_temp]
  crds_temp$animal_ID[i] <- as.character(animal_ID)
}

# Format data as sp object
data<-data.frame(ID = as.vector(crds_temp$animal_ID), X = as.vector(crds_temp$x), 
                 Y = as.vector(crds_temp$y))

spdata <- SpatialPointsDataFrame(coords = cbind(data$X, data$Y), 
                                 data = subset(data, select = ID),
                                 proj4string = CRS("+init=epsg:3857"))

# Check data
head(spdata)
mapview(spdata)


#---------------------- CREATE A MAP OF AREA OF INTEREST ----------------------


# Load world polygon data
world <- getMap(resolution = "less islands")

# Define extent
clipper <- as(extent(50, 170, -50, 50), "SpatialPolygons")

# Set CRS of extent to match world map
proj4string(clipper) <- CRS(proj4string(world))

# Subset world map to extent
world_clip <- raster::intersect(world, clipper)

# Visualise clipped map
plot(world_clip)


#-------------------------- IDENTIFY POINTS ON LAND ---------------------------


# Re project coordinates to match locations projection
land_mask <- spTransform(world_clip, CRS("+init=epsg:3857"))

# Identify points on land
pts_on_land <- over(spdata, land_mask, fn = NULL) 

# Check results
summary(pts_on_land)

# Create vector of 0/1 for points off/on land
pts_on_land_vector <- ifelse(is.na(pts_on_land$LON),0,1)

# Add the vector to the coordinates file
crds_temp$pts_on_land <- pts_on_land_vector


#------------- CALCULATE DISTANCE TO NEAREST POINT ON WORLD POLYGON -----------
#-------------------- AND ADD CORRESPONDING COORDINATES -----------------------


# Create an index of point on land in crds file
index_pts_on_land <- which(crds_temp$pts_on_land==1)

# Merging polygons (needed for proper later use of dist2Line function)
Simple_land <-  unionSpatialPolygons(world_clip, world_clip@data$REGION)

# Comparing the two maps
plot(world_clip)
plot(Simple_land)

# Subset crds file for points on land
temp <- crds_temp[index_pts_on_land, c("Lon_smooth","Lat_interp")]

# Calculate distance to line & retrieve corresponding coordinates
d <- as.data.frame(dist2Line(temp, Simple_land))

# Add distance and new lat/lon to the crds file based on pts on land incex
crds_temp$distance_inland <- NA
crds_temp$Lon_adj <- crds_temp$Lon_smooth
crds_temp$Lat_adj <- crds_temp$Lat_interp
crds_temp[index_pts_on_land,"distance_inland"] <- d$distance/1000
crds_temp[index_pts_on_land,"Lon_adj"] <- d$lon
crds_temp[index_pts_on_land,"Lat_adj"] <- d$lat


# Check distance data for potential points to remove (>180km)
summary(crds_temp$distance_inland)

# Plot new points
plot(Simple_land)
points(cbind(d$lon,d$lat), col='blue', pch=20)

# Convert coordinates to xy
proj_crds <- project(as.matrix(crds_temp[,17:18]), "+init=epsg:3857")

# Add converted coordinates to file
crds_temp$x_adj <- proj_crds[,1]
crds_temp$y_adj <- proj_crds[,2]

# Visualise all data with adjusted xy
spdata2 <- SpatialPointsDataFrame(coords = cbind(crds_temp$x_adj, crds_temp$y_adj), 
                                  data = subset(crds_temp, select = ID),
                                  proj4string = CRS("+init=epsg:3857"))
mapview(spdata2)


# Save file with adjusted coordinates
write.csv(crds_temp, file=paste0("crds_all_pts_land_adj_with_CI.csv"), 
          row.names = FALSE, quote=FALSE)



















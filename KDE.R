###############################################################################
#                           UTILISATION DISTRIBUTIONS
#                        FROM KERNEL DENSITY ESTIMATION
###############################################################################


rm(list=ls())

library(geosphere)
library(rgeos)
library(sp)
library(rworldmap)
library(raster)
library(rgdal)
library(adehabitatHR)
library(mapview)
library(dplyr)
library(tictoc)


############################### LOAD / FORMAT DATA ############################

# Load data
crds_temp <- read.csv(file = "crds_all_pts_land_adj.4.csv", header = TRUE, 
                      sep = ",", na = "NA") # Points on CI already removed
GLS_index <- read.csv("Abbotts_GLS_Index.csv", header = TRUE)

# Remove incomplete trips
crds_temp <- crds_temp[!crds_temp$ID %in% c("V905026","22584"),]

# Format data
data<-data.frame(ID = as.vector(crds_temp$animal_ID), X = as.vector(crds_temp$x_adj), 
                 Y = as.vector(crds_temp$y_adj), timestamp = paste0(
                   crds_temp$Date, " ", crds_temp$Time), 
                 sex=as.vector(crds_temp$sex), trip_ID = as.vector(crds_temp$ID))
data$timestamp <- as.POSIXct(data$timestamp)
summary(data)
head(data)

spdata <- SpatialPointsDataFrame(coords = cbind(data$X, data$Y), 
                                 data = subset(data, select = c(ID, sex, trip_ID, 
                                                                timestamp)),
                                 proj4string = CRS("+init=epsg:3857"))

# Check columns of the data to use in kernelUD function
head(spdata)
plot(spdata)


############################ CREATE 50*50km GRID ##############################


# Define 50*50km grid
cell_size <- 50*1000
x <- seq(xmin(spdata)-(50*cell_size),xmax(spdata)+(50*cell_size),by=cell_size)
y <- seq(ymin(spdata)-(50*cell_size),ymax(spdata)+(50*cell_size),by=cell_size)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy) 
projection(xy) <- CRS("+init=epsg:3857")

# View grid and data points
plot(xy)
plot(spdata, add=TRUE)

rm(x,y)

#################### IDENTIFY SINGLE/MULTIPLE TRIPS PER INDIV #################


ID_index <- as.data.frame(unique(crds_temp$ID))
colnames(ID_index)[1]<-"ID"
ID_index$ID_GLS <- ID_index$ID
ID_index$ID_GLS <- gsub("_1","",ID_index$ID_GLS)
ID_index$ID_GLS <- gsub("_2","",ID_index$ID_GLS)
for(i in 1:nrow(ID_index)){
  ID_temp <- ID_index$ID_GLS[i]
  ID_index$ID_animal[i] <- GLS_index$Animal_ID_ABBBS[GLS_index$ID==ID_temp]
}

n_trips <- count(ID_index, ID_animal)
multi_trips <- n_trips[!n_trips$n==1,]
single_trips <- n_trips[n_trips$n==1,]

multi_trips
single_trips


###############################################################################
#                          KDE w/ SAME h FOR EACH TRIP
###############################################################################


smooth_h <- 186*1000

#----------------- CALCULATE FOR INDIV WITH MULTIPLE TRIPS --------------------

multi_trips
multi_trips_list <- vector(mode = "list", length = nrow(multi_trips))
for(i in 1:nrow(multi_trips)){
  ID_animal_temp <- multi_trips$ID_animal[i]
  multi_trip_indiv_temp <- spdata[spdata$ID==ID_animal_temp,]
  # Calculate kernel
  k_multi <- kernelUD(multi_trip_indiv_temp[,3], grid=xy, extent=1, h=smooth_h)
  # Mean of values per cell
  k_multi_temp <- stack(lapply(k_multi, raster))
  udm <- mean(k_multi_temp)
  udm@data@names <- ID_animal_temp
  # Add raster to list
  multi_trips_list[[i]] <- udm
}

#------------------- CALCULATE FOR INDIV WITH SINGLE TRIPS --------------------

single_trips
single_trips_list <- vector(mode = "list", length = nrow(single_trips))
for(i in 1:nrow(single_trips)){
  ID_animal_temp <- single_trips$ID_animal[i]
  single_trip_temp <- spdata[spdata$ID==ID_animal_temp,]
  # Calculate kernel
  k_single <- kernelUD(single_trip_temp[,3], grid=xy, extent=1, h=smooth_h)
  #image(k_single)
  # Mean of values per cell
  k_single_temp <- stack(lapply(k_single, raster))
  udm <- mean(k_single_temp)
  udm@data@names <- ID_animal_temp
  # Add raster to list
  single_trips_list[[i]] <- udm
}

# ------------------------------ MERGE UDs ------------------------------------

# Combine list of rasters
list_combined <- c(multi_trips_list, single_trips_list)
# Stack the rasters
combined_temp <- stack(list_combined)
# Mean of values per cell
udm_all <- mean(combined_temp)

# create a single estUD object (assigning same ID to all records)
spdata_temp <- spdata
spdata_temp$ID <- "1" 
k_all <- kernelUD(spdata_temp[,1], grid=xy, extent=1, h=186*1000)
image(k_all)

# Replace the pixels values of the estUD object with the mean values from the raster
k_all[[1]]@data$ud <- udm_all@data@values[k_all[[1]]@grid.index] 
#image(k_all)

area_all <- kernel.area(k_all, percent = c(50,95,99.9), unin="m", unout="km2", 
                        standardize = FALSE)

# Calculate isopleths
pol99.9_all <- getverticeshr(k_all,99.9)
pol95_all <- getverticeshr(k_all,95)
pol50_all <- getverticeshr(k_all,50)

mapview(pol99.9_all)+mapview(pol95_all)+mapview(pol50_all)

rm(list_combined, combined_temp, udm_all, udm, k_single_temp, k_single, 
   ID_animal_temp, i, single_trips_list, multi_trips_list, multi_trip_indiv_temp, 
   k_multi, k_multi_temp)

#--------------- Identify overlap with land masses ----------------------------


# Get vertices for 99.9% UD
pol99.9_all <- getverticeshr(k_all,99.9)

# Get world polygon 
world <- getMap(resolution = "low")
# Subset area of interest
clipper <- as(extent(50, 170, -50, 50), "SpatialPolygons")
proj4string(clipper) <- CRS(proj4string(world))
world_clip <- raster::intersect(world, clipper)
# Reproject land mask and UD
land_mask <- spTransform(world_clip, CRS("+init=epsg:3857"))
UD_reproj <- spTransform(pol99.9_all, CRS("+init=epsg:3857"))
# View result
mapview(land_mask) + mapview(pol99.9_all)

# Retrieve parts of UD on land
int<-gIntersection(UD_reproj,land_mask,byid=F)
mapview(int)

# Calculate area on land and proportion of 99.9% UD
areas_overlap <- NULL
for(i in 1:length(int@polygons[[1]]@Polygons)){
  areas_overlap <- c(areas_overlap, int@polygons[[1]]@Polygons[[i]]@area)
}
overlap_area <- sum(areas_overlap)*1.0E-6
UD_area <- kernel.area(k_all, 99.9, unin = "m", unout = "km2", standardize = FALSE)
overlap_area/UD_area

rm(overlap_area, UD_area, int, UD_reproj, land_mask, clipper, world_clip, 
   areas_overlap, world)


#-------------- Calculate points distance from 99.9% UD contour ---------------


# Identify points in countour line
contour_temp <- as(pol99.9_all, "SpatialPolygons")
pts_in_99.9 <- over(spdata, contour_temp, fn = NULL) 
summary(pts_in_99.9)

# Calculate distance to line
temp_pol_coord <- UD_reproj@polygons[[1]]@Polygons[[1]]@coords
temp_pol <- project(temp_pol_coord, "+init=epsg:3857", inv=TRUE)
crds_dist <- crds_temp[,c("Lon_adj","Lat_adj")]
d <- as.data.frame(dist2Line(crds_dist, temp_pol))

# Check proportion of loc within 400km from the 99.9% UD contour line
within_400 <- d$distance[d$distance<400000]
length(within_400)/nrow(crds_temp)

rm(pts_in_99.9, contour_temp, temp_pol_coord, temp_pol, crds_dist)


###############################################################################
#                                 KDE PER SEX
###############################################################################


smooth_h <- 186*1000
sex_list <- c("F","M")

spdata2 <- SpatialPointsDataFrame(coords = cbind(data$X, data$Y), 
                                  data = subset(data, select = c(ID, sex)),
                                  proj4string = CRS("+init=epsg:3857"))
k_mf <- kernelUD(spdata2[,2], grid=xy, extent=1, h=smooth_h)
head(spdata2)


#---------------------------- !!!! LOOP !!!! ----------------------------------
for(i in 1:length(sex_list)){   
  sex <- sex_list[i]
  # Subset data per sex
  subset_sex <- spdata[spdata$sex==sex,]
  #----------------- CALCULATE FOR INDIV WITH MULTIPLE TRIPS ------------------
  # Add sex to list of multi trips
  for(j in 1:nrow(multi_trips)){
    ID_temp <- multi_trips$ID_animal[j]
    multi_trips$sex[j] <- GLS_index$Sex[which(GLS_index$Animal_ID_ABBBS==ID_temp)] 
  }
  # Subset list of multi trip birds per sex
  multi_trips_sex <- multi_trips[multi_trips$sex==sex,]  
  # Create empty list to store rasters
  multi_trips_list <- vector(mode = "list", length = nrow(multi_trips_sex))
  # Run KDE per bird
  for(k in 1:nrow(multi_trips_sex)){
    ID_animal_temp <- multi_trips_sex$ID_animal[k]
    multi_trip_indiv_temp <- subset_sex[subset_sex$ID==ID_animal_temp,]
    # Calculate kernel
    k_multi <- kernelUD(multi_trip_indiv_temp[,3], grid=xy, extent=1, h=smooth_h)
    # Mean of values per cell
    k_multi_temp <- stack(lapply(k_multi, raster))
    udm <- mean(k_multi_temp)
    udm@data@names <- ID_animal_temp
    # Add raster to list
    multi_trips_list[[k]] <- udm
  }
  #------------------- CALCULATE FOR INDIV WITH SINGLE TRIPS ------------------
  # Add sex to list of single trips
  for(m in 1:nrow(single_trips)){
    ID_temp <- single_trips$ID_animal[m]
    single_trips$sex[m] <- GLS_index$Sex[which(GLS_index$Animal_ID_ABBBS==ID_temp)] 
  }
  # Subset single trip birds per sex
  single_trips_sex <- single_trips[single_trips$sex==sex,]
  # Create empty list to store rasters
  single_trips_list <- vector(mode = "list", length = nrow(single_trips_sex))
  # Run KDE per bird
  for(j in 1:nrow(single_trips_sex)){
    ID_animal_temp <- single_trips_sex$ID_animal[j]
    single_trip_temp <- subset_sex[subset_sex$ID==ID_animal_temp,]
    # Calculate kernel
    k_single <- kernelUD(single_trip_temp[,3], grid=xy, extent=1, h=smooth_h)
    # Mean of values per cell
    k_single_temp <- stack(lapply(k_single, raster))
    udm <- mean(k_single_temp)
    udm@data@names <- ID_animal_temp
    # Add raster to list
    single_trips_list[[j]] <- udm
  }
  #------------- MERGE MUTLI & SINGLE LISTS + CONVERT INTO estUDm -------------
  # Combine list of rasters
  list_combined <- c(multi_trips_list, single_trips_list)
  # Stack the rasters
  combined_temp <- stack(list_combined)
  # Mean of values per cell
  plot(udm_all <- mean(combined_temp))
  # Replace values 
  k_mf[[i]]@data$ud <- udm_all@data@values[k_mf[[i]]@grid.index]
}
#------------------------ !!!! END LOOP !!!! ----------------------------------

rm(list_combined, udm, combined_temp, k_single_temp, k_single, ID_animal_temp, 
   single_trip_temp, single_trips_list, single_trips_sex, ID_temp, multi_trips_list,
   multi_trips_sex, k_multi_temp, k_multi, multi_trip_indiv_temp, sex, i, j, k, m)

# View KDE per sex
image(k_mf)

# Calculate isopleths
pol99_mf <- getverticeshr(k_mf,99)
pol50_mf <- getverticeshr(k_mf,50)
mapview(pol99_mf)+mapview(pol50_mf)

# Quantify overlap (Battacharyya's affinity)
overlap_95_mf <- kerneloverlaphr(k_mf, meth="BA", conditional=TRUE, percent=95)
overlap_50_mf <- kerneloverlaphr(k_mf, meth="BA", conditional=TRUE, percent=50)


###############################################################################
#                           KDE PER SEX RANDOMISED
###############################################################################


iterations <- 1000
smooth_h <- 186*1000
sex_list <- c("F","M")

# Load data
crds_rdm_sex <- crds_temp

sex_index <- GLS_index[,c("Animal_ID_ABBBS", "GLS_ID", "Sex")]
sex_index <- sex_index[sex_index$Animal_ID_ABBBS %in% spdata$ID,]
sex_index <- sex_index[!duplicated(sex_index$Animal_ID_ABBBS),]
head(sex_index)
unique(sex_index$Animal_ID_ABBBS)
length(sex_index$Animal_ID_ABBBS)

sex_vector <- ifelse(sex_index$Sex=="M",1,0)
sex_vector
length(sex_vector)

overlap_matrix <- NULL


#------------------------------- !!!! LOOP !!!! -------------------------------

for(l in 1:iterations){
  if(l==1){tic()}else{}
  # Print iterations to follow progress
  print(paste0(l,"/",iterations))
  
  #------------------------------ Randomise sex -------------------------------
  
  # Create a vector of random 0-1s
  rdm_sex <- sample(sex_vector,size=nrow(sex_index),replace=FALSE)
  # Add random vector to sex index
  sex_index$random <- rdm_sex
  # Assign sex based on 0-1s
  for(i in 1:nrow(sex_index)){
    if(sex_index$random[i]==1){
      sex_index$rdm_sex[i]<-"M"
    }else{sex_index$rdm_sex[i]<-"F"}
  }
  # Add the defined random sex to the crds file
  for(i in 1:nrow(crds_rdm_sex)){
    ID_animal_temp <- crds_rdm_sex$animal_ID[i]
    crds_rdm_sex$rdm_sex[i] <- 
      sex_index$rdm_sex[sex_index$Animal_ID_ABBBS==ID_animal_temp]  
  }
  
  #-------------------------- Format data for KDE -----------------------------
  
  # Reformat data
  data_rdm_sex <-data.frame(ID = as.vector(crds_rdm_sex$animal_ID), 
                            X = as.vector(crds_rdm_sex$x_adj), 
                            Y = as.vector(crds_rdm_sex$y_adj), 
                            sex=as.vector(crds_rdm_sex$rdm_sex), 
                            trip_ID = as.vector(crds_rdm_sex$ID))
  # Convert data into sp object
  spdata_rdm_sex <- SpatialPointsDataFrame(
    coords = cbind(data_rdm_sex$X, data_rdm_sex$Y), 
    data = subset(data_rdm_sex, select = c(ID, sex, trip_ID)),
    proj4string = CRS("+init=epsg:3857"))
  # Create estUDm object (will replace cell values later on)
  k_mf_rdm <- kernelUD(spdata_rdm_sex[,2], grid=xy, extent=1, h=smooth_h)
  
  #--------------------------------- KDE PER SEX ------------------------------
  
  for(i in 1:length(sex_list)){   
    sex <- sex_list[i]
    # Subset data per sex
    subset_sex <- spdata_rdm_sex[spdata_rdm_sex$sex==sex,]
    
    #----------------- CALCULATE FOR INDIV WITH MULTIPLE TRIPS ----------------
    
    # Add sex to list of multi trips
    for(j in 1:nrow(multi_trips)){
      ID_temp <- multi_trips$ID_animal[j]
      multi_trips$sex[j] <- spdata_rdm_sex$sex[which(spdata_rdm_sex$ID==ID_temp)] 
    }
    # Subset list of multi trip birds per sex
    multi_trips_sex <- multi_trips[multi_trips$sex==sex,]
    # Create empty list to store rasters
    multi_trips_list <- vector(mode = "list", length = nrow(multi_trips_sex))
    # Run KDE only if any bird with multiple trips
    if(nrow(multi_trips_sex)>=1){
      # Run KDE per bird
      for(k in 1:nrow(multi_trips_sex)){
        ID_animal_temp <- multi_trips_sex$ID_animal[k]
        multi_trip_indiv_temp <- subset_sex[subset_sex$ID==ID_animal_temp,]
        # Calculate kernel
        k_multi <- kernelUD(multi_trip_indiv_temp[,3], grid=xy, extent=1, h=smooth_h)
        # Mean of values per cell
        k_multi_temp <- stack(lapply(k_multi, raster))
        udm <- mean(k_multi_temp)
        udm@data@names <- ID_animal_temp
        # Add raster to list
        multi_trips_list[[k]] <- udm
        }
    }else{}
    
    #------------------- CALCULATE FOR INDIV WITH SINGLE TRIPS ----------------
    
    # Add sex to list of single trips
    for(m in 1:nrow(single_trips)){
      ID_temp <- single_trips$ID_animal[m]
      single_trips$sex[m] <- spdata_rdm_sex$sex[which(spdata_rdm_sex$ID==ID_temp)] 
    }
    # Subset single trip birds per sex
    single_trips_sex <- single_trips[single_trips$sex==sex,]
    # Create empty list to store rasters
    single_trips_list <- vector(mode = "list", length = nrow(single_trips_sex))
    # Run KDE per bird
    for(j in 1:nrow(single_trips_sex)){
      ID_animal_temp <- single_trips_sex$ID_animal[j]
      single_trip_temp <- subset_sex[subset_sex$ID==ID_animal_temp,]
      # Calculate kernel
      k_single <- kernelUD(single_trip_temp[,3], grid=xy, extent=1, h=smooth_h)
      # Mean of values per cell
      k_single_temp <- stack(lapply(k_single, raster))
      udm <- mean(k_single_temp)
      udm@data@names <- ID_animal_temp
      # Add raster to list
      single_trips_list[[j]] <- udm
    }
    
    #------------- MERGE MUTLI & SINGLE LISTS + CONVERT INTO estUDm -------------
    
    # Combine list of rasters
    list_combined <- c(multi_trips_list, single_trips_list)
    # Stack the rasters
    combined_temp <- stack(list_combined)
    # Mean of values per cell
    udm_all <- mean(combined_temp)
    # Replace values 
    k_mf_rdm[[i]]@data$ud <- udm_all@data@values[k_mf_rdm[[i]]@grid.index]
  }
  
  #------------------------ Quantify + store overlap --------------------------
  
  # Quantify overlap (Battacharyya's affinity)
  overlap_95_mf_rdm <- kerneloverlaphr(k_mf_rdm, meth="BA", conditional=TRUE, 
                                       percent=95)
  overlap_50_mf_rdm <- kerneloverlaphr(k_mf_rdm, meth="BA", conditional=TRUE, 
                                       percent=50)
  # Store overlap values
  overlap_temp <- data.frame(l,overlap_95_mf_rdm[1,2],overlap_50_mf_rdm[1,2])
  colnames(overlap_temp) <- c("iteration","95%","50%")
  overlap_matrix <- rbind(overlap_matrix,overlap_temp)
  # Print running time of the loop
  if(l==iterations){
    timer <- toc()
    elapsed <- (timer$toc-timer$tic)/60
    print(paste0(round(elapsed,3)," min elapsed"))
  }else{}
}

# Save overlap matrix
write.csv(overlap_matrix, file=paste0("mf_overlap_random.csv"), 
          row.names = FALSE, quote=FALSE)


#---------------------- Determine significance of overlap ---------------------


overlap_matrix <- read.csv(file = "HR/mf_overlap_random.csv", header = TRUE, 
                           sep = ",", na = "NA")

summary(overlap_matrix)
sd(overlap_matrix$X95.)
sd(overlap_matrix$X50.)

# p-values were determined as the proportion of randomized overlaps that were 
# smaller than the observed (Clay et al. 2017)

p_50 <- nrow(overlap_matrix[overlap_matrix$X50.< 
                              overlap_50_mf[2,1],])/nrow(overlap_matrix)
p_95 <- nrow(overlap_matrix[overlap_matrix$X95.< 
                              overlap_95_mf[2,1],])/nrow(overlap_matrix)

# Histograms
par(mfrow=c(1,2))
hist(overlap_matrix$X50., breaks=20, main="50%", xlab="overlap (BA)", 
     xlim = c(0,0.5))
text(x = 0.1, y = 100, paste0("p=",round(p_50,2)))
abline(v = overlap_50_mf[2,1], col="red")
hist(overlap_matrix$X95., breaks=10, main="95%", xlab="overlap (BA)", 
     xlim = c(0,1))
text(x = 0.2, y = 160, paste0("p=",round(p_95,2)))
abline(v = overlap_95_mf[2,1], col="red")


###############################################################################
#                                 KDE PER BR
###############################################################################


smooth_h <- 186*1000
BR_list <- sort(unique(crds_temp$BR))

# Compute UD per year
k_BR <- kernelUD(spdata[,4], grid=xy, extent=1, h=smooth_h)
image(k_BR)

for(i in 1:length(BR_list)){   
  BR <- BR_list[i]
  # Subset data per year
  subset_BR <- spdata[spdata$BR==BR,]
  # Create empty list to store rasters
  BR_bird_list <- vector(mode = "list", length = 
                           length(unique(subset_BR$ID_animal)))
  # Run KDE per bird
  for(j in 1:length(unique(subset_BR$ID_animal))){
    ID_animal_temp <- unique(subset_BR$ID_animal)[j]
    # Subset data per bird ID
    subset_ID_animal <- subset_BR[subset_BR$ID_animal==ID_animal_temp,]
    # Create empty list to store rasters
    BR_trips_list <- vector(mode = "list", length = 
                              length(unique(subset_ID_animal$ID)))
    # Merge track UDs per bird (if multiple trips)
    for(h in 1:length(unique(subset_ID_animal$ID))){
      ID_temp <- unique(subset_ID_animal$ID)[h]
      # Subset data per track ID
      ID_trip_spdata <- subset_ID_animal[subset_ID_animal$ID==ID_temp,]
      # Calculate kernel
      k_sub_BR <- kernelUD(ID_trip_spdata[,3], grid=xy, extent=1, h=smooth_h)
      # Mean of values per cell
      k_temp <- stack(lapply(k_sub_BR, raster))
      udm <- mean(k_temp)
      udm@data@names <- paste0(ID_temp,"/",ID_animal_temp)
      # Add raster to list
      BR_trips_list[[h]] <- udm
    }
    # Stack the rasters
    combined_birds <- stack(BR_trips_list)
    # Mean of values per cell
    udm_birds <- mean(combined_birds)
    # Add raster to list
    BR_bird_list[[j]] <- udm_birds
  }
  # Stack the rasters
  combined_temp <- stack(BR_bird_list)
  # Mean of values per cell
  udm_all <- mean(combined_temp)
  # Replace values 
  k_BR[[i]]@data$ud <- udm_all@data@values[k_BR[[i]]@grid.index]
  # Remove temporary files from environment
  if(i==length(BR_list)){
    rm(i,j,BR_trips_list,combined_temp,udm_all,k_temp,k_sub_BR,ID_trip_spdata,
       ID_temp,subset_BR,udm,BR,ID_animal_temp,ID_track,udm_birds, 
       subset_ID_animal)
  }
}

image(k_BR)

# Get isopleths
pol99_BR <- getverticeshr(k_BR,99)
pol50_BR <- getverticeshr(k_BR,50)
mapview(pol99_BR)
mapview(pol50_BR)

# Calculate UDs overlap between BR
BA_95 <- kerneloverlaphr(k_BR, meth="BA", conditional=TRUE, percent=95)
BA_50 <- kerneloverlaphr(k_BR, meth="BA", conditional=TRUE, percent=50)

###############################################################################
#                      KDE PER BR RANDOMISED (LOOP)
###############################################################################


iterations <- 1000
smooth_h <- 186*1000
BR_list <- sort(unique(crds_temp$BR))

# Load data

BR_index <- GLS_summary_data[,c("ID", "ID_animal", "BR")]
BR_index <- BR_index[BR_index$ID %in% spdata$ID,]
BR_vector <- BR_index$BR


overlap_matrix <- NULL


#------------------------------- !!!! LOOP !!!! -------------------------------

for(l in 1:iterations){
  if(l==1){tic()}else{}
  # Print iterations to follow progress
  print(paste0(l,"/",iterations))
  #------------------------------ Randomise BR ------------------------------
  # Create a vector of random BRs
  rdm_BR <- sample(BR_vector,size=length(BR_vector),replace=FALSE)
  # Add random vector to BR index
  BR_index$rdm <- rdm_BR
  #--------------------------------- KDE PER BR -----------------------------
  # Create an sp object with randomised BRs
  spdata_rdm_BR <- spdata
  for(i in 1:nrow(spdata_rdm_BR)){
    ID_temp <- spdata_rdm_BR$ID[i]
    spdata_rdm_BR$rdm_BR[i] <- BR_index$rdm[BR_index$ID==ID_temp]  
  }
  # Compute UD per BR
  k_BR_rdm <- kernelUD(spdata[,4], grid=xy, extent=1, h=smooth_h)
  for(i in 1:length(BR_list)){   
    BR <- BR_list[i]
    # Subset data per BR
    subset_BR <- spdata_rdm_BR[spdata_rdm_BR$rdm_BR==BR,]
    # Create empty list to store rasters
    BR_bird_list <- vector(mode = "list", length = 
                             length(unique(subset_BR$ID_animal)))
    # Run KDE per bird
    for(j in 1:length(unique(subset_BR$ID_animal))){
      ID_animal_temp <- unique(subset_BR$ID_animal)[j]
      # Subset data per bird ID
      subset_ID_animal <- subset_BR[subset_BR$ID_animal==ID_animal_temp,]
      # Create empty list to store rasters
      BR_trips_list <- vector(mode = "list", length = 
                                length(unique(subset_ID_animal$ID)))
      # Merge track UDs per bird (if multiple trips)
      for(h in 1:length(unique(subset_ID_animal$ID))){
        ID_temp <- unique(subset_ID_animal$ID)[h]
        # Subset data per track ID
        ID_trip_spdata <- subset_ID_animal[subset_ID_animal$ID==ID_temp,]
        # Calculate kernel
        k_sub_BR <- kernelUD(ID_trip_spdata[,3], grid=xy, extent=1, h=smooth_h)
        # Mean of values per cell
        k_temp <- stack(lapply(k_sub_BR, raster))
        udm <- mean(k_temp)
        udm@data@names <- paste0(ID_temp,"/",ID_animal_temp)
        # Add raster to list
        BR_trips_list[[h]] <- udm
      }
      # Stack the rasters
      combined_birds <- stack(BR_trips_list)
      # Mean of values per cell
      udm_birds <- mean(combined_birds)
      # Add raster to list
      BR_bird_list[[j]] <- udm_birds
    }
    # Stack the rasters
    combined_temp <- stack(BR_bird_list)
    # Mean of values per cell
    udm_all <- mean(combined_temp)
    # Replace values 
    k_BR_rdm[[i]]@data$ud <- udm_all@data@values[k_BR_rdm[[i]]@grid.index]
    # Remove temporary files from environment
    if(i==length(BR_list)){
       rm(i,j,BR_trips_list,combined_temp,udm_all,k_temp,k_sub_BR,ID_trip_spdata,
          ID_temp,subset_BR,udm,BR,ID_animal_temp,ID_track,udm_birds, 
          subset_ID_animal)
    }
  }
  #------------------------ Quantify + store overlap --------------------------
  # Quantify overlap (Battacharyya's affinity)
  BA_95_rdm <- kerneloverlaphr(k_BR_rdm, meth="BA", conditional=TRUE, percent=95)
  BA_50_rdm <- kerneloverlaphr(k_BR_rdm, meth="BA", conditional=TRUE, percent=50)
  # Store overlap values
  overlap_temp <- as.data.frame(t(as.matrix(c(l,BA_95_rdm[1,2],BA_50_rdm[1,2]))))
  overlap_matrix <- rbind(overlap_matrix,overlap_temp)
  # Print running time of the loop + final formatting of dataframes
  if(l==iterations){
    timer <- toc()
    elapsed <- (timer$toc-timer$tic)/60
    print(paste0(round(elapsed,3)," min elapsed"))
    colnames(overlap_matrix) <- c("iteration","95%","50%")
    rm(i,k,j,l, BA_50_rdm, BA_95_rdm, BR_trips_list, combined_temp, udm_all, k_temp, 
       k_sub_BR, ID_trip_spdata, ID_temp,subset_BR,udm,BR, rdm_BR, overlap_50_rdm,
       overlap_95_rdm)
  }
}

# Save overlap matrices
write.csv(overlap_matrix, file=paste0("overlap_rdm_BR.csv"), 
          row.names = FALSE, quote=FALSE)


#---------------------- Determine significance of overlap ---------------------


overlap_matrix <-  read.csv(file = "overlap_rdm_BR.csv", header = TRUE, 
                            sep = ",", na = "NA")

summary(overlap_matrix)
BA_50
BA_95

# p values
p_50 <- nrow(overlap_matrix[overlap_matrix$X50.< BA_50[2,1],])/nrow(overlap_matrix)
p_95 <- nrow(overlap_matrix[overlap_matrix$X95.< BA_95[2,1],])/nrow(overlap_matrix)

# Histograms
par(mfrow=c(1,2))
hist(overlap_matrix$X50., breaks=20, main="50%", xlab="overlap (BA)", 
     xlim = c(0,0.5), ylim=c(0,200))
text(x = 0.1, y = 80, paste0("p=",round(p_50,2)))
abline(v = BA_50[2,1], col="red", lwd=2)
hist(overlap_matrix$X95., breaks=10, main="95%", xlab="overlap (BA)", 
     xlim = c(0,1))
text(x = 0.2, y = 160, paste0("p=",round(p_95,2)))
abline(v = BA_95[2,1], col="red")



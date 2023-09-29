###############################################################################
#                          ADD LOC ESTIMATES SD FROM SGAT
#                           AND CALCULATE MEAN ACCURACY
###############################################################################

rm(list=ls())

library(fossil)
library(dplyr)

# Load data
crds_temp <- read.csv(file = "crds_all_pts_land_adj.4.csv", header = TRUE, 
                      sep = ",", na = "NA")
GLS_index <- read.csv("Abbotts_GLS_Index.csv", header = TRUE)

#---------------- RETRIEVE SD OF LOC ESTIMATES FROM MCMC (SGAT) ---------------

for(i in 1:length(unique(crds_temp$ID))){
  ID_temp <- unique(crds_temp$ID)[i]
  ID_GLS <- gsub("_1","",ID_temp)
  ID_GLS <- gsub("_2","",ID_GLS)
  ID_subset_vector <- crds_temp$timestamp[crds_temp$ID==ID_temp]
  estelle_temp <- read.csv(file = paste0(
    "GLS_files/",ID_GLS,"/",ID_temp,"_2estelleSST.csv"), header = TRUE, 
    sep = ",", na = "NA")
  estelle_temp$Time <- as.POSIXct(estelle_temp$Time)
  lon_sd_vector <- NULL
  lat_sd_vector <- NULL
  for(j in 1:length(ID_subset_vector)){
    lon_sd_vector <- c(lon_sd_vector, 
                       estelle_temp$Lon.sd[estelle_temp$Time==ID_subset_vector[j]])
    lat_sd_vector <- c(lat_sd_vector, 
                       estelle_temp$Lat.sd[estelle_temp$Time==ID_subset_vector[j]])
  }
  crds_temp$lat_sd[crds_temp$ID==ID_temp] <- lat_sd_vector
  crds_temp$lon_sd[crds_temp$ID==ID_temp] <- lon_sd_vector
}
rm(i,j, ID_GLS, ID_temp, ID_subset_vector, estelle_temp, lon_sd_vector, 
   lat_sd_vector)

# Check data
summary(crds_temp$lat_sd)
summary(crds_temp$lon_sd)

# Save file
write.csv(crds_temp, file="crds_all_pts_land_adj.4.csv", 
          row.names = FALSE, quote=FALSE)

#--------------------------- CALCULATE SD IN KM -------------------------------

crds_temp <- read.csv(file = "crds_all_pts_land_adj.4.csv", header = TRUE, 
                      sep = ",", na = "NA")


crds_temp$lon_sd_km <- deg.dist(crds_temp$Lon, crds_temp$Lat,
                                crds_temp$Lon-crds_temp$lon_sd, crds_temp$Lat)

crds_temp$lat_sd_km <- deg.dist(crds_temp$Lon, crds_temp$Lat,
                                crds_temp$Lon, crds_temp$Lat-crds_temp$lat_sd)

# Add GLS ID for grouping
crds_temp$GLS_ID <- gsub("_1","",crds_temp$ID)
crds_temp$G

# FOR LATITUDE

# Mean latitude OUTSIDE equinox period (+/- 3 weeks)
mean(crds_temp$lat_sd_km[crds_temp$eq_rmv==0])
sd(crds_temp$lat_sd_km[crds_temp$eq_rmv==0])
lat_sd_track_eq_0 <- crds_temp[crds_temp$eq_rmv==0,] %>% group_by(ID_GLS) %>% 
  summarise(Mean=mean(lat_sd_km), SD=sd(lat_sd_km))
mean(lat_sd_track_eq_0$Mean)
mean(lat_sd_track_eq_0$SD)

# Mean latitude AROUND equinox period (+/- 3 weeks)
mean(crds_temp$lat_sd_km[crds_temp$eq_rmv==1])
sd(crds_temp$lat_sd_km[crds_temp$eq_rmv==1])
lat_sd_track_eq_1 <- crds_temp[crds_temp$eq_rmv==1,] %>% group_by(ID_GLS) %>% 
  summarise(Mean=mean(lat_sd_km), SD=sd(lat_sd_km))
mean(lat_sd_track_eq_1$Mean)
mean(lat_sd_track_eq_1$SD)

# FOR LONGITUDE

# Mean longitude OUTSIDE equinox period (+/- 3 weeks)
mean(crds_temp$lon_sd_km[crds_temp$eq_rmv==0])
sd(crds_temp$lon_sd_km[crds_temp$eq_rmv==0])
lon_sd_track_eq_0 <- crds_temp[crds_temp$eq_rmv==0,] %>% group_by(ID_GLS) %>% 
  summarise(Mean=mean(lon_sd_km), SD=sd(lon_sd_km))
mean(lon_sd_track_eq_0$Mean)
mean(lon_sd_track_eq_0$SD)

# Mean longitude AROUND equinox period (+/- 3 weeks)
mean(crds_temp$lon_sd_km[crds_temp$eq_rmv==1])
sd(crds_temp$lon_sd_km[crds_temp$eq_rmv==1])
lon_sd_track_eq_1 <- crds_temp[crds_temp$eq_rmv==1,] %>% group_by(ID_GLS) %>% 
  summarise(Mean=mean(lon_sd_km), SD=sd(lon_sd_km))
mean(lon_sd_track_eq_1$Mean)
mean(lon_sd_track_eq_1$SD)


# --------------------------------PER MODEL

# LATITUDE

# Around equinox
# All points
mean(crds_temp$lat_sd_km[crds_temp$eq_rmv==1 & crds_temp$ID_GLS %in% 
                           GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]])
sd(crds_temp$lat_sd_km[crds_temp$eq_rmv==1 & crds_temp$ID_GLS %in% 
                         GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]])
mean(crds_temp$lat_sd_km[crds_temp$eq_rmv==1 & 
                           !(crds_temp$ID_GLS %in% 
                               GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"])])
sd(crds_temp$lat_sd_km[crds_temp$eq_rmv==1 & 
                         !(crds_temp$ID_GLS %in% 
                             GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"])])

# Averaged per GLS first
lat_sd_track_eq_1_MK7 <- 
  crds_temp[crds_temp$eq_rmv==1 & 
              crds_temp$ID_GLS %in% 
              GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"],] %>% 
  group_by(ID_GLS) %>% summarise(Mean=mean(lat_sd_km), SD=sd(lat_sd_km))
mean(lat_sd_track_eq_1_MK7$Mean)
mean(lat_sd_track_eq_1_MK7$SD)

lat_sd_track_eq_1_other <- 
  crds_temp[crds_temp$eq_rmv==1 & 
              !(crds_temp$ID_GLS %in% 
                  GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]),] %>% 
  group_by(ID_GLS) %>% summarise(Mean=mean(lat_sd_km), SD=sd(lat_sd_km))
mean(lat_sd_track_eq_1_other$Mean)
mean(lat_sd_track_eq_1_other$SD)

# Outside equinox
# All points
mean(crds_temp$lat_sd_km[crds_temp$eq_rmv==0 & 
                           crds_temp$ID_GLS %in% 
                           GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]])
sd(crds_temp$lat_sd_km[crds_temp$eq_rmv==0 & 
                         crds_temp$ID_GLS %in% 
                         GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]])
mean(crds_temp$lat_sd_km[crds_temp$eq_rmv==0 & 
                           !(crds_temp$ID_GLS %in% 
                               GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"])])
sd(crds_temp$lat_sd_km[crds_temp$eq_rmv==0 & 
                         !(crds_temp$ID_GLS %in% 
                             GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"])])

# Averaged per GLS first
lat_sd_track_eq_0_MK7 <- 
  crds_temp[crds_temp$eq_rmv==0 & 
              crds_temp$ID_GLS %in% 
              GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"],] %>% 
  group_by(ID_GLS) %>% summarise(Mean=mean(lat_sd_km), SD=sd(lat_sd_km))
mean(lat_sd_track_eq_0_MK7$Mean)
mean(lat_sd_track_eq_0_MK7$SD)

lat_sd_track_eq_0_other <- 
  crds_temp[crds_temp$eq_rmv==0 & 
              !(crds_temp$ID_GLS %in% 
                  GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]),] %>% 
  group_by(ID_GLS) %>% summarise(Mean=mean(lat_sd_km), SD=sd(lat_sd_km))
mean(lat_sd_track_eq_0_other$Mean)
mean(lat_sd_track_eq_0_other$SD)


# LONGITUDE

# Around equinox
# All points
mean(crds_temp$lon_sd_km[crds_temp$eq_rmv==1 & 
                           crds_temp$ID_GLS %in% 
                           GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]])
sd(crds_temp$lon_sd_km[crds_temp$eq_rmv==1 & 
                         crds_temp$ID_GLS %in% 
                         GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]])
mean(crds_temp$lon_sd_km[crds_temp$eq_rmv==1 & 
                           !(crds_temp$ID_GLS %in% 
                               GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"])])
sd(crds_temp$lon_sd_km[crds_temp$eq_rmv==1 & 
                         !(crds_temp$ID_GLS %in% 
                             GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"])])

# Averaged per GLS first
lon_sd_track_eq_1_MK7 <- 
  crds_temp[crds_temp$eq_rmv==1 & 
              crds_temp$ID_GLS %in% 
              GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"],] %>% 
  group_by(ID_GLS) %>% summarise(Mean=mean(lon_sd_km), SD=sd(lon_sd_km))
mean(lon_sd_track_eq_1_MK7$Mean)
mean(lon_sd_track_eq_1_MK7$SD)

lon_sd_track_eq_1_other <- 
  crds_temp[crds_temp$eq_rmv==1 & 
              !(crds_temp$ID_GLS %in% 
                  GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]),] %>% 
  group_by(ID_GLS) %>% summarise(Mean=mean(lon_sd_km), SD=sd(lon_sd_km))
mean(lon_sd_track_eq_1_other$Mean)
mean(lon_sd_track_eq_1_other$SD)

# Outside equinox
# All points
mean(crds_temp$lon_sd_km[crds_temp$eq_rmv==0 & 
                           crds_temp$ID_GLS %in% 
                           GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]])
sd(crds_temp$lon_sd_km[crds_temp$eq_rmv==0 & 
                         crds_temp$ID_GLS %in% 
                         GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]])
mean(crds_temp$lon_sd_km[crds_temp$eq_rmv==0 & 
                           !(crds_temp$ID_GLS %in% 
                               GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"])])
sd(crds_temp$lon_sd_km[crds_temp$eq_rmv==0 & 
                         !(crds_temp$ID_GLS %in% 
                             GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"])])

# Averaged per GLS first
lon_sd_track_eq_0_MK7 <- 
  crds_temp[crds_temp$eq_rmv==0 & 
              crds_temp$ID_GLS %in% 
              GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"],] %>% 
  group_by(ID_GLS) %>% summarise(Mean=mean(lon_sd_km), SD=sd(lon_sd_km))
mean(lon_sd_track_eq_0_MK7$Mean)
mean(lon_sd_track_eq_0_MK7$SD)

lon_sd_track_eq_0_other <- 
  crds_temp[crds_temp$eq_rmv==0 & 
              !(crds_temp$ID_GLS %in% 
                  GLS_index$GLS_ID[GLS_index$GLS_model=="MK7"]),] %>% 
  group_by(ID_GLS) %>% summarise(Mean=mean(lon_sd_km), SD=sd(lon_sd_km))
mean(lon_sd_track_eq_0_other$Mean)
mean(lon_sd_track_eq_0_other$SD)










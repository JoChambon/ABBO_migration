###############################################################################
#                            HABITAT MODELLING
#                             Data preparation
###############################################################################

rm(list=ls())

library(sp)
library(fossil)
library(raster)
library(rgdal)
library(rgeos)
library(rworldmap)
library(mapview)
library(ncdf4)
library(mgcv)
library(adehabitatHR)
library(lme4)

# Load locations and GLS summary
crds <- read.csv(file = "crds_all_pts_land_adj.4.csv", header = TRUE, 
                 sep = ",", na = "NA") # Points on CI already removed
GLS_summary_data <- read.csv(file="NB_metrics_all_tracks.4.csv", header = TRUE, sep = ",", na = "NA")
crds$timestamp <- as.POSIXct(crds$timestamp)

# Removing equinox data, points on land, and subset time spent at non-breeding area only
crds_sub <- subset(crds, eq_rmv==0 & phenophase=="NB" & pts_on_land==0)


#---------------------------- USE & RANDOM POINTS -----------------------------


# Set use/available ratio
n_rdm <- 50

# Load 99.9% UD
load(file = "RSF/HR99.9.RData")

# Remove land from polygon
world <- spTransform(getMap(resolution = "low"), CRS("+init=epsg:3857"))
UD99_adj <- gDifference(pol99_all, world)
mapview(UD99_adj)

# Keep tracks with at least 30 locations at NB
ID_keep <- as.data.frame(table(crds_sub$ID))
ID_keep <- subset(ID_keep, !(Freq<30))
crds_sub <- subset(crds_sub, ID %in% ID_keep$Var1)

# Create random points, merge with used points, and reformat
rsf_data <- NULL
for(i in 1:length(unique(crds_sub$ID))){
  ID_temp <- unique(crds_sub$ID)[i]
  crds_temp <- subset(crds_sub, ID==ID_temp)
  crds_temp <- crds_temp[,c("x","y","timestamp","Lon","Lat","ID","animal_ID")]
  crds_temp$case <- TRUE
  xy_random <- as.data.frame(spsample(UD99_adj, nrow(crds_temp)*n_rdm, "random")@coords)
  temp_timestamp <- NULL
  for(j in 1:n_rdm){
    temp_timestamp <- rbind(temp_timestamp, as.data.frame(crds_temp$timestamp))
  }
  xy_random$timestamp <- temp_timestamp[,1]
  pj <- as.data.frame(project(as.matrix(xy_random[,c(1,2)]), "+init=epsg:3857", inv=TRUE))
  plot(pj$x, pj$y, asp=1, col = "darkblue", pch = 19, cex = 0.5)
  points(crds_temp$Lon, crds_temp$Lat, pch = 19, col = "orange", cex = 0.5)
  xy_random$Lon <- pj$x
  xy_random$Lat <- pj$y
  xy_random$ID <- ID_temp
  xy_random$animal_ID <- unique(crds_temp$animal_ID)
  xy_random$case <- FALSE
  use_rdm_ID <- rbind(crds_temp, xy_random)
  rsf_data <- rbind(rsf_data, use_rdm_ID)
  if(i == length(unique(crds_sub$ID))){
    rm(xy_random, i, j, pj, use_rdm_ID,ID_temp, crds_temp, temp_timestamp)
  }
}


rm(world, crds_sub, UD99_adj, pol99_all, ID_keep)


#----------------------------------- SST --------------------------------------


#load(file="RSF/rsf_data_nrdm50_locvalue.RData")

# Load netcdf files of environment data
sst_daily <- nc_open("ENV_DATA/SST/OISST_daily_0.25x0.25/OISST_2008_2015.nc")

# Extract variables of interest
sst_vals <- ncvar_get(sst_daily, 'sst')

# EXTRACT SST
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$sst <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat/time in env data from location
  Ln <- which(abs(sst_daily$dim$lon$vals-rsf_data$Lon[j])==min(abs(sst_daily$dim$lon$vals-rsf_data$Lon[j])))
  Lt <- which(abs(sst_daily$dim$lat$vals-rsf_data$Lat[j])==min(abs(sst_daily$dim$lat$vals-rsf_data$Lat[j])))
  Tm <- which((as.POSIXct(strftime(as.POSIXct(sst_daily$dim$time$vals*86400, origin = "1800-01-01"), 
                                   format = "%Y-%m-%d", tz = "UTC"), format = "%Y-%m-%d", tz = "UTC"))  
              == as.POSIXct(strftime(rsf_data$timestamp[j],format = "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC"))
  rsf_data$sst[j] <- sst_vals[Ln, Lt, Tm]
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, Tm, j)
    print("sst")
    print(summary(rsf_data$sst))
    tictoc::toc()
  }
}

# Check distribution of values
par(mfrow=c(2,1))
hist(rsf_data$sst[rsf_data$case==TRUE], xlim = c(min(rsf_data$sst, na.rm = TRUE), max(rsf_data$sst, na.rm = TRUE)))
hist(rsf_data$sst[rsf_data$case==FALSE], xlim = c(min(rsf_data$sst, na.rm = TRUE), max(rsf_data$sst, na.rm = TRUE)),
     breaks=20)
par(mfrow=c(1,1))

# Remove env data files once satisfied of the data extraction
rm(sst_daily, sst_vals)
gc()


#----------------------------------- Chl --------------------------------------


#load(file="RSF/rsf_data_nrdm50_locvalue.RData")

# Load netcdf files of environment data
chl_daily <- nc_open("ENV_DATA/chl/cmems_mod_glo_bgc_my_0.25_P1D-m_CHL.nc")

# Extract variables of interest
chl_vals <- ncvar_get(chl_daily, 'chl')

# EXTRACT chl AND CHL A
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$chl <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat/time in env data from location
  Ln <- which(abs(chl_daily$dim$lon$vals-rsf_data$Lon[j])==min(abs(chl_daily$dim$lon$vals-rsf_data$Lon[j])))
  Lt <- which(abs(chl_daily$dim$lat$vals-rsf_data$Lat[j])==min(abs(chl_daily$dim$lat$vals-rsf_data$Lat[j])))
  Tm <- which((as.POSIXct(strftime(as.POSIXct(chl_daily$dim$time$vals*3600, origin = "1950-01-01"), 
                                   format = "%Y-%m-%d", tz = "UTC"), format = "%Y-%m-%d", tz = "UTC"))  
              == as.POSIXct(strftime(rsf_data$timestamp[j],format = "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC"))
  rsf_data$chl[j] <- chl_vals[Ln, Lt, Tm]
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, Tm, j)
    print("chl a")
    print(summary(rsf_data$chl))
    tictoc::toc()
  }
}

# Check distribution of values
par(mfrow=c(2,1))
hist(rsf_data$chl[rsf_data$case==TRUE], xlim = c(0, 0.3), breaks=50)
hist(rsf_data$chl[rsf_data$case==FALSE], xlim = c(0, 0.3), breaks=300)
par(mfrow=c(1,1))

# Remove env data files once satisfied of the data extraction
rm(chl_daily, chl_vals)
gc()


#-------------------------------- Salinity ------------------------------------


#load(file="RSF/rsf_data_nrdm50_locvalue.RData")

# Load netcdf files of environment data
sal_daily <- nc_open("ENV_DATA/sal/global-reanalysis-phy-001-031-grepv2-mnstd-daily_sal.nc")

# Extract variables of interest
sal_vals <- ncvar_get(sal_daily, 'so_mean')

# EXTRACT sal
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$sal <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat/time in env data from location
  Ln <- which(abs(sal_daily$dim$longitude$vals-rsf_data$Lon[j])==min(abs(sal_daily$dim$longitude$vals-rsf_data$Lon[j])))
  Lt <- which(abs(sal_daily$dim$latitude$vals-rsf_data$Lat[j])==min(abs(sal_daily$dim$latitude$vals-rsf_data$Lat[j])))
  Tm <- which((as.POSIXct(strftime(as.POSIXct(sal_daily$dim$time$vals*86400, origin = "1950-01-01"), 
                                   format = "%Y-%m-%d", tz = "UTC"), format = "%Y-%m-%d", tz = "UTC"))  
              == as.POSIXct(strftime(rsf_data$timestamp[j],format = "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC"))
  rsf_data$sal[j] <- sal_vals[Ln, Lt, Tm]
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, Tm, j)
    print("sal")
    print(summary(rsf_data$sal))
    tictoc::toc()
  }
}

# Check distribution of values
par(mfrow=c(2,1))
hist(rsf_data$sal[rsf_data$case==TRUE], xlim = c(min(rsf_data$sal, na.rm = TRUE), max(rsf_data$sal, na.rm = TRUE)))
hist(rsf_data$sal[rsf_data$case==FALSE], xlim = c(min(rsf_data$sal, na.rm = TRUE), max(rsf_data$sal, na.rm = TRUE)),
     breaks=50)
par(mfrow=c(1,1))

# Remove env data files once satisfied of the data extraction
rm(sal_daily, sal_vals)
gc()


#----------------------------- Depth + slope ----------------------------------


#load(file="RSF/rsf_data_nrdm50_locvalue.RData")

# Load netcdf files of environment data
bat_data <- nc_open("ENV_DATA/bat/GEBCO/gebco_2021_n25.0_s-36.0_w70.0_e170.0.nc")

# Extract variables of interest
bat_vals <- ncvar_get(bat_data, 'elevation')
lon <- ncvar_get(bat_data, "lon")
lat <- ncvar_get(bat_data, "lat", verbose = F)

# Annotate locations
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$bat <- NA
    rsf_data$slope <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat in env data from location
  Ln <- which(abs(bat_data$dim$lon$vals-rsf_data$Lon[j])== min(abs(bat_data$dim$lon$vals-rsf_data$Lon[j])))
  Lt <- which(abs(bat_data$dim$lat$vals-rsf_data$Lat[j])== min(abs(bat_data$dim$lat$vals-rsf_data$Lat[j])))
  rsf_data$bat[j] <- bat_vals[Ln, Lt]
  # Subset data around jth location to facilitate computation
  ln_array <- lon[seq((Ln-1200),(Ln+1200), by=24)] # 1200 for 5° as data at 15 arc second res, seq by 24 to resample at 0.1°
  lt_array <- lat[seq((Lt-1200),(Lt+1200), by=24)]
  bat_temp <- stack(as.data.frame(bat_vals[seq((Ln-1200),(Ln+1200), by=24), seq((Lt-1200),(Lt+1200), by=24)]))[1]
  bat_temp$lon <- matrix(ln_array,nrow = length(ln_array)*length(ln_array),ncol = 1)
  bat_temp$lat <- rep(lt_array, each=length(lt_array))
  bat_temp <- bat_temp[bat_temp$values<0,]
  # Calculate distance between grid cells centre and jth location
  for(i in 1:nrow(bat_temp)){
    bat_temp$dist[i] <- deg.dist(bat_temp$lon[i], bat_temp$lon[i], rsf_data$Lon[j], rsf_data$Lon[j])
  }
  # SD of bat in 385km radius from location
  rsf_data$slope[j] <- sd(bat_temp$values[bat_temp$dist<=385], na.rm=TRUE)
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, j, ln_array, lt_array, bat_temp, i)
    print("Bathymetry")
    print(summary(rsf_data$bat))
    print("Slope")
    print(summary(rsf_data$slope))
    tictoc::toc()
  }
}

# Check distribution of values
par(mfrow=c(2,1))
hist(rsf_data$bat[rsf_data$case==TRUE], xlim = c(min(rsf_data$bat, na.rm = TRUE), max(rsf_data$bat, na.rm = TRUE)))
hist(rsf_data$bat[rsf_data$case==FALSE], xlim = c(min(rsf_data$bat, na.rm = TRUE), max(rsf_data$bat, na.rm = TRUE)),
     breaks=20)
hist(rsf_data$slope[rsf_data$case==TRUE], xlim = c(min(rsf_data$slope, na.rm = TRUE), max(rsf_data$slope, na.rm = TRUE)))
hist(rsf_data$slope[rsf_data$case==FALSE], xlim = c(min(rsf_data$slope, na.rm = TRUE), max(rsf_data$slope, na.rm = TRUE)),
     breaks=15)
par(mfrow=c(1,1))


# Remove env data files once satisfied with the data extraction
rm(bat_data, bat_vals, lon, lat)
gc()


#----------------------------- SLOPE (diff radius) ----------------------------


#load(file="RSF/rsf_data_nrdm50_locvalue.RData")

# Load netcdf files of environment data
bat_data <- nc_open("ENV_DATA/bat/GEBCO/gebco_2021_n25.0_s-36.0_w70.0_e170.0.nc")

# Extract variables of interest
bat_vals <- ncvar_get(bat_data, 'elevation')
lon <- ncvar_get(bat_data, "lon")
lat <- ncvar_get(bat_data, "lat", verbose = F)

# Annotate locations
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$slope_200 <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat in env data from location
  Ln <- which(abs(bat_data$dim$lon$vals-rsf_data$Lon[j])== min(abs(bat_data$dim$lon$vals-rsf_data$Lon[j])))
  Lt <- which(abs(bat_data$dim$lat$vals-rsf_data$Lat[j])== min(abs(bat_data$dim$lat$vals-rsf_data$Lat[j])))
  # Subset data around jth location to facilitate computation
  ln_array <- lon[seq((Ln-1200),(Ln+1200), by=24)] # 1200 for 5° as data at 15 arc second res, seq by 24 to resample at 0.1°
  lt_array <- lat[seq((Lt-1200),(Lt+1200), by=24)]
  bat_temp <- stack(as.data.frame(bat_vals[seq((Ln-1200),(Ln+1200), by=24), seq((Lt-1200),(Lt+1200), by=24)]))[1]
  bat_temp$lon <- matrix(ln_array,nrow = length(ln_array)*length(ln_array),ncol = 1)
  bat_temp$lat <- rep(lt_array, each=length(lt_array))
  bat_temp <- bat_temp[bat_temp$values<0,]
  # Calculate distance between grid cells centre and jth location
  for(i in 1:nrow(bat_temp)){
    bat_temp$dist[i] <- deg.dist(bat_temp$lon[i], bat_temp$lon[i], rsf_data$Lon[j], rsf_data$Lon[j])
  }
  # SD of bat in 385km radius from location
  rsf_data$slope_200[j] <- sd(bat_temp$values[bat_temp$dist<=200], na.rm=TRUE)
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, j, ln_array, lt_array, bat_temp, i)
    print("Slope")
    print(summary(rsf_data$slope_200))
    tictoc::toc()
  }
}

# Check distribution of values
par(mfrow=c(2,1))
hist(rsf_data$slope_200[rsf_data$case==TRUE], xlim = c(min(rsf_data$slope_200, na.rm = TRUE), max(rsf_data$slope_200, na.rm = TRUE)))
hist(rsf_data$slope_200[rsf_data$case==FALSE], xlim = c(min(rsf_data$slope_200, na.rm = TRUE), max(rsf_data$slope_200, na.rm = TRUE)),
     breaks=15)
par(mfrow=c(1,1))


# Remove env data files once satisfied with the data extraction
rm(bat_data, bat_vals, lon, lat)
gc()


#-------------------------------- SST GRADIENT --------------------------------


#load(file="RSF/rsf_data_nrdm50_locvalue.RData")

# Load netcdf files of environment data
sst_daily <- nc_open("ENV_DATA/SST/OISST_daily_0.25x0.25/OISST_2008_2015.nc")

# Extract variables of interest
sst_vals <- ncvar_get(sst_daily, 'sst')
lon <- ncvar_get(sst_daily, "lon")
lat <- ncvar_get(sst_daily, "lat", verbose = F)
gc()

# EXTRACT SST
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$sst_grad <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat/time in env data from location
  Ln <- which(abs(sst_daily$dim$lon$vals-rsf_data$Lon[j])==min(abs(sst_daily$dim$lon$vals-rsf_data$Lon[j])))
  Lt <- which(abs(sst_daily$dim$lat$vals-rsf_data$Lat[j])==min(abs(sst_daily$dim$lat$vals-rsf_data$Lat[j])))
  Tm <- which((as.POSIXct(strftime(as.POSIXct(sst_daily$dim$time$vals*86400, origin = "1800-01-01"), 
                                   format = "%Y-%m-%d", tz = "UTC"), format = "%Y-%m-%d", tz = "UTC"))  
              == as.POSIXct(strftime(rsf_data$timestamp[j],format = "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC"))
  # Subset data around jth location to facilitate computation
  ln_array <- lon[seq((Ln-20),(Ln+20))] # 20 for 5° as data at 0.25°
  lt_array <- lat[seq((Lt-20),(Lt+20))]
  sst_temp <- stack(as.data.frame(sst_vals[seq((Ln-20),(Ln+20)), seq((Lt-20),(Lt+20)),Tm]))[1]
  sst_temp$lon <- matrix(ln_array,nrow = length(ln_array)*length(ln_array),ncol = 1)
  sst_temp$lat <- rep(lt_array, each=length(lt_array))
  sst_temp <- sst_temp[!is.na(sst_temp$values),]
  # Calculate distance between grid cells centre and jth location
  for(i in 1:nrow(sst_temp)){
    sst_temp$dist[i] <- deg.dist(sst_temp$lon[i], sst_temp$lon[i], rsf_data$Lon[j], rsf_data$Lon[j])
  }
  # SD of bat in 385km radius from location
  rsf_data$sst_grad[j] <- sd(sst_temp$values[sst_temp$dist<=200], na.rm=TRUE)
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, Tm, j)
    print("sst")
    print(summary(rsf_data$sst))
    tictoc::toc()
  }
}

# Check distribution of values
par(mfrow=c(2,1))
hist(rsf_data$sst_grad[rsf_data$case==TRUE], xlim = c(min(rsf_data$sst_grad, na.rm = TRUE), max(rsf_data$sst_grad, na.rm = TRUE)))
hist(rsf_data$sst_grad[rsf_data$case==FALSE], xlim = c(min(rsf_data$sst_grad, na.rm = TRUE), max(rsf_data$sst_grad, na.rm = TRUE)),
     breaks=20)
par(mfrow=c(1,1))

# Remove env data files once satisfied of the data extraction
rm(sst_daily, sst_vals, lon, lat)
gc()


#-------------------------- MIXED LAYER THICKNESS -----------------------------


# Load netcdf files of environment data
mlt_daily <- nc_open("ENV_DATA/MLT/global-reanalysis-phy-001-031-grepv2-mnstd-daily_mlt.nc")

# Extract variables of interest
mlt_vals <- ncvar_get(mlt_daily, 'mlotst_mean')
gc()

# EXTRACT MLT
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$mlt <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat/time in env data from location
  Ln <- which(abs(mlt_daily$dim$longitude$vals-rsf_data$Lon[j])==min(abs(mlt_daily$dim$longitude$vals-rsf_data$Lon[j])))
  Lt <- which(abs(mlt_daily$dim$latitude$vals-rsf_data$Lat[j])==min(abs(mlt_daily$dim$latitude$vals-rsf_data$Lat[j])))
  Tm <- which((as.POSIXct(strftime(as.POSIXct(mlt_daily$dim$time$vals*86400, origin = "1950-01-01"), 
                                   format = "%Y-%m-%d", tz = "UTC"), format = "%Y-%m-%d", tz = "UTC"))  
              == as.POSIXct(strftime(rsf_data$timestamp[j],format = "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC"))
  rsf_data$mlt[j] <- mlt_vals[Ln, Lt, Tm]
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, Tm, j)
    print("mlt")
    print(summary(rsf_data$mlt))
    tictoc::toc()
  }
}

# Check distribution of values
par(mfrow=c(2,1))
hist(rsf_data$mlt[rsf_data$case==TRUE], xlim = c(min(rsf_data$mlt, na.rm = TRUE), max(rsf_data$mlt, na.rm = TRUE)))
hist(rsf_data$mlt[rsf_data$case==FALSE], xlim = c(min(rsf_data$mlt, na.rm = TRUE), max(rsf_data$mlt, na.rm = TRUE)),
     breaks=50)
par(mfrow=c(1,1))

# Remove env data files once satisfied of the data extraction
rm(mlt_daily, mlt_vals)
gc()


#---------------- GEOSTROPHIC VELOCITY ANOMALIES (zonal - u) ------------------


# Load netcdf files of environment data
gva_u_daily <- nc_open("ENV_DATA/EKE/GVA_u/GVA_u.nc")

# Extract variables of interest
gva_u_vals <- ncvar_get(gva_u_daily, 'ugosa')
gc()

# EXTRACT GVA_u
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$gva_u <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat/time in env data from location
  Ln <- which(abs(gva_u_daily$dim$longitude$vals-rsf_data$Lon[j])==min(abs(gva_u_daily$dim$longitude$vals-rsf_data$Lon[j])))
  Lt <- which(abs(gva_u_daily$dim$latitude$vals-rsf_data$Lat[j])==min(abs(gva_u_daily$dim$latitude$vals-rsf_data$Lat[j])))
  Tm <- which((as.POSIXct(strftime(as.POSIXct(gva_u_daily$dim$time$vals*86400, origin = "1950-01-01"), 
                                   format = "%Y-%m-%d", tz = "UTC"), format = "%Y-%m-%d", tz = "UTC"))  
              == as.POSIXct(strftime(rsf_data$timestamp[j],format = "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC"))
  rsf_data$gva_u[j] <- gva_u_vals[Ln, Lt, Tm]
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, Tm, j)
    print("gva_u")
    print(summary(rsf_data$gva_u))
    tictoc::toc()
  }
}

# Check distribution of values
par(mfrow=c(2,1))
hist(rsf_data$gva_u[rsf_data$case==TRUE], xlim = c(min(rsf_data$gva_u, na.rm = TRUE), max(rsf_data$gva_u, na.rm = TRUE)))
hist(rsf_data$gva_u[rsf_data$case==FALSE], xlim = c(min(rsf_data$gva_u, na.rm = TRUE), max(rsf_data$gva_u, na.rm = TRUE)),
     breaks=50)
par(mfrow=c(1,1))

# Remove env data files once satisfied of the data extraction
rm(gva_u_daily, gva_u_vals)
gc()


#------------ GEOSTROPHIC VELOCITY ANOMALIES (meridional - u) -----------------


# Load netcdf files of environment data
gva_v_daily <- nc_open("ENV_DATA/EKE/gva_v/gva_v.nc")

# Extract variables of interest
gva_v_vals <- ncvar_get(gva_v_daily, 'vgosa')
gc()

# EXTRACT gva_v
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$gva_v <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat/time in env data from location
  Ln <- which(abs(gva_v_daily$dim$longitude$vals-rsf_data$Lon[j])==min(abs(gva_v_daily$dim$longitude$vals-rsf_data$Lon[j])))
  Lt <- which(abs(gva_v_daily$dim$latitude$vals-rsf_data$Lat[j])==min(abs(gva_v_daily$dim$latitude$vals-rsf_data$Lat[j])))
  Tm <- which((as.POSIXct(strftime(as.POSIXct(gva_v_daily$dim$time$vals*86400, origin = "1950-01-01"), 
                                   format = "%Y-%m-%d", tz = "UTC"), format = "%Y-%m-%d", tz = "UTC"))  
              == as.POSIXct(strftime(rsf_data$timestamp[j],format = "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC"))
  rsf_data$gva_v[j] <- gva_v_vals[Ln, Lt, Tm]
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, Tm, j)
    print("gva_v")
    print(summary(rsf_data$gva_v))
    tictoc::toc()
  }
}

# Check distribution of values
par(mfrow=c(2,1))
hist(rsf_data$gva_v[rsf_data$case==TRUE], xlim = c(min(rsf_data$gva_v, na.rm = TRUE), max(rsf_data$gva_v, na.rm = TRUE)))
hist(rsf_data$gva_v[rsf_data$case==FALSE], xlim = c(min(rsf_data$gva_v, na.rm = TRUE), max(rsf_data$gva_v, na.rm = TRUE)),
     breaks=50)
par(mfrow=c(1,1))

# Remove env data files once satisfied of the data extraction
rm(gva_v_daily, gva_v_vals)
gc()


#----------------------------  Calculation of EKE -----------------------------


rsf_data$eke <- 1/2*((rsf_data$gva_u)^2+(rsf_data$gva_v)^2)


# Check distribution of values
par(mfrow=c(2,1))
#hist(rsf_data$eke[rsf_data$case==TRUE], xlim = c(min(rsf_data$eke, na.rm = TRUE), max(rsf_data$eke, na.rm = TRUE)))
hist(log(rsf_data$eke[rsf_data$case==TRUE]), xlim=c(-15,0), breaks=200)
hist(log(rsf_data$eke[rsf_data$case==FALSE]), xlim=c(-15,0), breaks=1000)
par(mfrow=c(1,1))


#------------------------- BARRIER LAYED THICKNESS -----------------------------


# Load netcdf files of environment data
BLT_5days <- nc_open("ENV_DATA/BLT/hawaii_soest_40a7_cb5c_7a11_a207_21d2_8665.nc")

# Extract variables of interest
mlt_vals <- ncvar_get(BLT_5days, 'mlt')
mlp_vals <- ncvar_get(BLT_5days, 'mlp')
gc()
date_vec <- as.POSIXct(strftime(as.POSIXct(BLT_5days$dim$time$vals, origin = "1970-01-01"), 
                                format = "%Y-%m-%d", tz = "UTC"), format = "%Y-%m-%d", tz = "UTC")

# EXTRACT MLT & MLP
for(j in 1:nrow(rsf_data)){
  if(j == 1){
    tictoc::tic("Running time")
    rsf_data$blt_mlt <- NA
    rsf_data$blt_mlp <- NA
  } else {}
  print(paste0("Iteration ", j, "/",nrow(rsf_data)))
  # Find closest lon/lat/time in env data from location
  Ln <- which(abs(BLT_5days$dim$longitude$vals-rsf_data$Lon[j])==min(abs(BLT_5days$dim$longitude$vals-rsf_data$Lon[j])))
  Lt <- which(abs(BLT_5days$dim$latitude$vals-rsf_data$Lat[j])==min(abs(BLT_5days$dim$latitude$vals-rsf_data$Lat[j])))
  Tm <- which(abs(date_vec-as.POSIXct(strftime(rsf_data$timestamp[j],format = "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC")) 
              == min(abs(date_vec - as.POSIXct(strftime(rsf_data$timestamp[j],format = "%Y-%m-%d"), format = "%Y-%m-%d", tz = "UTC"))))
  rsf_data$blt_mlt[j] <- mlt_vals[Ln, Lt, Tm]
  rsf_data$blt_mlp[j] <- mlp_vals[Ln, Lt, Tm]
  if(j == nrow(rsf_data)){
    rm(Lt, Ln, Tm, j)
    print("mlt")
    print(summary(rsf_data$blt_mlt))
    print("mlp")
    print(summary(rsf_data$blt_mlp))
    tictoc::toc()
  }
}


# Calculate BLT
rsf_data$blt <- rsf_data$blt_mlt-rsf_data$blt_mlp
summary(rsf_data$blt)

# Check distribution of BLT value
par(mfrow=c(2,1))
hist(rsf_data$blt[rsf_data$case==TRUE], xlim = c(min(rsf_data$blt, na.rm = TRUE), 75))
hist(rsf_data$blt[rsf_data$case==FALSE], xlim = c(min(rsf_data$blt, na.rm = TRUE), 75),
     breaks=50)
par(mfrow=c(1,1))

# Check negative BLT values
BLT_neg <- rsf_data[which(rsf_data$blt<0),]

# Replace negative values with NA
rsf_data$blt[which(rsf_data$blt<0)] <- NA

# Re-check distribution of BLT value
par(mfrow=c(2,1))
hist(rsf_data$blt[rsf_data$case==TRUE], xlim = c(min(rsf_data$blt, na.rm = TRUE), 75), breaks=50)
hist(rsf_data$blt[rsf_data$case==FALSE], xlim = c(min(rsf_data$blt, na.rm = TRUE), 75), breaks=200)
par(mfrow=c(1,1))

# Remove env data files once satisfied of the data extraction
rm(BLT_5days, mlt_vals, mlp_vals)
gc()


#--------------------------- SAVE DATA READY FOR RSF --------------------------


save(rsf_data, file="RSF/rsf_data_nrdm50_locvalue_2.RData")
load(file="RSF/rsf_data_nrdm50_locvalue_2.RData")



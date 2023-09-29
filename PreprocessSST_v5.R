###############################################################################
#       Preprocess light data from GLS tags for estimation of locations 
#                             using geolocation 
###############################################################################

# Parts of the script by Kirsty Franklin


# clear console and environment
remove(list = ls())
cat("\014")

# Load packages
#install.packages("devtools")
#devtools::install_github("SWotherspoon/SGAT")
library(SGAT)
#devtools::install_github("SWotherspoon/BAStag")
library(BAStag)
library(raster)
library(maptools)
#devtools::install_github("SLisovski/TwGeos")
library(TwGeos) 
library(GeoLight)
library(ggplot2)
library(cowplot)
library(lubridate)

source("GLS_processing/SGAT/colours.r")
source("GLS_processing/SGAT/mask.r")

# Boundary box (for mapping)
xlim0 <- c(70,150)
ylim0 <- c(-40,20)

# tag recording frequency (in seconds)
dt <- 600
# time offset for light data plotting
offset <- 3


#------------------------------------------------------------------------------
#                     TAG METADATA + GENERAL PARAMETERS 
#------------------------------------------------------------------------------


# Set the tag id and the base filename for data files 
# Tagdata different from tag if original file fragmented in separate tracks

tag <- "0585"
tagdata <- "0585"

# Load index file with metadata for all GLS
GLS_index <- read.csv("Abbotts_GLS_Index.csv", header = TRUE)

# Retrieve model
GLS_model <- GLS_index[GLS_index$GLS_ID==tagdata,"GLS_model"]

# Retrieve deployment period (for the tagdata)
Start.Date <- as.POSIXct(paste(paste0(
  GLS_index[GLS_index$GLS_ID==tagdata,"Deploy_Date"], 
  GLS_index[GLS_index$GLS_ID==tagdata,"Deploy_GMT"])), 
  format="%d/%m/%Y %H:%M", tz="GMT")
End.Date <- as.POSIXct(paste(paste0(
  GLS_index[GLS_index$GLS_ID==tagdata,"Retrieval_Date"], 
  GLS_index[GLS_index$GLS_ID==tagdata,"Retrieval_GMT"])), 
  format="%d/%m/%Y %H:%M", tz="GMT")
# Retrieve calibration period (for the tag)
Calib.Start <- as.POSIXct(paste(paste0(
  GLS_index[GLS_index$GLS_ID==tagdata,"Calib_start_date"], 
  GLS_index[GLS_index$GLS_ID==tagdata,"Calib_start_time_GMT"])), 
  format="%d/%m/%Y %H:%M", tz="GMT")
Calib.End <- as.POSIXct(paste(paste0(
  GLS_index[GLS_index$GLS_ID==tagdata,"Calib_stop_date"], 
  GLS_index[GLS_index$GLS_ID==tagdata,"Calib_stop_time_GMT"])), 
  format="%d/%m/%Y %H:%M", tz="GMT")

# If prefer to set dates of deployment start and end manually
# Start.Date<-as.POSIXct("2011-07-20 12:00:00", "GMT")
# End.Date<-as.POSIXct("2012-06-27 12:00", "GMT")


# Define known locations
calib <- if (is.na(GLS_index[GLS_index$GLS_ID==tagdata,"lat_calib"])) {
  matrix(c(105.65549,-10.56142), 1, 2, byrow=T)
} else { matrix(c(GLS_index[GLS_index$GLS_ID==tagdata,"lon_calib"],
                  GLS_index[GLS_index$GLS_ID==tagdata,"lat_calib"]), 
                1, 2, byrow=T)} # calibration site 
locations <- if (is.na(GLS_index[GLS_index$GLS_ID==tagdata,"lat_deploy"])) {
  matrix(c(105.65549,-10.56142), 1, 2, byrow=T)
} else { matrix(c(GLS_index[GLS_index$GLS_ID==tagdata,"lon_deploy"],
                  GLS_index[GLS_index$GLS_ID==tagdata,"lat_deploy"]), 
                1, 2, byrow=T)} # Deployment/retrieval site


# Twilight error distribution (alpha), used in Essie and Estelle models
# Set to resemble the distribution of an open habitat species (Merkel et al. 2016)
alpha <- c(2.49,0.94) 


# Bird speed (beta), set gamma distribution (shape, rate)
beta <- c(3.72, 0.11) 
# Visualise gamma distribution of bird speed
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", 
        xlab = "km/h")

# ADD INFO THERE
threshold <- if(GLS_model %in% c("MK4","MK7")){2.5
} else {if(GLS_model %in% c("MK15","MK3006")){1} else {} # Adjust as needed
}


#------------------------------------------------------------------------------
#                            TWILIGHT ANNOTATION 
#------------------------------------------------------------------------------


# Read light data
d.lig <- readLig(paste0("GLS_files/", tag, "/", tag, "_000.lig"))

# Crop any erroneous data
d.lig <- subset(d.lig, !is.na(Date))

# Visualise data
lightImage(d.lig, offset = offset, zlim = c(0,64))
tsimageDeploymentLines(d.lig$Date, locations[1,1], locations[1,2], 
                       offset = offset, lwd = 2, 
                       col = adjustcolor("orange", alpha.f = 0.6))

# Select deployment period
d.lig <- subset(d.lig, Date >= Start.Date & Date <= End.Date)

# Visualise subset
lightImage(d.lig, offset = offset, zlim = c(0,64))
tsimageDeploymentLines(d.lig$Date, locations[1,1], locations[1,2], 
                       offset = offset, lwd = 2, 
                       col = adjustcolor("orange", alpha.f = 0.6))

# How many days of data? 
length(sort(unique(as.Date(d.lig$Date))))

# Process light data ("twilight annotation") - INTERACTIVE PROCESS!
# see ?preprocessLight for details of the process and commands
twl <- preprocessLight(d.lig,threshold,offset=offset, lmax = 64)

head(twl)

## Investigate twilights less than 3 hrs apart
which(diff(as.numeric(twl$Twilight)) < 3*60*60)

# Investigate rise/set out of sequence
which(diff(twl$Rise)==0)

## Remove deleted
twl <- twl[!twl$Deleted,]

## Store twilights
save(twl,file=paste0("GLS_files/", tag, "/", tagdata,"twl_SGAT.RData"))


#------------------------------------------------------------------------------
#                                CALIBRATION 
#------------------------------------------------------------------------------


#------------- ZENITH FOR INTIAL PATH WITH SIMPLE THRESHOLD METHOD ------------

# Load twilight (object name = twl)
load(file=paste0("GLS_files/", tag, "/", tagdata,"twl_SGAT.RData"))
twl <- twilightAdjust(twl,600)

# Set zenith angles to test
zeniths <- if(GLS_model=="MK4"){seq(94,98,by=0.3)} else {
  if(GLS_model=="MK7"){seq(92,96,by=0.3)} else {
    if(GLS_model=="MK15"){seq(91,95,by=0.3)} else {if(GLS_model=="MK3006"){
      seq(92,96,by=0.3)} else {"ISSUE"}}} #!ADJUST IF NEEDED!
}

# Map path with simple threshold method at different zenith angles
par(mfrow=c(3,5))
data(wrld_simpl)
for(i in 1:length(zeniths)){
  zenith_temp <- zeniths[i]
  path <- thresholdPath(twl$Twilight,twl$Rise,zenith=zenith_temp)
  plot(path$x,type="n", xlim=c(90,150), ylim=c(-50,30), 
       main=paste0("zenith = ", zenith_temp), cex.main = 2.5)
  plot(wrld_simpl,add=T,col=map1.col,border=map2.col)
  plot(elide(wrld_simpl,shift=c(360,0)),add=T,col=map1.col,border=map2.col)
  lines(path$x,col=trk.col)
  points(path$x,pch=16,cex=0.5,col=trk.col)
  points(locations, pch=16, cex=1.25)
  box()
}
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste0("GLS: ",tagdata,"\n","Threshold: ",threshold), 
     cex = 2.5, col = "black")
par(mfrow=c(1,1))

# Plot latitude over time
par(mfrow=c(3,5))
for(i in 1:length(zeniths)){
  zenith_temp <- zeniths[i]
  path <- thresholdPath(twl$Twilight,twl$Rise,zenith=zenith_temp)
  plot(path$time, path$x[,2], main=paste0("zenith = ", zenith_temp), 
       cex.main = 2.5, xlab = "time", ylab = "lat", ylim=c(-50,30))
  abline(h=locations[2])
  date_temp <- as.data.frame(path$time)
  colnames(date_temp)[1] <- "time" 
  date_temp$yday <- yday(as.Date(date_temp$time))
  date_temp$eq[date_temp$yday %in% c(79,265)] <- 1 
  date_temp$eq2[date_temp$yday %in% c(79-21,79+21,265-21,265+21)] <- 1
  abline(v=date_temp$time[date_temp$eq==1], col="red")
  abline(v=date_temp$time[date_temp$eq2==1], col="blue")
}
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste0("GLS: ",tagdata,"\n","Threshold: ",threshold), 
     cex = 2.5, col = "black")
par(mfrow=c(1,1))

# Define zenith based on map and latitude over time 
ZEN_map <- 93.4 # ADJUST!
ZEN_lat <- 93.4 # ADJUST!
zenith <- round(mean(c(ZEN_map, ZEN_lat)),1)

#----------------- ON-BIRD CALIBRATION: ZENITH ZERO DEVIATION -----------------


#### SUBSET PERIOD FOR ON-BIRD CALIBRATION ####


# Read light and activity data
d.lig <- readLig(paste0("GLS_files/", tag, "/", tag, "_000.lig"))
d.act <- readAct2(paste0("GLS_files/", tag, "/", tag, "_000.act"))

# Crop any erroneous data
d.lig <- subset(d.lig, !is.na(Date))
d.act <- subset(d.act, !is.na(Date) & !is.na(Activity))

# Subset deployment period
d.lig <- subset(d.lig, Date >= Start.Date & Date <= End.Date)

# Check the light data
lightImage(d.lig, offset = offset, zlim = c(0,64))
tsimageDeploymentLines(d.lig$Date, locations[1,1], locations[1,2], 
                       offset = offset, lwd = 2, 
                       col = adjustcolor("orange", alpha.f = 0.6))

# Define on-bird calib period (adjust to period of interest)
subset_dates <- as.POSIXct(c("2008-06-18","2008-09-30"), tz="GMT")
abline(v = subset_dates, lwd = 2, lty = 2, col = "orange")

# Subset light data 
d.lig.subset <- subset(d.lig, Date>=subset_dates[1] & Date<=subset_dates[2])

# Define twilight events and Check night activity data (is it mostly dry?)
# IF NO IMMERSION DATA GO TO LINE 342
temp_twilight <- preprocessLight(d.lig.subset, threshold=threshold, offset=offset)

temp_twilight <- twilightAdjust(temp_twilight, 10*60)
temp_twl_gl <- export2GeoLight(temp_twilight) 
temp_twl_gl <- temp_twl_gl[temp_twl_gl$type==2,]

night_act <- NULL
for(i in 1:nrow(temp_twl_gl)){
  temp_night <- data.frame(temp_twl_gl$tFirst[i], temp_twl_gl$tSecond[i])
  subset_act <- d.act[d.act$Date >= temp_night[1,1] & d.act$Date <= temp_night[1,2],]
  temp_night$activity <- sum(subset_act$Activity)
  night_act <- rbind(night_act, as.data.frame(temp_night))
}
colnames(night_act)[1:2] <- c("tFirst", "tSecond")

# Plot night activity data over time (most points must be 0)
plot(night_act$tFirst, night_act$activity*3/60/60, main="Night immersion", 
     xlab="Time", ylab="Immersion (hours)", pch=16, cex=0.7, 
     col=ifelse(night_act$activity > 0, "red", "black"))
(dry_nights <- length(which(night_act$activity == 0))/nrow(night_act)*100)

# IF PERCENTAGE OF DRY NIGHTS TOO LOW (VALUE?) THEN RECONSIDER DATE RANGE, 
# REPLOT LIGHT AND ACTIVITY DATA, AND RESUBSET LIGHT DATA


#### MARK TWILIGHTS ASSOCIATED WITH WET NIGHTS ####


# Sum immersion data for night locations (SKIP SECTION IF NO IMMERSION DATA)
for(i in 1:nrow(temp_twilight)){
  if(temp_twilight$Rise[i] == FALSE){
    temp_night <- c(temp_twilight$Twilight[i], temp_twilight$Twilight[i+1])
    subset_act <- d.act[d.act$Date >= temp_night[1] & d.act$Date <= temp_night[2],]
    temp_twilight$activity[i] <- sum(subset_act$Activity)
  } else { temp_twilight$activity[i] <- NA}
}
# Mark wet nights
for(i in 1:nrow(temp_twilight)){
  if(temp_twilight$Rise[i] == FALSE){
    if(temp_twilight$activity[i] == 0){temp_twilight$wet_rmv[i] <- 0
    } else {temp_twilight$wet_rmv[i] <- 1}
  } else {}
}
# Mark other coordinates derived from twilight of a wet night 
for(i in 1:nrow(temp_twilight)){
    if(temp_twilight$Rise[i]==TRUE){
      if((temp_twilight$wet_rmv[i-1]==0 && temp_twilight$wet_rmv[i+1]==0) ||
         (is.na(temp_twilight$wet_rmv[i-1]) && temp_twilight$wet_rmv[i+1]==0) ||
         (is.na(temp_twilight$wet_rmv[i+1]) && temp_twilight$wet_rmv[i-1]==0)){
        temp_twilight$wet_rmv[i] <- 0} else {temp_twilight$wet_rmv[i] <-1} 
    } else {}
}


#### DEFINE ZENITH 0 DEVIATION ####


# Index dates associated with wet nights (IF NO IMMERSION DATA GO TO LINE 349)
wet_night_index <- as.Date(temp_twilight[temp_twilight$wet_rmv==0, "Twilight"])

# Remove light data for these dates
d.lig_calib <- d.lig.subset[as.Date(d.lig.subset$Date) %in% wet_night_index,]

# Plot light data against corresponding zenith angles at deployment site
thresholdCalibrate(d.lig_calib,locations[1,1],locations[1,2],
                   xlim=c(85,105),ylim=c(0,10),pch=16, max.adjust= TRUE, type="b")
abline(h=threshold,v=zenith,lwd=2, col="blue")
zenith0_temp <- zenith+2.2 # ADJUST! 
# Black lines should be crossing on the furthest plausible line to the right.
abline(h=threshold,v=zenith0_temp,lwd=2)  
legend("topright", legend=paste0("z = ", round(zenith,2), " / z0 = ", 
                                 zenith0_temp), bty="n") 


# IF NO IMMERSION DATA:
thresholdCalibrate(d.lig.subset,locations[1,1],locations[1,2],
                   xlim=c(85,105),ylim=c(0,10),pch=16, max.adjust= TRUE, type="b")
abline(h=threshold,v=zenith, lwd=2, col="dark blue")
zenith0_temp <- zenith+2.4 # ADJUST! 
# Black lines should be crossing on the furthest plausible line to the right
abline(h=threshold,v=zenith0_temp, lwd=2)
legend("topright", legend=paste0("z = ", round(zenith,1), " / z0 = ", 
                                 zenith0_temp), bty="n")


###

zenith0 <- zenith0_temp


#------------------------------------------------------------------------------
#                            SET FIXED LOCATIONS 
#------------------------------------------------------------------------------

# Load twilight file (object name = twl)
load(file=paste0("GLS_files/", tag, "/", tagdata,"twl_SGAT.RData"))

# Define cases that are fixed locations at the colony: these are the first 
# and last locations in the file, assuming the bird was deployed on and 
# recovered there. If the tag packed up before the bird was recovered
# do not change the end location to 1 (i.e. skip that bit by hashing it out)

twl$Marker[1]<-1              # fix start location
twl$Marker[nrow(twl)]<-1      # fix end location
fixed<-twl$Marker==1          # create logical vector of fixed (TRUE or FALSE)

# Change the twilight times for the fixed locations to that predicted for the
# colony on that day by the astronomical algorithms rather than that logged
# by the tag
twl$Twilight[fixed] <- twilight(twl$Twilight[fixed],locations[1,1], 
                                locations[1,2],twl$Rise[fixed],zenith = zenith)


#------------------------------------------------------------------------------
#                              INITIAL PATH 
#------------------------------------------------------------------------------


# Adjust sunset times by 10 mins: see BAStrack manual for explanation of why
twl <- twilightAdjust(twl,600,fixed=fixed) # fixed twilights unmodified

# Compute initial path
path <- thresholdPath(twl$Twilight,twl$Rise,zenith=zenith)

# Set marked points on path based on start and end locations
n<-nrow(path$x)
path$x[1,1]<-locations[1,1]  #set start lon
path$x[1,2]<-locations[1,2]  #set start lat
path$x[n,1]<-locations[1,1]  #set end lon
path$x[n,2]<-locations[1,2]  #set end lat

# Mapping of the initial path
png(paste0("GLS_files/", tag, "/", tagdata, "_initial_path.png"), width = 10, 
    height = 8, units = "in", res = 300)
data(wrld_simpl)
plot(path$x,type="n", xlim=xlim0, ylim=ylim0, main=paste0("zenith = ", zenith))
plot(wrld_simpl,add=T,col=map1.col,border=map2.col)
plot(elide(wrld_simpl,shift=c(360,0)),add=T,col=map1.col,border=map2.col)
lines(path$x,col=trk.col)
points(path$x,pch=16,cex=0.5,col=trk.col)
points(locations, pch=16, cex=1.25)
box()
dev.off()

#------------------------------------------------------------------------------
#                                 SST DATA 
#------------------------------------------------------------------------------


# Load SST functions
source("SGAT/SST.r")

# Import the temperature data
d.tem <- readTem2(paste0("GLS_files/", tag, "/", tag, "_000.tem"))

# Check for temperature data
any(is.na(d.tem$Temp))

# Crop any erroneous data
d.tem <- subset(d.tem, !is.na(Temp))
d.tem <- subset(d.tem, !Valid=="SUSPECT")

# Select deployment period
d.tem <- subset(d.tem, Date > Start.Date & Date < End.Date)

# This function selects which SST values are unreliable by computing the
# median SST for the interval spanning a number of points around each
# point, and then marking points that lie limit degrees above the
# median as suspect.
del <- sstFilter(d.tem,before=3*24,after=3*24,limit=3)

# The set of observations marked for deletion by sstFilter can be
# interactively edited with the selectData function.
del <- selectData(d.tem$Date,d.tem$Temp,deleted=del,pch=16)

# This function searches the SST time series and returns mean SST in
# the interval surrounding each twilight, or NA if there are none.
twl$SST <- twilightData(d.tem$Date[!del],d.tem$Temp[!del],twl$Twilight,
                        before=4,after=4)

# Store twilights with SST
save(twl,file=paste0("GLS_files/", tag, "/", tagdata,"_twlSST_SGAT.RData"))


#------------------------------------------------------------------------------
#                         WRITE CONFIGURATION FILES 
#------------------------------------------------------------------------------


# The final stage is to output all the information needed to analyse
# the preprocessed twilights.  Three files are produced

# path.csv     - the initial path
# twl.csv      - the estimated twilights
# config.r     - an R script to set parameters

# Write path and twilights
colnames(path$x) <- c("lon","lat")
write.csv(cbind.data.frame(time=path$time,path$x),paste0(
  "GLS_files/", tag, "/", tagdata,"_path.csv"),row.names=FALSE)
write.csv(twl,paste0("GLS_files/", tag, "/", tagdata,"_twilight.csv"),
          row.names=FALSE)


# Estimate bounding box from the best estimated path
xlim <- c(floor(min(path$x[,1]))-10,ceiling(max(path$x[,1]))+10)
ylim <- c(max(floor(min(path$x[,2]))-10,-85),min(ceiling(max(path$x[,2]))+10,85))


# Write configuration data
cat("## Configuration\n",
    "tag <- ", deparse(tag),"\n",
    "tagdata <- ", deparse(tagdata),"\n",
    "alpha <- ", deparse(alpha),"\n",
    "beta <- ", deparse(beta),"\n",
    "offset <- ",offset,"\n",
    "threshold <- ",threshold,"\n",
    "zenith <- ", zenith,"\n",
    "zenith0 <- ", zenith0,"\n",
    "fixed.locations <- ",deparse(locations),"\n",
    "xlim <- ", deparse(xlim),"\n",
    "ylim <- ", deparse(ylim),"\n",
    "xlim0 <- ", deparse(xlim0), "\n",
    "ylim0 <- ", deparse(ylim0), "\n",
    file=paste0("GLS_files/", tag, "/", tagdata,"_config.r"),sep="")


###############################################################################
#                           GLS RELATIVE ACCURACY
#                      Following Halpin et al. (2021)
###############################################################################


rm(list=ls())

library(TwGeos)
library(lubridate)

# Load and format data
GLS_index <- read.csv(file="NB_metrics_all_tracks.4.csv", header = TRUE, 
                      sep = ",", na = "NA")
crds_temp <- read.csv(file = "crds_all_pts_land_adj.4.csv", header = TRUE, 
                      sep = ",", na = "NA")
crds_temp$Date <- as.POSIXct(crds_temp$Date)
GLS_index$dep <- as.POSIXct(GLS_index$dep)
GLS_index$ret <- as.POSIXct(GLS_index$ret)

# Set parameters for equinox periods
spring_equinox <- c((79-21):(79+21))
autumn_equinox <- c((265-21):(265+21))
eq_removal <- c(spring_equinox, autumn_equinox)


# Calculate the relative accuracy
relative_accuracy <- NULL
for(i in 1:nrow(GLS_index)){
  ID_temp <- GLS_index$ID[i]
  ID_GLS <- gsub("_1","",ID_temp)
  ID_GLS <- gsub("_2","",ID_GLS)
  # Load twilight file
  twilight_temp <- read.csv(file = paste0(
    "GLS_files/",ID_GLS,"/",ID_temp,"_twilight.csv"), header = TRUE, 
    sep = ",", na = "NA")
  # Subset twilight file for NB trip
  twilight_temp$Twilight <- as.POSIXct(twilight_temp$Twilight, tz="GMT")
  subset_data <- twilight_temp[twilight_temp$Twilight >= 
                                 crds_temp$Date[crds_temp$ID==ID_temp][1] & 
                                 twilight_temp$Twilight <=
                                 crds_temp$Date[crds_temp$ID==ID_temp]
                               [length(crds_temp$Date[crds_temp$ID==ID_temp])],]
  # Convert twilight file format
  twl.gl  <- export2GeoLight(subset_data)
  # Calculate time between two consecutive twilights
  twl.gl$d <- as.numeric(twl.gl$tSecond-twl.gl$tFirst)
  # Mark rows with at date within 21 days of an equinox
  for(i in 1:nrow(twl.gl)){
    if (yday(twl.gl$tFirst[i]) %in% eq_removal) { twl.gl$eq_rmv[i] <- 1} 
    else { twl.gl$eq_rmv[i] <- 0}
  }
  # Print min and max (check not >16 and <9)
  print(summary(twl.gl$d))
  # Calculate relative accuracy (formula from Halpin et al. (2021))
  twl.gl$r_acc <- exp(-0.5*((twl.gl$d - 12)/1.2)^2)
  # Plot relative accuracy over time
  plot(twl.gl$tFirst,twl.gl$r_acc, type="l")
  # Store mean and sd relative accuracy in and out of equinox periods
  relative_accuracy <- rbind(relative_accuracy, 
                             c(ID_temp, mean(twl.gl$r_acc[twl.gl$eq_rmv==0], 
                                             na.rm=TRUE),
                               sd(twl.gl$r_acc[twl.gl$eq_rmv==0], na.rm=TRUE),
                               mean(twl.gl$r_acc[twl.gl$eq_rmv==1], na.rm=TRUE),
                               sd(twl.gl$r_acc[twl.gl$eq_rmv==1], na.rm=TRUE)))
}
rm(i,ID_temp,subset_data,twl.gl, twilight_temp, ID_GLS)

# Reformat results
relative_acc_df <- as.data.frame(relative_accuracy)
rm(relative_accuracy)
colnames(relative_acc_df) <- c("ID","OUT_eq_mean", "OUT_eq_sd", 
                               "IN_eq_mean","IN_eq_sd")
# Means without SDs based only on one value (remove)
relative_acc_df$IN_eq_mean[is.na(relative_acc_df$IN_eq_sd)]<- NA 
relative_acc_df$OUT_eq_mean <- as.numeric(relative_acc_df$OUT_eq_mean)
relative_acc_df$OUT_eq_sd <- as.numeric(relative_acc_df$OUT_eq_sd)
relative_acc_df$IN_eq_mean <- as.numeric(relative_acc_df$IN_eq_mean)
relative_acc_df$IN_eq_sd <- as.numeric(relative_acc_df$IN_eq_sd)


# Calculate global mean
mean(relative_acc_df$OUT_eq_mean)
mean(relative_acc_df$OUT_eq_sd)
mean(relative_acc_df$IN_eq_mean, na.rm=TRUE)
mean(relative_acc_df$IN_eq_sd, na.rm=TRUE)



# Save results
write.csv(relative_acc_df, file="GLS_relative_accuracy.csv", 
          row.names = FALSE, quote=FALSE)












###############################################################################
#                           KDE PER BR + RANDOMISATION
###############################################################################


rm(list=ls())


library(sp)
library(rworldmap)
library(raster)
library(rgdal)
library(adehabitatHR)
library(mapview)
library(dplyr)
library(tictoc)
library(ggplot2)


#------------------------------ LOAD / FORMAT DATA ----------------------------


# Load data
crds_temp <- read.csv(file = "crds_all_pts_land_adj.4.csv", header = TRUE, sep = ",", na = "NA") # Points on CI already removed
GLS_summary_data <- read.csv(file="NB_metrics_all_tracks.4.csv", header = TRUE, sep = ",", na = "NA")
GLS_index <- read.csv("Abbotts_GLS_Index.csv", header = TRUE)

# If removing equinox data:
# crds_temp <- crds_temp[crds_temp$eq_rmv==0,]

# If removing incomplete trips
crds_temp <- crds_temp[!crds_temp$ID %in% c("V905026","22584"),]

for(i in 1:nrow(crds_temp)){
  ID_temp <- crds_temp$ID[i]
  crds_temp$BR[i] <- GLS_summary_data$BR[GLS_summary_data$ID==ID_temp]
}
rm(i, ID_temp)

# Check how many tracks per breeding outcome
n_BR <- count(GLS_summary_data[GLS_summary_data$ID %in% crds_temp$ID,], BR)
n_BR
rm(n_BR)

# Removing unknown breeding outcome (n=4)
crds_temp <- crds_temp[!crds_temp$BR=="Unknown",]
unique(crds_temp$BR)

# Format data
data<-data.frame(ID_animal = as.vector(crds_temp$animal_ID), X = as.vector(crds_temp$x_adj), 
                 Y = as.vector(crds_temp$y_adj), timestamp = paste0(crds_temp$Date, " ", crds_temp$Time), 
                 sex = as.vector(crds_temp$sex), ID = as.vector(crds_temp$ID),
                 BR = as.vector(crds_temp$BR))
data$timestamp <- as.POSIXct(data$timestamp)
summary(data)
head(data)

spdata <- SpatialPointsDataFrame(coords = cbind(data$X, data$Y), 
                                 data = subset(data, select = c(ID_animal, sex, ID,
                                                                BR, timestamp)),
                                 proj4string = CRS("+init=epsg:3857"))

# Check columns of the data to use in kernelUD function
head(spdata)
# Check spatial distribution of the data
plot(spdata)


#------------------------------ CREATE 50*50km GRID ---------------------------


# Define 50*50km grid
cell_size <- 50*1000
x <- seq(xmin(spdata)-(50*cell_size),xmax(spdata)+(50*cell_size),by=cell_size)
y <- seq(ymin(spdata)-(50*cell_size),ymax(spdata)+(50*cell_size),by=cell_size)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy) 
projection(xy) <- CRS("+init=epsg:3857")

# Check grid compared to data
plot(xy)
plot(spdata, add=TRUE)

rm(x,y, cell_size)


#------------------- IDENTIFY SINGLE/MULTIPLE TRIPS PER INDIV -----------------


ID_index <- as.data.frame(unique(crds_temp$ID))
colnames(ID_index)[1]<-"ID"
ID_index$ID_GLS <- ID_index$ID
ID_index$ID_GLS <- gsub("_1","",ID_index$ID_GLS)
ID_index$ID_GLS <- gsub("_2","",ID_index$ID_GLS)
for(i in 1:nrow(ID_index)){
  ID_temp <- ID_index$ID_GLS[i]
  ID_index$ID_animal[i] <- GLS_index$Animal_ID_ABBBS[GLS_index$ID==ID_temp]
  ID_track <- ID_index$ID[i]
  ID_index$BR[i] <- GLS_summary_data$BR[GLS_summary_data$ID==ID_track]
}
rm(i, ID_temp)

# Check if bird with several trips in same BR
n_trips_BR <- ID_index %>% count(BR, ID_animal)
n_trips_BR


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
  BR_bird_list <- vector(mode = "list", length = length(unique(subset_BR$ID_animal)))
  # Run KDE per bird
  for(j in 1:length(unique(subset_BR$ID_animal))){
    ID_animal_temp <- unique(subset_BR$ID_animal)[j]
    # Subset data per bird ID
    subset_ID_animal <- subset_BR[subset_BR$ID_animal==ID_animal_temp,]
    # Create empty list to store rasters
    BR_trips_list <- vector(mode = "list", length = length(unique(subset_ID_animal$ID)))
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
       ID_temp,subset_BR,udm,BR,ID_animal_temp,ID_track,udm_birds, subset_ID_animal)
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

BA_95
BA_50


# Save shapefiles
# writeOGR(obj=to_adjust, dsn="GIS_temp", 
#          layer=paste0("to adjust",smooth_h/1000), driver="ESRI Shapefile")
# writeOGR(obj=to_adjust, dsn="GIS_temp", 
#          layer=paste0("to adjust",smooth_h/1000), driver="ESRI Shapefile")
# 

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
    BR_bird_list <- vector(mode = "list", length = length(unique(subset_BR$ID_animal)))
    # Run KDE per bird
    for(j in 1:length(unique(subset_BR$ID_animal))){
      ID_animal_temp <- unique(subset_BR$ID_animal)[j]
      # Subset data per bird ID
      subset_ID_animal <- subset_BR[subset_BR$ID_animal==ID_animal_temp,]
      # Create empty list to store rasters
      BR_trips_list <- vector(mode = "list", length = length(unique(subset_ID_animal$ID)))
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
    # if(i==length(BR_list)){
    #   rm(i,j,BR_trips_list,combined_temp,udm_all,k_temp,k_sub_BR,ID_trip_spdata,
    #      ID_temp,subset_BR,udm,BR,ID_animal_temp,ID_track,udm_birds, subset_ID_animal)
    #}
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
    # rm(i,k,j,l, BA_50_rdm, BA_95_rdm, BR_trips_list, combined_temp, udm_all, k_temp, 
    #    k_sub_BR, ID_trip_spdata, ID_temp,subset_BR,udm,BR, rdm_BR, overlap_50_rdm,
    #    overlap_95_rdm)
  }
}

# Save overlap matrices
write.csv(overlap_matrix, file=paste0("overlap_rdm_BR.csv"), 
          row.names = FALSE, quote=FALSE)


###############################################################################
#                             p-value + plotting
###############################################################################


overlap_matrix <-  read.csv(file = "overlap_rdm_BR.csv", header = TRUE, sep = ",", na = "NA")

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
hist(overlap_matrix$X95., breaks=10, main="95%", xlab="overlap (BA)", xlim = c(0,1))
text(x = 0.2, y = 160, paste0("p=",round(p_95,2)))
abline(v = BA_95[2,1], col="red")

summary(overlap_matrix)
sd(overlap_matrix$X50.)

# Plot for poster


(BR_overlap_plot <- ggplot(overlap_matrix, aes(X50.))+
    theme_bw()+
    geom_histogram(alpha = 1, bins = 25, position = 'identity', fill= "#B7C5E4")+
    geom_vline(alpha=1, xintercept = BA_50[2,1], col = "red", linetype = "solid", size=2.5)+
    ylab(" \n")+
    xlab("\nOverlap (BA)")+
    xlim(0,0.50)+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 25, colour="white"),
          plot.margin = unit(c(0.75,0.75,0.75,0.75), units = , "cm"),
          axis.text.y = element_text(size = 25, colour="white"),
          axis.title.y = element_text(size = 28, colour = "white"),
          axis.title.x = element_text(size = 28, colour = "white"),
          strip.text = element_text(size = 13),
          axis.ticks = element_line(colour = "white", size=1.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "white", size=1.5),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
    )
)


ggsave(BR_overlap_plot, filename = "overlap_BR_trans.png", bg = "transparent", dpi = 500)


# Plot for thesis

(BR_overlap_plot <- ggplot(overlap_matrix, aes(X50.))+
    theme_bw()+
    geom_histogram(alpha = 1, bins = 25, fill= "grey")+
    geom_vline(alpha=1, xintercept = BA_50[2,1], col = "red", linetype = "solid", size=1.5)+
    ylab("Frequency \n")+
    xlab("\nOverlap (BA)")+
    ggtitle("50%")+
    xlim(0,0.50)+
    ylim(0,200)+
    annotate("text", x = 0.1, y = 150, label = paste0("p = ",round(p_50,2)), size=6)+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 16, colour="black"),
          plot.margin = unit(c(0.75,0.75,0.75,0.75), units = , "cm"),
          axis.text.y = element_text(size = 16, colour="black"),
          axis.title.y = element_text(size = 18, colour = "black"),
          axis.title.x = element_text(size = 18, colour = "black"),
          axis.ticks = element_line(colour = "black", size=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", size=1),
          plot.title = element_text(size=22, hjust = 1)
    )
)

(BR_overlap_plot <- ggplot(overlap_matrix, aes(X95.))+
    theme_bw()+
    geom_histogram(alpha = 1, bins = 50, fill= "grey")+
    geom_vline(alpha=1, xintercept = BA_95[2,1], col = "red", linetype = "solid", size=1.5)+
    ylab(" \n")+
    xlab("\nOverlap (BA)")+
    ggtitle("95%")+
    xlim(0,1)+
    ylim(0,300)+
    annotate("text", x = 0.25, y = 230, label = paste0("p = ",round(p_95,2)), size=6)+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 16, colour="black"),
          plot.margin = unit(c(0.75,0.75,0.75,0.75), units = , "cm"),
          axis.text.y = element_text(size = 16, colour="black"),
          axis.title.y = element_text(size = 18, colour = "black"),
          axis.title.x = element_text(size = 18, colour = "black"),
          axis.ticks = element_line(colour = "black", size=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", size=1),
          plot.title = element_text(size=22, hjust = 1)
    )
)


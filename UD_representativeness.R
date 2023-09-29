###############################################################################
#               UTILISATION DISTRIBUTION - REPRESENTATIVENESS
#                    Following Lascelles et al. (2016)
###############################################################################


rm(list=ls())

library(sp)
library(rworldmap)
library(raster)
library(rgdal)
library(adehabitatHR)
library(dplyr)
library(ggplot2)


#------------------------------ LOAD / FORMAT DATA ----------------------------


# Load data
crds_temp <- read.csv(file = "crds_all_pts_land_adj.4.csv", header = TRUE, 
                      sep = ",", na = "NA") 
GLS_index <- read.csv("Abbotts_GLS_Index.csv", header = TRUE)

# Remove incomplete trips
crds_temp <- crds_temp[!crds_temp$ID %in% c("V905026","22584"),]

# Format data
data<-data.frame(ID = as.vector(crds_temp$animal_ID), X = as.vector(crds_temp$x_adj), 
                 Y = as.vector(crds_temp$y_adj), trip_ID = as.vector(crds_temp$ID))

spdata <- SpatialPointsDataFrame(coords = cbind(data$X, data$Y), 
                                 data = subset(data, select = c(ID, trip_ID)),
                                 proj4string = CRS("+init=epsg:3857"))
head(spdata)
plot(spdata)


#--------------------------- CREATE 50*50km GRID ------------------------------


# Define 50*50km grid
cell_size <- 50*1000
x <- seq(xmin(spdata)-(50*cell_size),xmax(spdata)+(50*cell_size),by=cell_size)
y <- seq(ymin(spdata)-(50*cell_size),ymax(spdata)+(50*cell_size),by=cell_size)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy) 
projection(xy) <- CRS("+init=epsg:3857")

plot(xy)
plot(spdata, add=TRUE)


#----------------------- COUNT NUMBER OF TRIPS PER INDIV ----------------------


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


###############################################################################
#                        REPRESENTATIVENESS ANALYSIS
###############################################################################


iterations <- 100
smooth_h <- 186*1000

inclusion_data <- data.frame(matrix(ncol = iterations, nrow = nrow(n_trips)-1))


#------------------------- CALCULATE INCLUSION VALUES -------------------------


for(i in 1:iterations){
  print(paste0(i,"/",iterations))
  for(j in 1:(nrow(n_trips)-1)){
    # randomly sampling j number of indiv
    n_rdm <- sample(n_trips$ID_animal,size=j,replace=FALSE)
    list_rasters <- vector(mode = "list", length = j) # list to add indiv UDs in
    spdata_temp <- spdata[spdata$ID %in% n_rdm,]
    spdata_temp$ID <- "1" 
    k_temp <- kernelUD(spdata_temp[,1], grid=xy, extent=1, h=smooth_h)
    # Calculate UDs for sampled indiv
    for(k in 1:length(n_rdm)){
      spdata_subset <- spdata[spdata$ID==n_rdm[k],] 
      k_multi <- kernelUD(spdata_subset[,3], grid=xy, extent=1, h=smooth_h)
      k_multi_temp <- stack(lapply(k_multi, raster))
      udm <- mean(k_multi_temp)
      udm@data@names <- n_rdm[k]
      list_rasters[[k]] <- udm
    }
    # Average across indiv
    list_stacked <- stack(list_rasters)
    udm_all <- mean(list_stacked)
    k_temp[[1]]@data$ud <- udm_all@data@values[k_temp[[1]]@grid.index]
    iso50 <- getverticeshr(k_temp,50)
    # Subset data not used for UD calculation
    spdata_unsampled <- spdata[!spdata$ID %in% n_rdm,]
    # Calculate proportion of unsampled data within 50% isopleth
    pts_in_50 <- over(spdata_unsampled, iso50, fn = NULL)
    pts_in_50 <- ifelse(is.na(pts_in_50$id),0,1)
    prop_temp <- sum(pts_in_50)/nrow(spdata_unsampled)
    # Add data to data
    inclusion_data[j,i] <- prop_temp 
  }
}

rm(prop_temp, pts_in_50, spdata_unsampled, k_temp, iso50, list_rasters, udm, 
   k_multi_temp, k_multi, spdata_subset, k, i, j, udm_all, list_stacked, n_rdm,
   spdata_temp)


# Save results
write.csv(inclusion_data, file="inclusion_random_100_iter.csv", 
          row.names = FALSE, quote=FALSE)


#--------------------- CALCULATE REPRESENTATIVENESS ---------------------------


inclusion_data <- read.csv(file = "inclusion_random_100_iter.csv", 
                           header = TRUE, sep = ",", na = "NA")


# Calculate mean and sd for each sample size
for(i in 1:nrow(inclusion_data)){
  inclusion_data$mean[i] <- mean(unname(unlist(inclusion_data[i,1:iterations])))
  inclusion_data$sd[i] <- sd(unname(unlist(inclusion_data[i,1:iterations])))
}
# Add sample size
inclusion_data$sample_size <- seq(1,21,1)

# Fit non linear regression
x_nls <- inclusion_data$sample_size
y_nls <- inclusion_data$mean
n<-nls(y_nls~(a*x_nls)/(b+x_nls))  # Michaelis-Menten equation - "a" is asymptote

# Check goodness of fit
cor(y_nls,predict(n))
plot(y_nls, predict(n))

# Check results
summary(n)
summary(predict(n))

# Calculate predicted value for our whole dataset (22 birds)
pred <- (summary(n)$coefficients[1]*nrow(n_trips))/(summary(n)$coefficients[2]+nrow(n_trips))

# Calculate representativeness: highest predicted value / asymptote of non-linear regression model
represent <- as.numeric(pred/summary(n)$coefficients[1])*100

# Plot mean inclusion data and prediction from non linear model
(represent_plot <- ggplot(inclusion_data, aes(sample_size, mean))+
    geom_point(shape=1, size=2)+
    geom_ribbon(aes(ymin = (mean-sd), ymax = (mean+sd)), alpha = 0.15)+
    geom_line(aes(sample_size, predict(n)), size=1)+
    theme_bw()+
    ylab("Mean inclusion value\n")+
    xlab("\nNumber of individuals sampled")+
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
          plot.title = element_text(size=22, hjust = 1))+
    annotate("text", x = 12, y = 0.3, 
             label = paste0(round(represent,2),"%"), color="black", 
             size=6)
)  

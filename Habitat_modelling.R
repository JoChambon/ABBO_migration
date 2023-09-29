###############################################################################
#                         Resource Selection Function
#                                 w/ GAMMs
###############################################################################

rm(list=ls())

library(ggplot2)
library(ggpubr)
library(mgcv)
library(nlme)
library(gratia)
library(pROC)


load(file="RSF/rsf_data_nrdm50_locvalue_2.RData")


#------------------------------------------------------------------------------


# For all data
data_temp <- subset(rsf_data, !is.na(sst) & !is.na(chl) & !is.na(sal) & bat<0 & !is.na(slope_200) 
                    & !is.na(mlt) & !is.na(sst_grad) & !is.na(eke) & !is.na(slope) & !is.na(blt))
row.names(data_temp) <- NULL


# ----------------------- FINDING MOST PREDICTIVE VAR -------------------------


# SST
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sst_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sst, k=4, bs="cs")+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sst_CV <- rbind(AUC_sst_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sst_mean <- mean(AUC_sst_CV[,1])
  }else{}
}

# CHLOROPHYLL a
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_chl_log_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(chl_log, k=4, bs="cs")+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_chl_log_CV <- rbind(AUC_chl_log_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_chl_log_mean <- mean(AUC_chl_log_CV[,1])
  }else{}
}

# CHLOROPHYLL a
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_chl_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(chl, k=4, bs="cs")+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_chl_CV <- rbind(AUC_chl_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_chl_mean <- mean(AUC_chl_CV[,1])
  }else{}
}

#BATHYMETRY
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_bat_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(bat, k=4, bs="cs")+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_bat_CV <- rbind(AUC_bat_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_bat_mean <- mean(AUC_bat_CV[,1])
  }else{}
}

# SALINITY
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sal, k=4, bs="cs")+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_CV <- rbind(AUC_sal_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_mean <- mean(AUC_sal_CV[,1])
  }else{}
}

# MLT
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_mlt_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(mlt, k=4, bs="cs")+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_mlt_CV <- rbind(AUC_mlt_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_mlt_mean <- mean(AUC_mlt_CV[,1])
  }else{}
}

# EKE
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_eke_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(eke, k=4, bs="cs")+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_eke_CV <- rbind(AUC_eke_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_eke_mean <- mean(AUC_eke_CV[,1])
  }else{}
}

# SLOPE 200 km
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_slope_200_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_slope_200_CV <- rbind(AUC_slope_200_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_slope_200_mean <- mean(AUC_slope_200_CV[,1])
  }else{}
}

# BLT
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_blt_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(blt, k=4, bs="cs")+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_blt_CV <- rbind(AUC_blt_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_blt_mean <- mean(AUC_blt_CV[,1])
  }else{}
}


# Saving AUC values from CV for most predictive variable
#save(AUC_sal_CV, file="RSF/AUC_sal_k4.RData")


#----------------------------- TESTING ADDING SEX -----------------------------


#SLOPE + SEX
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_slope_sex_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope, k=5, bs="cs")+sex+s(animal_ID, bs = 're'), data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_slope_sex_CV <- rbind(AUC_slope_sex_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_slope_sex_mean <- mean(AUC_slope_sex_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_slope_CV, AUC_slope_sex_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_slope_CV, AUC_slope_sex_CV, paired=T, alternative = c("less"))


#---------------------------- TESTING 2ND VARIABLE ----------------------------


load(file="RSF/AUC_sal_k4.RData")

#SST
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_sst_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_sst_CV <- rbind(AUC_sal_sst_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_sst_mean <- mean(AUC_sal_sst_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_CV, AUC_sal_sst_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_CV, AUC_sal_sst_CV, paired=T, alternative = c("less"))

#sst_grad grad
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_sst_grad_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sal, k=4, bs="cs")+s(sst_grad, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_sst_grad_CV <- rbind(AUC_sal_sst_grad_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_sst_grad_mean <- mean(AUC_sal_sst_grad_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_CV, AUC_sal_sst_grad_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_CV, AUC_sal_sst_grad_CV, paired=T, alternative = c("less"))

#CHLOROPHYLL a
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_chl_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sal, k=4, bs="cs")+s(chl, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_chl_CV <- rbind(AUC_sal_chl_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_chl_mean <- mean(AUC_sal_chl_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_CV, AUC_sal_chl_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_CV, AUC_sal_chl_CV, paired=T, alternative = c("less"))


#BATHYMETHRY
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_bat_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sal, k=4, bs="cs")+s(bat, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_bat_CV <- rbind(AUC_sal_bat_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_bat_mean <- mean(AUC_sal_bat_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_CV, AUC_sal_bat_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_CV, AUC_sal_bat_CV, paired=T, alternative = c("less"))


#SLOPE (200 km radius)
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sal, k=4, bs="cs")+s(slope_200, k=5, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_CV <- rbind(AUC_sal_slope_200_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_mean <- mean(AUC_sal_slope_200_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_CV, AUC_sal_slope_200_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_CV, AUC_sal_slope_200_CV, paired=T, alternative = c("less"))


#MIXED LAYER DEPTH
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_mlt_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sal, k=4, bs="cs")+s(mlt, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_mlt_CV <- rbind(AUC_sal_mlt_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_mlt_mean <- mean(AUC_sal_mlt_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_CV, AUC_sal_mlt_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_CV, AUC_sal_mlt_CV, paired=T, alternative = c("less"))


#EDDY KINETIC ENERGY (EKE)
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_eke_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sal, k=4, bs="cs")+s(eke, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_eke_CV <- rbind(AUC_sal_eke_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_eke_mean <- mean(AUC_sal_eke_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_CV, AUC_sal_eke_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_CV, AUC_sal_eke_CV, paired=T, alternative = c("less"))


# BLT
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_blt_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(sal, k=4, bs="cs")+s(blt, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_blt_CV <- rbind(AUC_sal_blt_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_blt_mean <- mean(AUC_sal_blt_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_CV, AUC_sal_blt_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_CV, AUC_sal_blt_CV, paired=T, alternative = c("less"))


# Saving AUC values from CV for most predictive variable
save(AUC_sal_slope_200_CV, file="RSF/AUC_sal_slope_200_CV_adj.RData")


#-------------------------- TESTING 3RD VARIABLE ------------------------------


load(file="RSF/AUC_sal_slope_200_CV_adj.RData")


#SST
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_CV <- rbind(AUC_sal_slope_200_sst_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_mean <- mean(AUC_sal_slope_200_sst_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_sst_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_sst_CV, paired=T, alternative = c("less"))

#SST gradient
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_grad_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst_grad, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_grad_CV <- rbind(AUC_sal_slope_200_sst_grad_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_grad_mean <- mean(AUC_sal_slope_200_sst_grad_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_sst_grad_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_sst_grad_CV, paired=T, alternative = c("less"))

#CHLOROPHYLL a
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_chl_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(chl, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_chl_CV <- rbind(AUC_sal_slope_200_chl_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_chl_mean <- mean(AUC_sal_slope_200_chl_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_chl_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_chl_CV, paired=T, alternative = c("less"))

#BATHYMETRY
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_bat_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(bat, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_bat_CV <- rbind(AUC_sal_slope_200_bat_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_bat_mean <- mean(AUC_sal_slope_200_bat_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_bat_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_bat_CV, paired=T, alternative = c("less"))

#EDDY KINETIC ENERGY (EKE)
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_eke_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(eke, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_eke_CV <- rbind(AUC_sal_slope_200_eke_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_eke_mean <- mean(AUC_sal_slope_200_eke_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_eke_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_eke_CV, paired=T, alternative = c("less"))

#MIXED LAYER DEPTH
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_mlt_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(mlt, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_mlt_CV <- rbind(AUC_sal_slope_200_mlt_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_mlt_mean <- mean(AUC_sal_slope_200_mlt_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_mlt_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_mlt_CV, paired=T, alternative = c("less"))


# BLT
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_blt_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(blt, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_blt_CV <- rbind(AUC_sal_slope_200_blt_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_blt_mean <- mean(AUC_sal_slope_200_blt_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_blt_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_CV, AUC_sal_slope_200_blt_CV, paired=T, alternative = c("less"))


save(AUC_sal_slope_200_sst_CV, file="RSF/AUC_sal_slope_200_sst_CV_adj.RData")
save(AUC_sal_slope_200_chl_CV, file="RSF/AUC_sal_slope_200_chl_CV_adj.RData")
save(AUC_sal_slope_200_bat_CV, file="RSF/AUC_sal_slope_200_bat_CV_adj.RData")
save(AUC_sal_slope_200_eke_CV, file="RSF/AUC_sal_slope_200_eke_CV_adj.RData")
save(AUC_sal_slope_200_mlt_CV, file="RSF/AUC_sal_slope_200_mlt_CV_adj.RData")


#-------------------------- TESTING 4th VARIABLE ------------------------------


load(file="RSF/AUC_sal_slope_200_sst_CV_adj.RData")


#CHLOROPHYLL a
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_chl_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(chl, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_chl_CV <- rbind(AUC_sal_slope_200_sst_chl_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_chl_mean <- mean(AUC_sal_slope_200_sst_chl_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_chl_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_chl_CV, paired=T, alternative = c("less"))

#BATHYMETRY
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_bat_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(bat, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_bat_CV <- rbind(AUC_sal_slope_200_sst_bat_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_bat_mean <- mean(AUC_sal_slope_200_sst_bat_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_bat_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_bat_CV, paired=T, alternative = c("less"))

#EDDY KINETIC ENERGY (EKE)
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_eke_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(eke, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_eke_CV <- rbind(AUC_sal_slope_200_sst_eke_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_eke_mean <- mean(AUC_sal_slope_200_sst_eke_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_eke_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_eke_CV, paired=T, alternative = c("less"))

#MIXED LAYER DEPTH
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_mlt_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(mlt, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_mlt_CV <- rbind(AUC_sal_slope_200_sst_mlt_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_mlt_mean <- mean(AUC_sal_slope_200_sst_mlt_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_mlt_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_mlt_CV, paired=T, alternative = c("less"))


# BLT
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_blt_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(blt, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_blt_CV <- rbind(AUC_sal_slope_200_sst_blt_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_blt_mean <- mean(AUC_sal_slope_200_sst_blt_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_blt_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_CV, AUC_sal_slope_200_sst_blt_CV, paired=T, alternative = c("less"))


save(AUC_sal_slope_200_sst_chl_CV, file="RSF/AUC_sal_slope_200_sst_chl_CV_adj.RData")
save(AUC_sal_slope_200_sst_bat_CV, file="RSF/AUC_sal_slope_200_sst_bat_CV_adj.RData")
save(AUC_sal_slope_200_sst_eke_CV, file="RSF/AUC_sal_slope_200_sst_eke_CV_adj.RData")
save(AUC_sal_slope_200_sst_mlt_CV, file="RSF/AUC_sal_slope_200_sst_mlt_CV_adj.RData")


#-------------------------- TESTING 5th VARIABLE ------------------------------


load(file="RSF/AUC_sal_slope_200_sst_mlt_CV_adj.RData")


#CHLOROPHYLL a
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_mlt_chl_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(mlt, k=4, bs="cs")+s(chl, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_mlt_chl_CV <- rbind(AUC_sal_slope_200_sst_mlt_chl_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_mlt_chl_mean <- mean(AUC_sal_slope_200_sst_mlt_chl_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_sal_slope_200_sst_mlt_chl_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_sal_slope_200_sst_mlt_chl_CV, paired=T, alternative = c("less"))

#BATHYMETRY
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_mlt_bat_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(mlt, k=4, bs="cs")+s(bat, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_mlt_bat_CV <- rbind(AUC_sal_slope_200_sst_mlt_bat_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_mlt_bat_mean <- mean(AUC_sal_slope_200_sst_mlt_bat_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_sal_slope_200_sst_mlt_bat_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_sal_slope_200_sst_mlt_bat_CV, paired=T, alternative = c("less"))

#EDDY KINETIC ENERGY (EKE)
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_mlt_eke_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(mlt, k=4, bs="cs")+s(eke, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_mlt_eke_CV <- rbind(AUC_sal_slope_200_sst_mlt_eke_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_mlt_eke_mean <- mean(AUC_sal_slope_200_sst_mlt_eke_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_sal_slope_200_sst_mlt_eke_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_sal_slope_200_sst_mlt_eke_CV, paired=T, alternative = c("less"))


# BLT
for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_sal_slope_200_sst_mlt_blt_CV <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(mlt, k=4, bs="cs")+s(blt, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_sal_slope_200_sst_mlt_blt_CV <- rbind(AUC_sal_slope_200_sst_mlt_blt_CV, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_sal_slope_200_sst_mlt_blt_mean <- mean(AUC_sal_slope_200_sst_mlt_blt_CV[,1])
    #rm()
  }else{}
}

t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_sal_slope_200_sst_mlt_blt_CV, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_sal_slope_200_sst_mlt_blt_CV, paired=T, alternative = c("less"))


save(AUC_sal_slope_200_sst_mlt_chl_CV, file="RSF/AUC_sal_slope_200_sst_mlt_chl_CV_adj.RData")
save(AUC_sal_slope_200_sst_mlt_bat_CV, file="RSF/AUC_sal_slope_200_sst_mlt_bat_CV_adj.RData")
save(AUC_sal_slope_200_sst_mlt_eke_CV, file="RSF/AUC_sal_slope_200_sst_mlt_eke_CV_adj.RData")


#------------------------------- FINAL MODEL ----------------------------------


m_final <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
               +s(mlt, k=4, bs="cs")+s(animal_ID, bs = 're'), 
               data=data_temp, method = "REML", family = binomial())

save(m_final, file="RSF/model_final_adj.RData")

load(file="RSF/model_final_adj.RData")

summary(m_final)
gam.check(m_final)

par(mfrow=c(2,3))
draw(m_final)
par(mfrow=c(1,1))


# Check collinearity between variable using Variance Inflation Factor (VIF)

library(car)

test_conc <- concurvity(m_final)
draw(test_conc)
concurvity(m_final, full=FALSE)
vif(m_final)



###############################################################################
#                           BIC for all combinations
###############################################################################


globalmodel <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
                   +s(mlt, k=4, bs="cs")+s(chl, k=4, bs="cs")+s(bat, k=4, bs="cs")
                   +s(eke, k=4, bs="cs")+s(animal_ID, bs = 're'), 
                   data=data_temp, method = "ML", family = binomial(), na.action = "na.fail")

save(globalmodel, file="RSF/globalmodel.RData")

combinations <- dredge(globalmodel, rank = "BIC", fixed = ~s(animal_ID, bs = 're'), trace=2)

save(combinations, file="RSF/model_comp_MuMIn.RData")

print(combinations)
coefTable(combinations)


# Check AUC of best models


for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_BIC_model <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(mlt, k=4, bs="cs")+s(chl, k=4, bs="cs")+s(bat, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_BIC_model <- rbind(AUC_BIC_model, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_BIC_model_mean <- mean(AUC_BIC_model[,1])
  }else{}
}

save(AUC_BIC_model, file="RSF/AUC_top_model_BIC.RData")

for(i in 1:length(unique(data_temp$animal_ID))){
  if(i == 1){AUC_BIC_model_2 <- NULL}else{}
  print(paste0(i,"/",length(unique(data_temp$animal_ID))))
  ID_temp <- unique(data_temp$animal_ID)[i]
  # Set training and testing datasets
  train.dat <-  subset(data_temp, !animal_ID==ID_temp)
  test.dat <- subset(data_temp, animal_ID==ID_temp)
  # Fit model and predict with train dataset
  fit <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
             +s(mlt, k=4, bs="cs")+s(chl, k=4, bs="cs")+s(bat, k=4, bs="cs")
             +s(eke, k=4, bs="cs")+s(animal_ID, bs = 're'), 
             data=train.dat, method = "REML", family = binomial())
  predicted <- predict(fit, test.dat, type = "response")
  # ROC curve and AUC
  roc_temp <- roc(test.dat$case ~ predicted)
  AUC_temp <- as.numeric(roc_temp$auc)
  AUC_BIC_model_2 <- rbind(AUC_BIC_model_2, AUC_temp)
  if(i==length(unique(data_temp$animal_ID))){
    AUC_BIC_model_2_mean <- mean(AUC_BIC_model_2[,1])
  }else{}
}

save(AUC_BIC_model_2, file="RSF/AUC_top_model_2_BIC.RData")



# Check best model


m_best <- gam(case~s(slope_200, k=5, bs="cs")+s(sal, k=4, bs="cs")+s(sst, k=4, bs="cs")
              +s(mlt, k=4, bs="cs")+s(chl, k=4, bs="cs")+s(bat, k=4, bs="cs")+s(animal_ID, bs = 're'), 
              data=data_temp, method = "REML", family = binomial())

save(m_best, file="RSF/model_best_BIC.RData")

draw(m_best)


# Compare best model BIC with best model forward selection


load(file="RSF/AUC_top_model_BIC.RData")
load(file="RSF/AUC_top_model_2_BIC.RData")
load(file="RSF/AUC_sal_slope_200_sst_mlt_CV_adj.RData")

t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_BIC_model, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_BIC_model, paired=T, alternative = c("less"))

t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_BIC_model_2, paired=T, alternative = c("two.sided"))
t.test(AUC_sal_slope_200_sst_mlt_CV, AUC_BIC_model_2, paired=T, alternative = c("less"))

t.test(AUC_BIC_model, AUC_BIC_model_2, paired=T, alternative = c("two.sided"))
t.test(AUC_BIC_model, AUC_BIC_model_2, paired=T, alternative = c("less"))




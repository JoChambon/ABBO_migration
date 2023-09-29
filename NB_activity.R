###############################################################################
# ANALYSES OF ACTIVITY
# DURING NON-BREEDING MIGRATION
###############################################################################


rm(list=ls())
library(dplyr)
library(ggplot2)
library(fitdistrplus)
library(lme4)
library(AICcmodavg)
library(brms)
library(boot)
library(DHARMa)


#------------------------------------------------------------------------------
# LOAD AND FORMAT DATA
#------------------------------------------------------------------------------


# Load data
GLS_index <- read.csv("Abbotts_GLS_Index.csv", header = TRUE)
GLS_summary_data <- read.csv(file="NB_metrics_all_tracks.4.csv", header = TRUE,
                             sep = ",", na = "NA")
GLS_act <- read.csv(file = "Actave/actave-data_ABBOTTS.pheno.csv", header = TRUE,
                    sep = ",", na = "NA")

# Add GLS model to GLS_summary_data and NB_act
for(i in 1:nrow(GLS_summary_data)){
  ID_temp <- GLS_summary_data$ID[i]
  GLS_summary_data$GLS_model[i] <- GLS_index$GLS_model[GLS_index$ID==ID_temp]
}
for(i in 1:nrow(GLS_act)){
  ID_temp <- GLS_act$ID[i]
  GLS_act$GLS_model[i] <- GLS_index$GLS_model[GLS_index$ID==ID_temp]
}
rm(i, ID_temp)

# Remove MK7s (MK7 activity data overflowed before/during NB trip)
GLS_act <- GLS_act[!GLS_act$GLS_model=="MK7",]

# Remove GLS without non-breeding trips
GLS_act <- GLS_act[!GLS_act$ID %in% c("22583_2", "22588_2", "22591_2"),]

# Check data
sort(unique(GLS_act$ID))
length(unique(GLS_act$ID))
sort(unique(GLS_summary_data$ID[!GLS_summary_data$GLS_model=="MK7"]))
length(unique(GLS_summary_data$ID[!GLS_summary_data$GLS_model=="MK7"]))
unique(GLS_summary_data$ID[!GLS_summary_data$ID %in% unique(GLS_act$ID) &
                             !GLS_summary_data$GLS_model=="MK7"])

# Subset non-breeding trip data
GLS_act$date <- as.Date(GLS_act$date)
GLS_summary_data$dep <- as.Date(GLS_summary_data$dep)
GLS_summary_data$ret <- as.Date(GLS_summary_data$ret)
str(GLS_act)
NB_act <- NULL
for(i in 1:length(unique(GLS_act$ID))){
  ID_temp <- unique(GLS_act$ID)[i]
  GLS_subset <- GLS_act[GLS_act$ID==ID_temp,]
  dep_date <- GLS_summary_data$dep[GLS_summary_data$ID==ID_temp]
  ret_date <- GLS_summary_data$ret[GLS_summary_data$ID==ID_temp]
  if(is.na(ret_date)==TRUE){
    NB_subset <- GLS_subset[GLS_subset$date > dep_date,]
  }else{
    NB_subset <- GLS_subset[GLS_subset$date > dep_date &
                              GLS_subset$date < ret_date, ]
  }
  NB_act <- rbind(NB_act, NB_subset)
}
rm(NB_subset, GLS_subset,i,ID_temp,dep_date,ret_date)

# Checking if the departure date is after the last record of activity
missing_GLS <- unique(GLS_act$ID)[!unique(GLS_act$ID) %in% unique(NB_act$ID)]
for(i in 1:length(missing_GLS)){
  max_act_date <- max(GLS_act$date[GLS_act$ID==missing_GLS[i]])
  dep <- GLS_summary_data[GLS_summary_data$ID==missing_GLS[i], "dep"]
  result <- ifelse(max_act_date < dep, "NO DATA", "ISSUE")
  print(paste0(missing_GLS[i]," - ",result))
}
rm(dep, max_act_date, result, i, missing_GLS)

# Add explanatory variables to the dataset
data_act <- NB_act
for(i in 1:nrow(data_act)){
  ID_temp <- data_act$ID[i]
  data_act$bird_ID[i] <- GLS_summary_data$ID_animal[GLS_summary_data$ID==ID_temp]
  data_act$sex[i] <- GLS_summary_data$sex[GLS_summary_data$ID==ID_temp]
  data_act$year[i] <- GLS_summary_data$year[GLS_summary_data$ID==ID_temp]
  data_act$BR[i] <- GLS_summary_data$BR[GLS_summary_data$ID==ID_temp]
}
rm(i, ID_temp)

# Reformat data
data_act$sex <- as.factor(data_act$sex)
data_act$year <- as.factor(data_act$year)
data_act$bird_ID <- as.factor(data_act$bird_ID)
data_act$BR <- as.factor(data_act$BR)
data_act$ID <- as.factor(data_act$ID)
data_act$phenophase <- as.factor(data_act$phenophase)
data_act$GLS_model <- as.factor(data_act$GLS_model)
data_act <- droplevels(data_act)

# Check dataset
str(data_act)
unique(data_act$bird_ID)
length(unique(data_act$bird_ID))
unique(data_act$ID)
length(unique(data_act$ID))


#---------------------------- CALCULATE DRY/WET PROP --------------------------


# Calculate proportion of time spent on water (daylight/night/dusk/dawn)
data_act$prop_day <- data_act$daynight.day.wet/
  (data_act$daynight.day.dry+data_act$daynight.day.wet)
data_act$prop_night <- data_act$daynight.night.wet/
  (data_act$daynight.night.dry+data_act$daynight.night.wet)
data_act$prop_dawn <- data_act$daynight.dusk.wet/
  (data_act$daynight.dusk.wet+data_act$daynight.dusk.dry)
data_act$prop_dusk <- data_act$daynight.dawn.wet/
  (data_act$daynight.dawn.wet+data_act$daynight.dawn.dry)

# Calculate proportion of different immersion activities
data_act$prop_forage <- data_act$act.foraging/24
hist(data_act$prop_forage)
data_act$prop_flightland <- data_act$act.flightland/24
hist(data_act$prop_flightland)
data_act$prop_onwater <- data_act$act.onwater/24
hist(data_act$prop_onwater)


#------------------------------------------------------------------------------
# MODELING
#------------------------------------------------------------------------------


#------------------------- Prop dry/wet per day -------------------------------


data_act_6 <- data_act[!is.na(data_act$prop_day),]

# Remove tracks with unknown breeding outcome
data_act_6 <- data_act_6[!data_act_6$BR=="Unknown",]

# Check the spread of the data
summary(data_act_6$prop_day) # 0s and 1s so need zero-one-inflated beta regression

# MODELS

# Pheno * sex * BR
mod_prop_day_inter <-
  brm(prop_day ~ sex * BR * phenophase + (1 | bird_ID) + (1 | year),
      data = data_act_6, family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_inter, file= "mod_prop_day_inter.RData")

# Pheno + sex + BR
mod_prop_day <-
  brm(prop_day ~ sex + BR + phenophase + (1 | bird_ID) + (1 | year),
      data = data_act_6, family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day, file= "mod_prop_day.RData")

# Pheno
mod_prop_day_pheno <-
  brm(prop_day ~ phenophase + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_pheno, file= "mod_prop_day_pheno.RData")

# Sex
mod_prop_day_sex <-
  brm(prop_day ~ sex + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_sex, file= "mod_prop_day_sex.RData")

# BR
mod_prop_day_BR <-
  brm(prop_day ~ BR + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_BR, file= "mod_prop_day_BR.RData")

# Pheno + BR
mod_prop_day_pheno_BR <-
  brm(prop_day ~ phenophase + BR + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_pheno_BR, file= "mod_prop_day_pheno_BR.RData")

# Pheno * BR
mod_prop_day_pheno_BR_inter <-
  brm(prop_day ~ phenophase *BR + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_pheno_BR_inter, file= "mod_prop_day_pheno_BR_inter.RData")

# Pheno + sex
mod_prop_day_pheno_sex <-
  brm(prop_day ~ phenophase + sex + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_pheno_sex, file= "mod_prop_day_pheno_sex.RData")

# Pheno * sex
mod_prop_day_pheno_sex_inter <-
  brm(prop_day ~ phenophase * sex + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_pheno_sex_inter, file= "mod_prop_day_pheno_sex_inter.RData")

# Sex + BR
mod_prop_day_sex_BR <-
  brm(prop_day ~ sex + BR + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_sex_BR, file= "mod_prop_day_sex_BR.RData")

#Sex * BR
mod_prop_day_sex_BR_inter <-
  brm(prop_day ~ sex * BR + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_sex_BR_inter, file= "mod_prop_day_sex_BR_inter.RData")

# NULL
mod_prop_day_null <-
  brm(prop_day ~ 1 + (1 | bird_ID) + (1 | year), data = data_act_6,
      family = zero_one_inflated_beta(), cores=2,
      control = list(adapt_delta = 0.999, max_treedepth=15))
save(mod_prop_day_null, file= "mod_prop_day_null.RData")

# Model comparison
mod1 <- loo(mod_prop_day)
mod2 <- loo(mod_prop_day_BR)
mod3 <- loo(mod_prop_day_null)
mod4 <- loo(mod_prop_day_pheno)
mod5 <- loo(mod_prop_day_pheno_sex_inter)
mod6 <- loo(mod_prop_day_sex)
mod7 <- loo(mod_prop_day_inter)
mod8 <- loo(mod_prop_day_pheno_sex)
mod9 <- loo(mod_prop_day_pheno_BR)
mod10 <- loo(mod_prop_day_pheno_BR_inter)
mod11 <- loo(mod_prop_day_sex_BR)
mod12 <- loo(mod_prop_day_sex_BR_inter)
comp_prop_forage <- loo_compare(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8,
                                mod9, mod10, mod11, mod12)
comp_prop_forage

# load best fit model
load(file = "mod_prop_day_pheno_sex_inter.RData")
mod_prop_day_pheno_sex_inter

# plot conditional effects (brms plotting function)
temp_conditions <- make_conditions(mod_prop_day_pheno_sex_inter,
                                   c("sex", "phenophase"))
plot(conditional_effects(mod_prop_day_pheno_sex_inter, effects = "phenophase:sex"),
     points= TRUE, point_args = list(width = 0.025, col="grey", alpha=0),
     errorbar_args = list(width = 0.1), theme = theme_classic(), ylab("test"),
     xlab("test"))

# Retrieve back transformed model estimates
temp <- conditional_effects(mod_prop_day_pheno_sex_inter)
temp$`phenophase:sex`

# Reformat model conditional effects for neater plotting
plot_data_day <- data.frame(matrix(nrow=6, ncol=5))
colnames(plot_data_day) <- c("phenophase", "Sex", "estimate","l_CI","u_CI")
plot_data_day$phenophase <- c("Back", "Back", "NB", "NB", "Out", "Out")
plot_data_day$Sex <- c("F", "M", "F", "M", "F", "M")
plot_data_day$estimate <- temp$`phenophase:sex`$estimate__*100
plot_data_day$l_CI <- temp$`phenophase:sex`$lower__*100
plot_data_day$u_CI <- temp$`phenophase:sex`$upper__*100
plot_data_day$phenophase <- factor(plot_data_day$phenophase,
                                   levels = c("Out", "NB", "Back"),
                                   ordered = TRUE )
plot_data_day$Sex <- as.factor(plot_data_day$Sex)

# plot of reformated conditional effects
ggplot(plot_data_day, aes(phenophase, estimate, col=Sex))+
  geom_point(size=4, position = position_dodge(width=0.2))+
  ylim(0,100)+
  geom_errorbar(aes(ymin=l_CI, ymax=u_CI), position = position_dodge(width=0.2),
                width=.1)+
  ylab(bquote("Time spent on water (%)\n"))+
  xlab("\nPhenophase")+
  theme_classic()+
  annotate("text",x="Back", y=100, label="Day", size =6)+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size=12))


###############################################################################


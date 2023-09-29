###############################################################################
#                          NB METRICS ANALYSES
###############################################################################

rm(list=ls())

library(fitdistrplus)
library(nlme)
library(AICcmodavg)
library(ggplot2)


# load data
GLS_summary_data <- read.csv(file="NB_metrics_all_tracks.4.csv", header = TRUE, 
                             sep = ",", na = "NA")

temp <- GLS_summary_data

temp$sex <- as.factor(temp$sex)
temp$BR <- as.factor(temp$BR)
temp$year <- as.factor(temp$year)


#---------------------------- DEPARTURE FROM CI -------------------------------


# Remove tracks with breeding success unknown
data_dep_doy <- temp[!temp$BR=="Unknown",]
# Check distribution of data
fit_norm <- fitdist(data_dep_doy$dep_DOY_adj, "norm")
plot(fit_norm)

# MODELS

# Sex + BR
mod_dep_doy_sex_BR <- lme(dep_DOY_adj ~ sex + BR, data = data_dep_doy, 
                          random = ~ 1 | year, method = "ML")

# Sex * BR
mod_dep_doy_sex_BR_inter <- lme(dep_DOY_adj ~ sex * BR, data = data_dep_doy, 
                                random = ~ 1 | year, method = "ML")

# Sex
mod_dep_doy_sex <- lme(dep_DOY_adj ~ sex, data = data_dep_doy, 
                       random = ~ 1 | year, method = "ML")

# BR
mod_dep_doy_BR <- lme(dep_DOY_adj ~ BR, data = data_dep_doy, 
                      random = ~ 1 | year, method = "ML")

# Null
mod_dep_doy_null <- lme(dep_DOY_adj ~ 0, data = data_dep_doy, 
                        random = ~ 1 | year, method = "ML")


# Compare models
models <- list(mod_dep_doy_BR,mod_dep_doy_null,mod_dep_doy_sex,mod_dep_doy_sex_BR,
               mod_dep_doy_sex_BR_inter)
model.names <- c("BR", "NULL", "Sex", "Sex + BR","Sex * BR")
aictab(models, model.names)


# Check best fit model (with restricted maximum likelihood)
mod_dep_doy_BR <- lme(dep_DOY_adj ~ BR, data = data_dep_doy, 
                      random = ~ 1 | year, method = "REML")

shapiro.test(resid(mod_dep_doy_BR))
qqnorm(resid(mod_dep_doy_BR))
qqline(resid(mod_dep_doy_BR))
plot(mod_dep_doy_BR)
summary(mod_dep_doy_BR)

# Plot results of best fit model
temp_plot <- as.data.frame(coef(summary(mod_dep_doy_BR)))
temp_plot$Value[2:nrow(temp_plot)] <- 
  (temp_plot$Value+temp_plot$Value[1])[2:nrow(temp_plot)]
temp_plot$`t-value`[2:nrow(temp_plot)] <- 
  (temp_plot$`t-value`+temp_plot$`t-value`[1])[2:nrow(temp_plot)]
temp_plot$l_CI <- temp_plot$Value - 1.96*temp_plot$Std.Error
temp_plot$u_CI <- temp_plot$Value + 1.96*temp_plot$Std.Error
temp_plot$cat1 <- c("Failed","Successful") 
temp_plot

ggplot(temp_plot, aes(cat1, Value))+
  geom_point(size=2.5, position = position_dodge(width=0.2))+
  geom_errorbar(aes(ymin=l_CI, ymax=u_CI), width=.05,
                position = position_dodge(width=0.2))+
  theme_classic()+
  ylab(bquote("Departure from CI (day of the year)"))+
  xlab("Breeding success")


# SUPPORTING R CODE FOR
# Evidence for regional-scale drought-induced declines in carabid beetles in old lowland beech forests
# Fabio Weiss, Susanne Winter, Dirk Pflugmacher, Thomas Kolling, Andreas Linde








# RECALCULATING DECLINES FOR WEISS ET AL. 2024
# Original article at: https://doi.org/10.1111/ecog.07020 
# Data at: https://doi.org/10.48548/pubdata-46
# Original R code at: https://github.com/fabioweiss/forest_carabid_declines






#### Package list #####

library(lme4)
library(bbmle)
library(sjPlot)
library(gridExtra)
library(DHARMa)
library(gamm4)
library(glmmTMB)
library(cols4all)
library(ggplot2)
library(ggeffects)
library(ggpubr)
library(psych)
library(dplyr)




#### Loading data ####

# loading raw data
# Data available at:  https://doi.org/10.48548/pubdata-46
beetle_samples <- read.csv("EWcarabids1999-2022_rawdata.csv")

# loading carabid weights
# carabid weights were calculated using the approach by Weiss & Linde (2022, https://doi.org/10.1007/s10841-022-00391-6)
carabid_weights <- read.csv("EWcarabids1999-2022_species_biomass.csv")

# loading sampling meta data can be recalculated using code from https://github.com/fabioweiss/forest_carabid_declines
sampling_meta <- read.csv2("sampling_meta_dwd_spei.csv")




#### prepare data  #####

beetle_samples$species[beetle_samples$species== "Lonicera pilicornis"] <- "Loricera pilicornis"


beetle_samples$species_weight <- carabid_weights$weight[match(beetle_samples$species, carabid_weights$species)]


beetle_samples$species_weight[beetle_samples$species =="no carabids"] <- 0

# assign abundance and biomass

# abundance
beetle_samples$sampling_abundance <- beetle_samples$abundance

# biomass
beetle_samples$weight_sum <- beetle_samples$species_weight * beetle_samples$sampling_abundance




#### Aggregate abundance and biomass ####

# aggregate by sampleID*trap
agg_samples <- aggregate(cbind(beetle_samples$sampling_abundance, beetle_samples$corr_abundance2, beetle_samples$weight_sum, beetle_samples$corr_weight_sum), by= list(beetle_samples$sample_id, beetle_samples$sample_id2, beetle_samples$year,  beetle_samples$interval, beetle_samples$trap, beetle_samples$plot, beetle_samples$site), FUN=sum  )

agg_samples<-agg_samples %>% 
  rename(
    sample_id = Group.1,
    sample_id2 = Group.2,
    year = Group.3,
    interval = Group.4,
    trap = Group.5,
    plot= Group.6,
    site=Group.7,
    sampling_abundance = V1,
    sampling_biomass = V2
  )


#### Add information about N.brevicollis####

# aggregate data of N.brevicollis
brevicollis <- beetle_samples[beetle_samples$species == "Nebria brevicollis",]

brevicollis_agg <- aggregate(cbind(brevicollis$sampling_abundance, brevicollis$corr_abundance2, brevicollis$weight_sum, brevicollis$corr_weight_sum), by= list(brevicollis$sample_id, brevicollis$sample_id2, brevicollis$year, brevicollis$interval,brevicollis$trap, brevicollis$plot, brevicollis$site), FUN=sum  )

brevicollis_agg<-brevicollis_agg %>% 
  rename(
    sample_id = Group.1,
    sample_id2 = Group.2,
    year = Group.3,
    interval = Group.4,
    trap = Group.5,
    plot= Group.6,
    site=Group.7,
    sampling_abundance = V1,
    sampling_biomass = V2
  )

# match n brevicollis abundance with full data
agg_samples$sampling_abundance_nb <- brevicollis_agg$sampling_abundance[match(agg_samples$sample_id2, brevicollis_agg$sample_id2)]
agg_samples$sampling_abundance_nb[is.na(agg_samples$sampling_abundance_nb)] <- 0

agg_samples$sampling_biomass_nb <- brevicollis_agg$sampling_biomass[match(agg_samples$sample_id2, brevicollis_agg$sample_id2)]
agg_samples$sampling_biomass_nb[is.na(agg_samples$sampling_biomass_nb)] <- 0



#### Match sampling metadata with aggregated samples ####

sampling_meta$sample_id <- as.factor(sampling_meta$sample_id)
agg_samples$sample_id <- as.factor(agg_samples$sample_id)

sampling_meta <- sampling_meta %>% 
  rename("sampling_length"="Trapping_length" )


# add metadata by merging
agg_samples2<- merge(agg_samples, sampling_meta[,c(7,9:30)], by.x = "sample_id", by.y ="sample_id") 

# remove all samples of intervals with more or less than 4 traps
agg_samples2 <- agg_samples2[agg_samples2$sampling_effort == 4,]


#### Restructuring data ####

#rename data
LTabundance <-agg_samples2




# year variable
LTabundance$year <- as.integer(LTabundance$year)
LTabundance$year2 <- as.factor(LTabundance$year)
LTabundance$scaled_year <- c(scale(LTabundance$year))

# interval
LTabundance$scaled_effort <- c(scale(LTabundance$sampling_length))

LTabundance$interval[LTabundance$interval== 1] <- "may"
LTabundance$interval[LTabundance$interval== 2] <- "june"
LTabundance$interval[LTabundance$interval== 3] <- "july"

LTabundance$interval <- as.factor(LTabundance$interval)
LTabundance$interval <- relevel(LTabundance$interval, ref = "may")

# abundances and biomass without Nebria brevicollis
LTabundance$sampling_abundance2 <- LTabundance$sampling_abundance - LTabundance$sampling_abundance_nb
LTabundance$sampling_biomass2 <- LTabundance$sampling_biomass - LTabundance$sampling_biomass_nb

# new trap ID for random intercept
LTabundance$trapID <- paste(LTabundance$plot,LTabundance$year,LTabundance$trap, sep="_")
LTabundance$trapID <- as.factor(LTabundance$trapID)
LTabundance$plot <- as.factor(LTabundance$plot)
LTabundance$site <-as.factor( LTabundance$site )

# scale precipitation variables

# scale temp by interval
intervals <- c(unique(LTabundance$interval))

LTabundance$scaled_sampling_temp <- NA
LTabundance$scaled_sampling_rain <- NA

for(i in intervals){
  LTabundance$scaled_sampling_temp[LTabundance$interval == i] <- c(scale(LTabundance$sampling_temp[LTabundance$interval == i]))
  LTabundance$scaled_sampling_rain[LTabundance$interval == i] <- c(scale(LTabundance$sampling_rain[LTabundance$interval == i]))
}


# select data from 1999-2001 and 2020-2022 

LTabundance2 <- LTabundance[LTabundance$year %in% c(1999:2001, 2020:2022),] 

LTabundance2$period <- "past"
LTabundance2$period[LTabundance2$year %in% c(2020:2022)] <- "present"



## ABUNDANCE #########

ab_glmm <- glmmTMB(sampling_abundance2 ~ period + scaled_sampling_temp + scaled_sampling_rain  + scaled_effort +I(scaled_effort^2) + (1|year2) + (1|interval)+ (1|site/plot/trapID), data= LTabundance2, family = nbinom2)

# predict mean abundance based on model
ggpreds1 <- ggpredict(ab_glmm, terms= c("period [all]"), ci.lvl = .95, type="fixed", condition = c(scaled_effort = 0, scaled_sampling_temp = 0,scaled_sampling_rain =0 ))

# relative change
(ggpreds1[1,]$predicted - ggpreds1[2,]$predicted) / ggpreds1[1,]$predicted

# bootstrap CIs of relative change 
# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model, terms= c("period [all]"), ci.lvl = NA, type="fixed", condition = c( scaled_effort =0, scaled_sampling_rain = 0,scaled_sampling_temp =0 ))
return(ggpreds$predicted)
}

# quantile approach for CIs
abundance_boot <- bootMer(ab_glmm, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed= 128)

boot_rates <- c()
abundance_boot$t

for (i in c(1:1000)){
  boot_rates[i] <-(abundance_boot$t[i,1] - abundance_boot$t[i,2]) / abundance_boot$t[i,1]
}

#lwr   
quantile(boot_rates, probs=.025, na.rm=TRUE)

#upr
quantile(boot_rates, probs=.975, na.rm=TRUE)




## BIOMASS #########

bio_glmm <- glmmTMB(sampling_biomass2 ~ period + scaled_sampling_temp + scaled_sampling_rain  + scaled_effort +I(scaled_effort^2) + (1|year2) + (1|interval)+ (1|site/plot/trapID), data= LTabundance2, ziformula = ~ period, family = ziGamma(link="log"))

# predict mean biomass based on model
ggpreds4 <- ggpredict(bio_glmm, terms= c("period [all]"),ci.lvl = .95, type="fe.zi", condition = c( scaled_effort =0, scaled_sampling_rain = 0,scaled_sampling_temp =0 ))

# relative change
(ggpreds4[1,]$predicted - ggpreds4[2,]$predicted) / ggpreds4[1,]$predicted

# bootstrap CIs of relative change 
# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model, terms= c("period [all]"), ci.lvl = NA, type="fe.zi", condition = c( scaled_effort =0, scaled_sampling_rain = 0,scaled_sampling_temp =0 ))
return(ggpreds$predicted)
}

# quantile approach for CIs
biomass_boot <- bootMer(bio_glmm, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed= 129)

boot_rates <- c()
biomass_boot$t

for (i in c(1:10)){
  boot_rates[i] <-(biomass_boot$t[i,1] - biomass_boot$t[i,2]) / biomass_boot$t[i,1]
}

#lwr   
quantile(boot_rates, probs=.025, na.rm=TRUE)

#upr
quantile(boot_rates, probs=.975, na.rm=TRUE)






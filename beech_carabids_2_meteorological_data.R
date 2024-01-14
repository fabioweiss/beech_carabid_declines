# SUPPORTING R CODE FOR
# Evidence for regional-scale drought-induced declines in carabid beetles in old lowland beech forests
# Fabio Weiss, Susanne Winter, Dirk Pflugmacher, Thomas Kolling, Andreas Linde




# INTERPOLATING METEOROLOGICAL DATA
# DWD data available at: https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/


# RSK = daily sum of rain
# TMK = daily mean temp
# TXK = daily max temp
# TNK = daily min temp

# DWD station IDs
# 164	angermuende
# 400	buch
# 1345	feldberg_1
# 7351	feldberg_2
# 6109	grambow
# 1869	gruenow
# 7389	heckelberg
# 3509	menz
# 3376	muencheberg
# 3552	neuruppin
# 5109	trollenhagen
# 5349	waren
# 5643	wittstock
# 5745	zehdenick

dwd_angerm <- read.csv("angermuende_dwd.txt", sep=";")
dwd_buch <- read.csv("buch_dwd.txt", sep=";")
dwd_gruenow <- read.csv("gruenow_dwd.txt", sep=";")
dwd_muenche <- read.csv("muencheberg_dwd.txt", sep=";")
dwd_waren <- read.csv("waren_dwd.txt", sep=";")
dwd_wittstock <- read.csv("wittstock_dwd.txt", sep=";")
dwd_zehdenick <- read.csv("zehdenick_dwd.txt", sep=";")
dwd_feldberg1 <- read.csv("feldberg1_dwd.txt", sep=";")
dwd_feldberg2 <- read.csv("feldberg2_dwd.txt", sep=";")
dwd_menz <- read.csv("menz_dwd.txt", sep=";")
dwd_grambow <- read.csv("grambow_dwd.txt", sep=";")
dwd_neuruppin <- read.csv("neuruppin_dwd.txt", sep=";")
dwd_trollen <- read.csv("trollenhagen_dwd.txt", sep=";")
dwd_heckelberg <- read.csv("heckelberg_dwd.txt", sep=";")


# dwd station locations
dwd_locations <- read.csv2("dwd_coordinates.csv")

# study site locations
study_sites <- read.csv("site_coordinates.csv")


### prepare data ####
dwd <- rbind(dwd_angerm, dwd_buch, dwd_gruenow, dwd_muenche, dwd_trollen, dwd_waren, dwd_wittstock, dwd_zehdenick, dwd_feldberg1, dwd_feldberg2, dwd_menz, dwd_grambow, dwd_neuruppin, dwd_heckelberg)

dwd$MESS_DATUM <- as.Date(ymd(dwd$MESS_DATUM ))

dwd$lat <- dwd_locations$lat[match(dwd$STATIONS_ID, dwd_locations$dwd_ID)]
dwd$long <- dwd_locations$long[match(dwd$STATIONS_ID, dwd_locations$dwd_ID)]


dwd_weather <- dwd[,c(1,2,7,14,20,21)]



dwd_weather2 <- dwd_weather[dwd_weather$MESS_DATUM >= "1974-01-01",]
dwd_weather2$date <- dwd_weather2$MESS_DATUM 
dwd_weather2$station <- dwd_locations$station[match(dwd_weather2$STATIONS_ID, dwd_locations$dwd_ID)]
dwd_weather2$RSK[dwd_weather2$RSK == -999] <- NA
dwd_weather2$TMK[dwd_weather2$TMK == -999] <- NA

dwd_temp <- dwd_weather2[,c(8,7,4)]
dwd_prec <- dwd_weather2[,c(8,7,3)]


# reformat data
model_temp <- reshape(dwd_temp, idvar = "date", timevar = "station", direction = "wide")
model_prec <- reshape(dwd_prec, idvar = "date", timevar = "station", direction = "wide")



### backcasting for all stations using a loop #####

dwd_weather_new <- dwd_weather
dwd_weather_new$type <- "measured"

model_sites <- c("angermuende", "buch", "gruenow", "waren", "muencheberg", "neuruppin", "zehdenick",
                 "neuruppin")

sites <- c( "trollenhagen", "heckelberg", "wittstock", "grambow", "menz", "feldberg_2" )

## loop #####

for(i in sites){
  
  model_sites[9] <- i
  
  # temp
  dwd_temp2 <- dwd_temp[dwd_temp$station %in% model_sites,]
  dwd_temp2$station[dwd_temp2$station == i] <- "modelstation"
  dwd_temp3 <- reshape(dwd_temp2, idvar = "date", timevar = "station", direction = "wide")
  dwd_temp3$month <- as.factor(c(substr(dwd_temp3$date,6,7)))
  model_temp <- dwd_temp3[complete.cases(dwd_temp3), ]
  
  m_temp_loop <- lm(TMK.modelstation ~ TMK.angermuende + TMK.buch + TMK.gruenow + TMK.waren  + TMK.muencheberg + TMK.neuruppin + TMK.zehdenick + month, data=model_temp)
  
  model_temp2 <- dwd_temp3[is.na(dwd_temp3$TMK.modelstation),]
  
  temp_preds <- predict.lm(m_temp_loop, newdata = model_temp2, se.fit = FALSE, type = "response")
  model_temp2$TMK.modelstation <- round(temp_preds,1)
  
  # precip
  dwd_prec2 <- dwd_prec[dwd_prec$station %in% model_sites,]
  dwd_prec2$station[dwd_prec2$station == i] <- "modelstation"
  dwd_prec3 <- reshape(dwd_prec2, idvar = "date", timevar = "station", direction = "wide")
  dwd_prec3$month <- as.factor(c(substr(dwd_prec3$date,6,7)))
  model_prec <- dwd_prec3[complete.cases(dwd_prec3), ]
  
  m_prec_loop <- lm(RSK.modelstation ~ RSK.angermuende + RSK.buch + RSK.gruenow + RSK.waren + RSK.muencheberg + RSK.neuruppin + RSK.zehdenick, data=model_prec)
  
  model_prec2 <- dwd_prec3[is.na(dwd_prec3$RSK.modelstation),]
  
  prec_preds <- predict.lm(m_prec_loop, newdata = model_prec2, se.fit = FALSE, type = "response")
  model_prec2$RSK.modelstation <- round(prec_preds,1)
  model_prec2$RSK.modelstation[model_prec2$RSK.modelstation <=0] <- 0 
  
  
  # combine
  loop_results <- data.frame(STATIONS_ID = dwd_locations$dwd_ID[dwd_locations$station == i],
                             MESS_DATUM = as.Date(union(model_temp2$date,model_prec2$date)))
  
  loop_results$RSK <- model_prec2$RSK.modelstation[match(loop_results$MESS_DATUM, model_prec2$date)]
  loop_results$TMK <- model_temp2$TMK.modelstation[match(loop_results$MESS_DATUM, model_temp2$date)]
  loop_results$lat <- dwd_locations$lat[dwd_locations$station == i]
  loop_results$long <- dwd_locations$long[dwd_locations$station == i]
  loop_results$type <- "predicted"
  
  dwd_weather_new <- rbind(dwd_weather_new, loop_results)
  model_sites <- model_sites[-9]
  
  print(i)
}


dwd_weather_modelled <- dwd_weather_new


## IDW interpolation loop for full data with gstat ####

dwd_weather_modelled$date_chr <- as.character(dwd_weather_modelled$MESS_DATUM)
dwd_weather_modelled$RSK[dwd_weather_modelled$RSK == -999] <- NA
dwd_weather_modelled$TMK[dwd_weather_modelled$TMK == -999] <- NA

test <- dwd_weather_modelled[dwd_weather_modelled$TMK == -999,]

# from 1974 until today
period <- as.character(seq.Date(
  as.Date(ymd("19740101")),
  as.Date(ymd("20221231")), by= 1))


dwd_interpolation <- data.frame(  date = c(period[1]), 
                                  site = c("dummy"), 
                                  temp = c(-999),
                                  prec = c(-999))




for(i in period){
  
loop_data <- dwd_weather_modelled[dwd_weather_modelled$date_chr== i,]
loop_data1 <- loop_data[!(is.na(loop_data$TMK)),]
loop_data2 <- loop_data[!(is.na(loop_data$RSK)),]


sample_temp <- data.frame(lat = loop_data1$lat,long= loop_data1$long,TMK = loop_data1$TMK)
coordinates(sample_temp) = ~long+lat
proj4string(sample_temp) <- CRS("+proj=longlat +datum=WGS84")

sample_rain <- data.frame(lat = loop_data2$lat,long= loop_data2$long,RSK = loop_data2$RSK)
coordinates(sample_rain) = ~long+lat
proj4string(sample_rain) <- CRS("+proj=longlat +datum=WGS84")

loc <- data.frame(lat = study_sites$lat,long=  study_sites$long)
coordinates(loc)  <- ~long+lat
proj4string(loc) <- CRS("+proj=longlat +datum=WGS84")

idw_temp <- idw(formula=TMK ~ 1, locations = sample_temp, newdata = loc, idp = 2)

idw_prec <- idw(formula=RSK ~ 1, locations = sample_rain, newdata = loc, idp = 2)

loop_inter <- data.frame(
  date = i, 
  site = study_sites$site_short, 
  temp = round(idw_temp@data$var1.pred, digits=1),
  prec = round(idw_prec@data$var1.pred, digits=1)
    )

dwd_interpolation <- rbind(dwd_interpolation, loop_inter)

print(i)
}

dwd_full <- dwd_interpolation






## Calculating SPEI for the period of study ####

## aggregate data by month and plot ####

dwd_full$year_month <- substr(dwd_full$date,1,7)

agg_dwd <- group_by(dwd_full, year_month, site) %>% summarise(temp = mean(temp), prec = sum(prec))



# calculate potential-evapotransperation sensu Thornthwaite (1948)
# lat = N 52.822110 (decimal degree)

sites <- c(study_sites$site_short)
agg_dwd2 <- agg_dwd

dwd_spei <- data.frame(  year_month = c("1900-01"), 
                         site = c("dummy"), 
                         temp = c(-999),
                         prec = c(-999),
                         PE = c(-999),
                         CWB = c(-999),
                         spei72 = c(-999) )

for(i in sites){
  
  loopdata <- agg_dwd2[agg_dwd2$site == i,]
  
  latitude <- study_sites$lat[study_sites$site_short == i]
  
  loopdata$PE <- thornthwaite(Tave = loopdata$temp, lat= latitude)
  
  loopdata$CWB <- loopdata$prec - loopdata$PE
  
  spei <- spei(loopdata$CWB, scale = 72)
  loopdata$spei72 <- spei$fitted
  
  dwd_spei <- rbind(dwd_spei, loopdata)
  
  print(i)
}

dwd_spei <- dwd_spei[-1,]



## creating Fig. S2 ####

# colors
colors <- c4a("berlin", 9)
colors[10] <- adjustcolor(colors[2], alpha=0.7)
colors[11] <- adjustcolor(colors[4], alpha=0.7)
colors[12] <- adjustcolor(colors[4], alpha=0)

dwd_spei$year <- as.numeric(substr(dwd_spei$year_month,1,4))
dwd_spei$month <- as.numeric(substr(dwd_spei$year_month,6,7))


dwd_spei$year_dec <- dwd_spei$year + 1/12 * dwd_spei$month


ggplot(dwd_spei, aes(year_dec, spei72, color=site)) + 
  theme_classic() +
  theme(text = element_text(size = 35)) +
  labs(y = "SPEI", x = "Year")+ 
  geom_vline(xintercept = 2000+4*1/12, lwd=20, color = colors[10])+
  geom_vline(xintercept = 1999+4*1/12, lwd=1, color = colors[2])+
  geom_vline(xintercept = 2001+4*1/12, lwd=1, color = colors[2])+
  
  geom_vline(xintercept = 2021+4*1/12, lwd=20, color = colors[10])+
  geom_vline(xintercept = 2020+4*1/12, lwd=1, color = colors[2])+
  geom_vline(xintercept = 2022+4*1/12, lwd=1, color = colors[2])+
  
  xlim(c(1967,2024))+
  geom_line(lwd=2)+ 
  scale_color_manual(values = c(c4a("reds2", 15)))+  
  theme(legend.position="none")




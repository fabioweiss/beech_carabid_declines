# SUPPORTING R CODE FOR
# Evidence for regional-scale drought-induced declines in carabid beetles in old lowland beech forests
# Fabio Weiss, Susanne Winter, Dirk Pflugmacher, Thomas Kolling, Andreas Linde




# PREPARING RAW DATA FOR ABUNDANCE/BIOMASS ANALYSES


## load data ####

setwd("C:/Users/fweiss/Promotion/PhD - own research/projects/Repeated Survey/data/clean data")
carabids <- read.csv("raw_carabids_publish.csv", sep=";")
carabids <- carabids[,-1]




## calculate biomass ####

# import species traits table
species_sizes <- read.csv2("species_sizes_beech.csv")
species_sizes$size_mid <- as.numeric((species_sizes$size_min + species_sizes$size_max) / 2)

# size-weight equations can be found here: https://doi.org/10.1007/s10841-022-00391-6

szyszko_weight<- function(size){
  ln_size <- log(size)
  ln_weight <- -8.92804283 + 2.5554921* ln_size
  weight_mg <- exp(ln_weight)*1000
  return(weight_mg)
}

booij_weight<- function(size){
  log_size <- log10(size)
  log_weight <- (-1.3 + 2.95 *log_size)
  weight_mg <- 10^(log_weight)
  return(weight_mg)
}

species_sizes$weight[species_sizes$size_mid >= 11.8] <- szyszko_weight(species_sizes$size_mid[species_sizes$size_mid >= 11.8])
species_sizes$weight[species_sizes$size_mid < 11.8] <- booij_weight(species_sizes$size_mid[species_sizes$size_mid < 11.8])

# assign weights to data based on species and calculate biomass
carabids$species_weight <- species_sizes$weight[match(carabids$species, species_sizes$species_new)]
carabids$species_weight[carabids$species =="no carabids"] <- 0
carabids$biomass <- carabids$abundance * carabids$species_weight



## aggregate by sample ####
carabids$sample_id <- paste(carabids$site_short, carabids$trap, carabids$year, carabids$interval, sep = "_")

agg_carabids <- group_by(carabids, sample_id, site_short, trap, year, interval, collection_date, quantitative_data, period) %>% summarise(abundance = sum(abundance), biomass = sum(biomass))



## add temperature and precipitation data #####

# see previous code file for data
dwd <- dwd_full

agg_carabids$collection_date <- as.Date(agg_carabids$collection_date)

agg_carabids$temp <- NA
agg_carabids$precip <- NA

rows <- row.names(agg_carabids)

for (i in rows){
  
  date <- agg_carabids[i,]$collection_date
  
  duration <- 14
  
  period <- as.character(seq.Date(date-duration , date , by= 1))
  
  site <- agg_carabids[i,]$site_short
  
  dwd1 <- dwd[dwd$date %in% period,]
  
  dwd2 <- dwd1[dwd1$site == site,]
  
  agg_carabids[i,]$temp <- mean(dwd2$temp)
  agg_carabids[i,]$precip <- sum(dwd2$prec)
  
 print(i) 
}




full_data <- agg_carabids

full_data$trap_ID <- paste(full_data$site_short, full_data$trap, sep = "_") 

# quick fix
full_data[full_data$site_short == "r3" & full_data$year == 2020,]
full_data$year[full_data$site_short == "r3" & full_data$year == 2020] <- 2021


## plot data availability #####
plot_data <- full_data

plot_data$plot_year <- 1901
plot_data$plot_year[plot_data$year == 2000] <- 1902
plot_data$plot_year[plot_data$year == 2001] <- 1903
plot_data$plot_year[plot_data$year == 2020] <- 1905
plot_data$plot_year[plot_data$year == 2021] <- 1906
plot_data$plot_year[plot_data$year == 2022] <- 1907


plot_data$plot_date <- ymd(paste(plot_data$plot_year, substr(plot_data$collection_date, 5,10), sep = ""))

plot_data$trap_ID <- paste(plot_data$site_short, plot_data$trap, sep = "_")
plot_data$trap_ID<- fct_relevel(plot_data$trap_ID)


plot_data$quantitative_data[plot_data$year== 1999 & plot_data$quantitative_data== 0] <- 2
plot_data$quantitative_data[plot_data$interval == 0] <-2
plot_data$quantitative_data[plot_data$year== 2022 & plot_data$interval == 8 & plot_data$quantitative_data== 0] <- 2 
plot_data$quantitative_data <- as.factor(plot_data$quantitative_data)


plot_data2 <- plot_data[plot_data$site_short != "w1",]


plot_data2$period[plot_data2$year %in% c(1999,2000,2001)] <- "past"
plot_data2$period[plot_data2$year %in% c(2020,2021, 2022)]  <- "present"

sum(plot_data2$abundance)
sum(plot_data2$abundance[plot_data2$period == "present"])
sum(plot_data2$abundance[plot_data2$period == "past"])


colors <- c4a("renoir", 12) 
colors[13]<- adjustcolor(colors[2], alpha= 0.9)



ggplot(plot_data , aes(plot_date, trap_ID, color=quantitative_data))+
  theme(text = element_text(size = 20)) +
  scale_color_manual(values = c( colors[2],colors[10],colors[8]), breaks = c("1","2","0"), labels= c("standardized", "non-standardized", "disturbed") )+
  labs(y = "", x = "")+
  geom_point(pch=15, size=4)+
  scale_x_continuous( labels= c( "1999","2000","2001","2020", "2021","2022"), breaks = as.Date(c("1901-07-01","1902-07-01","1903-07-01","1905-07-01","1906-07-01","1907-07-01" )) )+
  theme(legend.title=element_blank(),
        legend.position = "bottom")





## check intersection of the two periods and maintain only data from traps that were sampled in both periods ############

# exclude not suitable data
# due to disturbance and diverging sampling duration
full_data$quantitative_data[full_data$interval %in% c(222,333)] <-0
full_data3 <- full_data[full_data$quantitative_data == 1,]

full_data3$trap_ID <- paste(full_data3$site_short, full_data3$trap, sep = "_")

full_data3$period[full_data3$year %in% c(1999,2000,2001)] <- "past"
full_data3$period[full_data3$year %in% c(2020,2021, 2022)]  <- "present"



# split and check intersection
past <- full_data3[full_data3$year %in% c(1999,2000,2001) ,]
present <- full_data3[full_data3$year %in% c(2020,2021, 2022) ,]

ref <- c(as.character(unique(past$trap_ID)))

new <- c(as.character(unique(present$trap_ID)))

full <- intersect(ref,new)

setdiff(ref, full)

setdiff(new, full)



full_data4 <- full_data3[full_data3$trap_ID %in% full,]




## scale precipitation and temperature (separately per interval)  ####

intervals <- c(unique(full_data4$interval))

full_data4[full_data4$interval == 0,]

full_data4$scaled_temp2 <- NA
full_data4$scaled_precip2 <- NA

for(i in intervals){
  
  full_data4$scaled_temp2[full_data4$interval == i] <- c(scale(full_data4$temp[full_data4$interval == i]))
  full_data4$scaled_precip2[full_data4$interval == i] <- c(scale(full_data4$precip[full_data4$interval == i]))
  
}


# remove one study site with previous timber harvest
full_data5 <- full_data4[full_data4$site_short != "w1",]





























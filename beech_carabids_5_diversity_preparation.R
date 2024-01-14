# SUPPORTING R CODE FOR
# Evidence for regional-scale drought-induced declines in carabid beetles in old lowland beech forests
# Fabio Weiss, Susanne Winter, Dirk Pflugmacher, Thomas Kolling, Andreas Linde




# PREPARING RAW DATA FOR DIVERSITY AND TRAIT ANALYSES



## load raw data ####

setwd("C:/Users/fweiss/Promotion/PhD - own research/projects/Repeated Survey/data/clean data")
carabids <- read.csv("raw_carabids_publish.csv", sep=";")
full_data <- carabids[,-1]



# quick fix
full_data[full_data$site_short == "r3" & full_data$year == 2020,]
full_data$year[full_data$site_short == "r3" & full_data$year == 2020] <- 2021

# adjust quantitative data
full_data$quantitative_data[full_data$year== 1999 & full_data$quantitative_data== 0] <- 2
full_data$quantitative_data[full_data$interval == 0] <-2
full_data$quantitative_data[full_data$year== 2022 & full_data$interval == 8 & full_data$quantitative_data== 0] <- 2 



## exclude unsuitable data ####
# only maintain data from plots that were sampled for a whole season

# exclude grumsin in 2020 (incomplete season)
full_data3 <- full_data[!(full_data$site_short=="k3" & full_data$year == 2020),]

# exclude one trap in fauler ort (incomplete season)
full_data3 <- full_data3[!(full_data3$trap =="P4" & full_data3$site_short == "r3" & full_data3$year == 2021),]

# exclude 1999 except fauler ort (incomplete season)
full_data3 <- full_data3[!( !(full_data3$site_short=="r3" & full_data3$trap %in% c("P1","P2","P3","P4","P5" ) ) & full_data3$year == 1999),]

# exclude lüttenhagen (recent timber harvest)
full_data3 <- full_data3[full_data3$site_short != "w1",]

# exclude serrahn in 2020 (incomplete season)
full_data3 <- full_data3[!(full_data3$site_short=="r1" & full_data3$year == 2020),]

# exclude two traps in grumsin 2021 (incomplete season)
full_data3 <- full_data3[!(full_data3$site_short=="k3" & full_data3$year == 2021 & full_data3$trap %in% c("23","28")),]



## aggregate per plot & year ##########
full_data3$trap_ID <- paste(full_data3$site_short, full_data3$trap, sep = "_")
full_data3$trap_year <- paste(full_data3$trap_ID, full_data3$year, sep="_")

agg_sampling <- group_by(full_data3, trap_year, species, site_short, trap, year) %>% summarise(abundance = sum(abundance))
agg_sampling <- agg_sampling[agg_sampling$species != "no carabids",]

## how large are samples? ####
agg_sampling2 <- aggregate(cbind(agg_sampling$abundance), by= list( agg_sampling$trap_year), FUN=sum  )

# keep only samples with sample size >15 to increase potential extrapolation
keep <- agg_sampling2$Group.1[agg_sampling2$V1 >= 15]

full_data4 <- full_data3[full_data3$trap_year %in% keep,]



## check intersection of the 2 datasets ############

full_data4$period[full_data4$year %in% c(1999,2000,2001)] <- "past"
full_data4$period[full_data4$year %in% c(2020,2021, 2022)]  <- "present"

past <- full_data4[full_data4$year %in% c(1999,2000,2001) ,]
present <- full_data4[full_data4$year %in% c(2020,2021, 2022) ,]


ref <- c(as.character(unique(past$trap_ID)))

new <- c(as.character(unique(present$trap_ID)))

full <- intersect(ref,new)

setdiff(ref, full)

setdiff(new, full)

full_data4 <- full_data4[full_data4$trap_ID %in% full,]





# plot suitable data

plot_data <- full_data4

plot_data$plot_year <- 1901
plot_data$plot_year[plot_data$year == 2000] <- 1902
plot_data$plot_year[plot_data$year == 2001] <- 1903
plot_data$plot_year[plot_data$year == 2020] <- 1905
plot_data$plot_year[plot_data$year == 2021] <- 1906
plot_data$plot_year[plot_data$year == 2022] <- 1907


plot_data$plot_date <- ymd(paste(plot_data$plot_year, substr(plot_data$collection_date, 5,10), sep = ""))

plot_data$trap_ID <- paste(plot_data$site_short, plot_data$trap, sep = "_")
plot_data$trap_ID<- fct_relevel(plot_data$trap_ID)
plot_data$quantitative_data <- as.factor(plot_data$quantitative_data)


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




## transform to wide format ###########
agg_sampling <- group_by(full_data4, trap_year,  species, site_short, trap, year,trap_ID) %>% summarise(abundance = sum(abundance))
agg_sampling <- agg_sampling[agg_sampling$species != "no carabids",]

# wide format
wide_sampling <- data.frame(agg_sampling[-1, c(1,2,7)])
wide_sampling  <- reshape(wide_sampling, idvar = "species", timevar = "trap_year", direction = "wide")
wide_sampling[is.na(wide_sampling)] <- 0

row.names(wide_sampling) <- wide_sampling$species
wide_sampling <- wide_sampling[,-1] 



#### STANDARDIZED DIVERSITY ##########

stand_richness <- iNEXT(wide_sampling, q=c(0), datatype="abundance",knots = 30, endpoint=30)

stand_richness <- data.frame(stand_richness$iNextEst$coverage_based)


# find max coverage 
ids <- c(unique(stand_richness$Assemblage))
max_SC <- data.frame(assemblage = stand_richness$Assemblage , max_SC=NA)
 
  for(i in ids){
    loopdata <- stand_richness[stand_richness$Assemblage == i,]
    max_SC$max_SC[max_SC$assemblage == i] <- max(loopdata$SC)
  }
min(max_SC$max_SC)  


# -> 0.719

# repeat for coverage 0.72
diversity <- group_by(full_data4, trap_year, site_short, trap, year) %>% summarise(abundance = sum(abundance))

richness72 <- estimateD(wide_sampling, q = c(0), datatype = "abundance", base="coverage", level=0.72, conf=0.95)
richness72$Assemblage <- substring(richness72$Assemblage, 11)

simpson72 <- estimateD(wide_sampling, q = c(2), datatype = "abundance", base="coverage", level=0.72, conf=0.95)
simpson72$Assemblage <- substring(simpson72$Assemblage, 11)

diversity$stand_richness <- richness72$qD[match(diversity$trap_year, richness72$Assemblage)]
diversity$stand_simpson <- simpson72$qD[match(diversity$trap_year, simpson72$Assemblage)]
diversity$stand_evenness <- diversity$stand_simpson / diversity$stand_richness

# for later use
stand_diversity <- diversity



#### TRAITS ###########

### prepare trait data #####

# load carabids.org data complement based on Müller-Motzfeld 2004. (FREUDE HARDE LOSE)
traits <- read.csv2("C:/Users/fweiss/Promotion/PhD - own research/projects/Repeated Survey/data/clean data/traits_modified.csv")

# convert to numerical trait variables

traits$predator <- 0
traits$predator[traits$guild == "Predator"] <- 1
traits$predator[is.na(traits$guild)] <- NA

traits$shortwinged <- 0
traits$shortwinged[traits$wings == "shortwinged"] <- 1
traits$shortwinged[is.na(traits$wings)] <- NA

traits$dimorphic <- 0
traits$dimorphic[traits$wings == "dimorphic"] <- 1
traits$dimorphic[is.na(traits$wings)] <- NA

traits$winged <- 0
traits$winged[traits$wings == "winged"] <- 1
traits$winged[is.na(traits$wings)] <- NA

traits$short_di <- 0
traits$short_di[traits$wings %in% c("shortwinged", "dimorphic")] <- 1
traits$short_di[is.na(traits$wings)] <- NA



### Community means #####

# use data as for standardized diversity (aggregated per plot and year)
# continue with full_data3 (before exclusion of small samples <15)

## check intersection of the 2 datasets
full_data3$trap_ID <- paste(full_data3$site_short, full_data3$trap, sep = "_")
full_data3$trap_year <- paste(full_data3$trap_ID, full_data3$year, sep="_")

full_data3$period[full_data3$year %in% c(1999,2000,2001)] <- "past"
full_data3$period[full_data3$year %in% c(2020,2021, 2022)]  <- "present"

past <- full_data3[full_data3$year %in% c(1999,2000,2001) ,]
present <- full_data3[full_data3$year %in% c(2020,2021, 2022) ,]


ref <- c(as.character(unique(past$trap_ID)))

new <- c(as.character(unique(present$trap_ID)))

full <- intersect(ref,new)

setdiff(ref, full)

setdiff(new, full)

full_data5 <- full_data3[full_data3$trap_ID %in% full,]


# aggregate per year and plot 
agg_sampling <- group_by(full_data5, trap_year,  species, site_short, trap, year,trap_ID) %>% summarise(abundance = sum(abundance))
agg_sampling <- agg_sampling[agg_sampling$species != "no carabids",]


# merge traits with data
agg_sampling$size <- traits$size[match(agg_sampling$species, traits$species)]
agg_sampling$predator <- traits$predator[match(agg_sampling$species, traits$species)]
agg_sampling$shortwinged <- traits$shortwinged[match(agg_sampling$species, traits$species)]
agg_sampling$dimorphic <- traits$dimorphic[match(agg_sampling$species, traits$species)]
agg_sampling$winged <- traits$winged[match(agg_sampling$species, traits$species)]
agg_sampling$short_di <- traits$short_di[match(agg_sampling$species, traits$species)]


# calculate CMs

# use a manual approach with dplyr (https://rpubs.com/CPEL/cwm)

cm_data <-   # New dataframe where we can inspect the result
  agg_sampling %>%   # First step in the next string of statements
  group_by(trap_year) %>%   # Groups the summary file by Plot number
  summarize(           # Coding for how we want our CWMs summarized
    size_cm = mean(na.omit(size))
  )


## add variables
cm_data$site_short <- agg_sampling$site_short[match(cm_data$trap_year, agg_sampling$trap_year)]
cm_data$trap <- agg_sampling$trap[match(cm_data$trap_year, agg_sampling$trap_year)]
cm_data$year <- agg_sampling$year[match(cm_data$trap_year, agg_sampling$trap_year)]
cm_data$trap_ID <- agg_sampling$trap_ID[match(cm_data$trap_year, agg_sampling$trap_year)]

cm_data$period <- NA
cm_data$period[cm_data$year %in% c(1999:2001)] <- "past"
cm_data$period[cm_data$year %in% c(2020:2022)] <- "present"



## save binomial data

agg_sampling$period <- NA
agg_sampling$period[agg_sampling$year %in% c(1999:2001)] <- "past"
agg_sampling$period[agg_sampling$year %in% c(2020:2022)] <- "present"

# for later use 
BIN_species <- agg_sampling




### Community WEIGHTED means ##########

# adjust data selection
# use all data and model per interval

# reload raw data
full_data <- carabids[,-1]

# quick fix
full_data[full_data$site_short == "r3" & full_data$year == 2020,]
full_data$year[full_data$site_short == "r3" & full_data$year == 2020] <- 2021

# adjust quantitative data
full_data$quantitative_data[full_data$year== 1999 & full_data$quantitative_data== 0] <- 2
full_data$quantitative_data[full_data$interval == 0] <-2
full_data$quantitative_data[full_data$year== 2022 & full_data$interval == 8 & full_data$quantitative_data== 0] <- 2 

# adjust intervals to make them comparable
full_data$interval[full_data$interval == 222] <- 8
full_data$interval[full_data$interval == 333] <- 9

full_data$interval[full_data$trap %in% c("P1","P2","P3","P4", "P5") & full_data$year == 1999 & full_data$interval == 7] <- 15
full_data$interval[full_data$trap %in% c("P1","P2","P3","P4", "P5") & full_data$year == 1999 & full_data$interval == 6] <- 13
full_data$interval[full_data$trap %in% c("P1","P2","P3","P4", "P5") & full_data$year == 1999 & full_data$interval == 5] <- 11
full_data$interval[full_data$trap %in% c("P1","P2","P3","P4", "P5") & full_data$year == 1999 & full_data$interval == 4] <- 9
full_data$interval[full_data$trap %in% c("P1","P2","P3","P4", "P5") & full_data$year == 1999 & full_data$interval == 3] <- 7
full_data$interval[full_data$trap %in% c("P1","P2","P3","P4", "P5") & full_data$year == 1999 & full_data$interval == 2] <- 5
full_data$interval[full_data$trap %in% c("P1","P2","P3","P4", "P5") & full_data$year == 1999 & full_data$interval == 1] <- 3


# exclude not suitable data
full_data <- full_data[full_data$quantitative_data %in% c(1,2), ]
full_data <- full_data[full_data$site_short != "w1",]



# check intersection
full_data$trap_ID <- paste(full_data$site_short, full_data$trap, sep = "_")

full_data$period[full_data$year %in% c(1999,2000,2001)] <- "past"
full_data$period[full_data$year %in% c(2020,2021, 2022)]  <- "present"

past <- full_data[full_data$year %in% c(1999,2000,2001) ,]
present <- full_data[full_data$year %in% c(2020,2021, 2022) ,]

ref <- c(as.character(unique(past$trap_ID)))

new <- c(as.character(unique(present$trap_ID)))

full <- intersect(ref,new)

setdiff(ref, full)

setdiff(new, full)

full_data6 <- full_data[full_data$trap_ID %in% full,]



## aggregating samples
full_data6$sample_ID <- paste( full_data6$site_short, full_data6$trap, full_data6$year, full_data6$interval, sep = "_")

agg_sampling <- group_by(full_data6, sample_ID, site_short, trap, interval, year,  species,) %>% summarise(abundance = sum(abundance))
agg_sampling <- agg_sampling[agg_sampling$species != "no carabids",]




# merge traits with data
agg_sampling$size <- traits$size[match(agg_sampling$species, traits$species)]
agg_sampling$predator <- traits$predator[match(agg_sampling$species, traits$species)]
agg_sampling$shortwinged <- traits$shortwinged[match(agg_sampling$species, traits$species)]
agg_sampling$dimorphic <- traits$dimorphic[match(agg_sampling$species, traits$species)]
agg_sampling$winged <- traits$winged[match(agg_sampling$species, traits$species)]
agg_sampling$short_di <- traits$short_di[match(agg_sampling$species, traits$species)]


## calculate CWM 

# use a manual approach with dplyr (https://rpubs.com/CPEL/cwm)

cwm_data <-   # New dataframe where we can inspect the result
  agg_sampling %>%   # First step in the next string of statements
  group_by(sample_ID) %>%   # Groups the summary file by Plot number
  summarize(           # Coding for how we want our CWMs summarized
    size_cwm = weighted.mean(na.omit(size), abundance)
  )


## add variables
cwm_data$site_short <- agg_sampling$site_short[match(cwm_data$sample_ID, agg_sampling$sample_ID)]
cwm_data$trap <- agg_sampling$trap[match(cwm_data$sample_ID, agg_sampling$sample_ID)]
cwm_data$year <- agg_sampling$year[match(cwm_data$sample_ID, agg_sampling$sample_ID)]
cwm_data$interval <- agg_sampling$interval[match(cwm_data$sample_ID, agg_sampling$sample_ID)]

cwm_data$trap_ID <- paste(cwm_data$site_short, cwm_data$trap, sep="_")   

cwm_data$period <- NA
cwm_data$period[cwm_data$year %in% c(1999:2001)] <- "past"
cwm_data$period[cwm_data$year %in% c(2020:2022)] <- "present"


## save binomial data

agg_sampling$period <- NA
agg_sampling$period[agg_sampling$year %in% c(1999:2001)] <- "past"
agg_sampling$period[agg_sampling$year %in% c(2020:2022)] <- "present"


agg_sampling2 <- agg_sampling %>%
  mutate(freq =abundance) %>%
  uncount(freq)
 
# for later use 
BIN_indiv <- agg_sampling2


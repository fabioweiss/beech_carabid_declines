
# SUPPORTING R CODE FOR
# Evidence for regional-scale drought-induced declines in carabid beetles in old lowland beech forests
# Fabio Weiss, Susanne Winter, Dirk Pflugmacher, Thomas Kolling, Andreas Linde




# DIVERSITY AND TRAIT ANALYSES
# for data see previous code file





# use data products from previous code files


#setwd("C:/Users/fweiss/Promotion/PhD - own research/projects/Repeated Survey/data/clean data")
#stand_diversity <- read.csv2("standardized_diversity.csv")
#site_trends1 <- read.csv2("site_trends1.csv")


stand_diversity

cwm_data 
cm_data 

BIN_indiv
BIN_species

site_trends1 <- site_trends


## colors
colors <- c4a("renoir", 12)
colors[13] <- adjustcolor(colors[11], alpha=0.9)
colors[14] <- adjustcolor(colors[8], alpha=0.9)
colors[15] <- adjustcolor(colors[2], alpha=0.9)

show_col(colors)



## prepare data ############
stand_diversity$period[stand_diversity$year %in% c("1999","2000","2001")] <- "past"
stand_diversity$period[stand_diversity$year %in% c("2020","2021", "2022")]  <- "present"

stand_diversity$trap_ID <- paste(stand_diversity$site_short, stand_diversity$trap, sep="_")




## RICHNESS ######
richness1 <- glmmTMB(stand_richness ~ period + (1|year) +(1|site_short/trap_ID), data= stand_diversity, family = gaussian(link = "identity"))

tab_model(richness1)

dharma_sim1 <- simulateResiduals(fittedModel = richness1, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)


# site specific model
richness2 <- glmmTMB(stand_richness ~ period*site_short + (1|year)+ (1|trap_ID), data= stand_diversity, family = gaussian(link = "identity"))

dharma_sim1 <- simulateResiduals(fittedModel = richness2, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)


# plotting
ggpreds1 <- ggpredict(richness1, terms= c("period [all]"),ci.lvl = .95, type="fixed")
ggpreds1$plotting_year <- c(2,4)

ggpreds2 <- ggpredict(richness2, terms= c("period [all]", "site_short [all]"),ci.lvl = NA, type="fixed")
ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4

plot1 <- ggplot(ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  ylim(c(2,8))+
  labs(y = "Richness (q=0)", x = "")+
  scale_y_continuous( limits = c(2,8), labels= c("2","4","6","     8"), breaks= c(2,4,6,8))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=12, color= colors[11])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[11])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "n.s.", y_position = 7),
                  textsize = 8, vjust = -0.2,size = 1,
                  manual = F) 

plot2 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(2,8), labels= c("2","4","6","     8"), breaks= c(2,4,6,8))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")



## extract site specific trends ###
past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 
site_trends1$richness <- (present$predicted - past$predicted) 



## INVERSE SIMPSON ##########
simpson1 <- glmmTMB(stand_simpson ~ period + (1|year) +(1|site_short/trap_ID), data= stand_diversity, family = gaussian(link = "identity"))

tab_model(simpson1)

dharma_sim1 <- simulateResiduals(fittedModel = simpson1, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)


# site specific model
simpson2 <- glmmTMB(stand_simpson ~ period*site_short + (1|year)+ (1|trap_ID), data= stand_diversity, family = gaussian(link = "identity"))

dharma_sim1 <- simulateResiduals(fittedModel = simpson2, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)


# plotting

ggpreds1 <- ggpredict(simpson1, terms= c("period [all]"),ci.lvl = .95, type="fixed")
ggpreds1$plotting_year <- c(2,4)

ggpreds2 <- ggpredict(simpson2, terms= c("period [all]", "site_short [all]"),ci.lvl = NA, type="fixed")
ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4


plot3 <- ggplot(ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "Simpson (q=2)", x = "")+
  scale_y_continuous( limits = c(1.5,4.5), labels= c("2","3","     4"), breaks= c(2,3,4))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=12, color= colors[11])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[11])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "n.s.", y_position = 4),
                  textsize = 8, vjust = -0.2,size = 1,
                  manual = F) 

plot4 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(1.5,4.55), labels= c("2","3","     4"), breaks= c(2,3,4))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")


## extract site specific trends ###

past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 
site_trends1$simpson <- (present$predicted - past$predicted) 




## EVENNESS ##########

evenness1 <- glmmTMB(stand_evenness ~ period + (1|year) +(1|site_short/trap_ID), data= stand_diversity, family = beta_family(link = "logit"))

dharma_sim1 <- simulateResiduals(fittedModel = evenness1, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)

tab_model(evenness1)

# site specific model
evenness2 <- glmmTMB(stand_evenness ~ period*site_short + (1|year)+ (1|trap_ID), data= stand_diversity, family = beta_family(link = "logit"))

dharma_sim1 <- simulateResiduals(fittedModel = evenness2, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)


# plotting

ggpreds1 <- ggpredict(evenness1, terms= c("period [all]"),ci.lvl = .95, type="fixed")
ggpreds1$plotting_year <- c(2,4)

ggpreds2 <- ggpredict(evenness2, terms= c("period [all]", "site_short [all]"),ci.lvl = NA, type="fixed")
ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4



plot5 <- ggplot(ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "Evenness", x = "")+
  scale_y_continuous( limits = c(0.5,0.9), labels= c("0.5","0.7","  0.9"), breaks= c(0.5,0.7,0.9))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=12, color= colors[11])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[11])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "n.s.", y_position = 0.85),
                  textsize = 8, vjust = -0.2,size = 1,
                  manual = F) 

plot6 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(0.5,0.9), labels= c("0.5","0.7","  0.9"), breaks= c(0.5,0.7,0.9))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")


## extract site specific trends ###

past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 
site_trends1$evenness <- (present$predicted - past$predicted) 




# plot
# Fig.3
grid.arrange(plot1, plot2,
             plot3, plot4,
             plot5, plot6, 
             widths= c(1,1), ncol=2)
#export 10x15 portrait



### TRAITS ######


### SIZE ############
cwm_data$trap_ID <- paste(cwm_data$site_short, cwm_data$trap, sep="_")

## CWM size  ########

cwm_size <- glmmTMB(size_cwm ~ period + (1|year) + (1|interval) + (1|site_short/trap_ID), data= cwm_data, family = gaussian(link = "identity"))

dharma_sim1 <- simulateResiduals(fittedModel = cwm_size, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)

tab_model(cwm_size)


# site specific model
cwm_size2 <- glmmTMB(size_cwm ~ period *site_short + (1|year) + (1|interval) + (1|trap_ID), data= cwm_data, family = gaussian(link = "identity"))

dharma_sim1 <- simulateResiduals(fittedModel = cwm_size2, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)



# plotting
ggpreds1 <- ggpredict(cwm_size, terms= c("period [all]"),ci.lvl = .95, type="fixed")
ggpreds1$plotting_year <- c(2,4)

ggpreds2 <- ggpredict(cwm_size2, terms= c("period [all]", "site_short [all]"),ci.lvl = NA, type="fixed")
ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4




plot1 <- ggplot(ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "Size (CWM)", x = "")+
  scale_y_continuous( limits = c(13,20), labels= c("14","16", "18", "     20"), breaks= c(14,16,18,20))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=12, color= colors[2])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[2])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "*", y_position = 19.5),
                  textsize = 8, vjust = -0.2,size = 1,
                  manual = F) 

plot2 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(13,20), labels= c("14","16", "18", "     20"), breaks= c(14,16,18,20))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")




## extract site specific trends ###

past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 
site_trends1$cwm_size <- (present$predicted - past$predicted) 

site_trends1$cwm_original <- past$predicted


## CM size  ########
cm_size <- glmmTMB(size_cm ~ period + (1|year) + (1|site_short/trap_ID), data= cm_data, family = gaussian(link = "identity"))

dharma_sim1 <- simulateResiduals(fittedModel = cm_size, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)

tab_model(cm_size)

# site specific model
cm_size2 <- glmmTMB(size_cm ~ period *site_short + (1|year) + (1|trap_ID), data= cm_data, family = gaussian(link = "identity"))

dharma_sim1 <- simulateResiduals(fittedModel = cm_size2, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)


# plotting

ggpreds1 <- ggpredict(cm_size, terms= c("period [all]"),ci.lvl = .95, type="fixed")
ggpreds1$plotting_year <- c(2,4)

ggpreds2 <- ggpredict(cm_size2, terms= c("period [all]", "site_short [all]"),ci.lvl = NA, type="fixed")
ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4





plot3 <- ggplot(ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "Size (CM)", x = "")+
  scale_y_continuous( limits = c(12,22), labels= c("14","16", "18", "  20"), breaks= c(14,16,18,20))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=12, color= colors[8])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[8])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "n.s.", y_position = 19),
                  textsize = 8, vjust = -0.2,size = 1,
                  manual = F) 

plot4 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(12,22), labels= c("14","16", "18", "  20"), breaks= c(14,16,18,20))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")



## extract site specific trends ###

past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 
site_trends1$cm_size <- (present$predicted - past$predicted) 


### WINGS ####

BIN_species$trap_ID <- paste(BIN_species$site_short, BIN_species$trap, sep="_")

BIN_indiv$trap_ID <- paste(BIN_indiv$site_short, BIN_indiv$trap, sep="_")



## individual based 

ind_wings1 <- glmmTMB(winged ~ period + (1|year) + (1|interval) + (1|site_short/trap_ID), data= BIN_indiv, family = "binomial")

dharma_sim1 <- simulateResiduals(fittedModel = ind_wings1, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)

tab_model(ind_wings1)


# site specific model
ind_wings2 <- glmmTMB(winged ~ period *site_short + (1|year) + (1|interval) + (1|trap_ID), data= BIN_indiv, family = "binomial")

dharma_sim1 <- simulateResiduals(fittedModel = ind_wings2, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)


# plotting

ggpreds1 <- ggpredict(ind_wings1, terms= c("period [all]"),ci.lvl = .95, type="fixed")
ggpreds1$plotting_year <- c(2,4)

ggpreds2 <- ggpredict(ind_wings2, terms= c("period [all]", "site_short [all]"),ci.lvl = NA, type="fixed")
ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4


plot5 <- ggplot(ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "P of flying individual", x = "")+
  scale_y_continuous( limits = c(0,0.2), labels= c("0.05","0.10", "0.15", "  0.20"), breaks= c(0.05,0.1,0.15,0.2))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=12, color= colors[2])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[2])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "*", y_position = 0.15),
                  textsize = 8, vjust = -0.2,size = 1,
                  manual = F) 

plot6 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(0,0.2), labels= c("0.05","0.10", "0.15", "  0.20"), breaks= c(0.05,0.1,0.15,0.2))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")



## extract site specific trends ###
past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 
site_trends1$ind_wings <- present$predicted - past$predicted
 


## species based


spec_wings1 <- glmmTMB(winged ~ period + (1|year) + (1|site_short/trap_ID), data= BIN_species, family = "binomial")

dharma_sim1 <- simulateResiduals(fittedModel = spec_wings1, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)

tab_model(spec_wings1)


# site specific model
spec_wings2 <- glmmTMB(winged ~ period *site_short + (1|year)  + (1|trap_ID), data= BIN_species, family = "binomial")

dharma_sim1 <- simulateResiduals(fittedModel = spec_wings2, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)


# plotting

ggpreds1 <- ggpredict(spec_wings1, terms= c("period [all]"),ci.lvl = .95, type="fixed")
ggpreds1$plotting_year <- c(2,4)

ggpreds2 <- ggpredict(spec_wings2, terms= c("period [all]", "site_short [all]"),ci.lvl = NA, type="fixed")
ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4



plot7 <- ggplot(ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "P of flying species", x = "")+
  scale_y_continuous( limits = c(0,0.4), labels= c("0.1","0.2", "0.3", "  0.4"), breaks= c(0.1,0.2,0.3,0.4))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=12, color= colors[8])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[8])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "n.s.", y_position = 0.3),
                  textsize = 8, vjust = -0.2,size = 1,
                  manual = F) 

plot8 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(0,0.4), labels= c("0.1","0.2", "0.3", "  0.4"), breaks= c(0.1,0.2,0.3,0.4))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")



## extract site specific trends ###
past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 
site_trends1$spec_wings <- present$predicted - past$predicted


### PREDATOR ######


## individual based 

ind_pred1 <- glmmTMB(predator ~ period + (1|year) + (1|interval) + (1|site_short/trap_ID), data= BIN_indiv, family = "binomial")

dharma_sim1 <- simulateResiduals(fittedModel = ind_pred1, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)

tab_model(ind_pred1)


# site specific model
ind_pred2 <- glmmTMB(predator ~ period *site_short + (1|year) + (1|interval) + (1|trap_ID), data= BIN_indiv, family = "binomial")

dharma_sim1 <- simulateResiduals(fittedModel = ind_pred2, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)



# plotting

ggpreds1 <- ggpredict(ind_pred1, terms= c("period [all]"),ci.lvl = .95, type="fixed")
ggpreds1$plotting_year <- c(2,4)

ggpreds2 <- ggpredict(ind_pred2, terms= c("period [all]", "site_short [all]"),ci.lvl = NA, type="fixed")
ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4


plot9 <- ggplot(ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "P of predatory individual", x = "")+
  scale_y_continuous( limits = c(0.955,1.005), labels= c("0.96","0.97","0.98","0.99", "  1.00"), breaks= c(0.96,0.97,0.98,0.99,1))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=12, color= colors[2])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[2])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "n.s.", y_position = 1.004),
                  textsize = 8, vjust = -0.2,size = 1,
                  manual = F) 

plot10 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(0.955,1.005), labels= c("0.96","0.97","0.98","0.99", "  1.00"), breaks= c(0.96,0.97,0.98,0.99,1))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")



## extract site specific trends ###
past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 
site_trends1$ind_pred <- present$predicted - past$predicted


## species based
  
  
spec_pred1 <- glmmTMB(predator ~ period + (1|year) + (1|site_short/trap_ID), data= BIN_species, family = "binomial")
  
  dharma_sim1 <- simulateResiduals(fittedModel = spec_pred1, re.form= NULL)
  plot(dharma_sim1)
  testDispersion(dharma_sim1)
  
  tab_model(spec_pred1)
  
  
# site specific model
spec_pred2 <- glmmTMB(predator ~ period *site_short + (1|year)  + (1|trap_ID), data= BIN_species, family = "binomial")
  
  dharma_sim1 <- simulateResiduals(fittedModel = spec_pred2, re.form= NULL)
  plot(dharma_sim1)
  testDispersion(dharma_sim1)
  

  
# plotting
  
ggpreds1 <- ggpredict(spec_pred1, terms= c("period [all]"),ci.lvl = .95, type="fixed")
ggpreds1$plotting_year <- c(2,4)
  
ggpreds2 <- ggpredict(spec_pred2, terms= c("period [all]", "site_short [all]"),ci.lvl = NA, type="fixed")
ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4
  

plot11 <- ggplot(ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "P of predatory species", x = "")+
  scale_y_continuous( limits = c(0.8,1), labels= c("0.85","0.90", "0.95", " 1.00"), breaks= c(0.85,0.9,0.95,1))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=12, color= colors[8])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[8])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "n.s.", y_position = 0.99),
                  textsize = 8, vjust = -0.2,size = 1,
                  manual = F) 

plot12 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(0.8,1), labels= c("0.85","0.90", "0.95", " 1.00"), breaks= c(0.85,0.9,0.95,1))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")


## extract site specific trends ###
past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 
site_trends1$spec_pred <- present$predicted - past$predicted



### full plot #####
# Fig. 4
#export 21x28

grid.arrange(plot1, plot2, plot3, plot4,
             plot5, plot6, plot7, plot8,
             plot9, plot10, plot11, plot12,
             widths= c(1,1,1,1), ncol=4)










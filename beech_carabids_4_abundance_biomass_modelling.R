
# SUPPORTING R CODE FOR
# Evidence for regional-scale drought-induced declines in carabid beetles in old lowland beech forests
# Fabio Weiss, Susanne Winter, Dirk Pflugmacher, Thomas Kolling, Andreas Linde




# ABUNDANCE/BIOMASS ANALYSES
# for data see previous code file



# colors
# multi-color palette

colors <- c4a("renoir", 12) 
colors[13]<- adjustcolor(colors[2], alpha= 0.9)


## prepare data ####
full_data5$scaled_interval <- c(scale(as.numeric(full_data5$interval)))
full_data5$year2 <- as.factor(factor(full_data5$year, ordered = F))
full_data5$interval1 <- as.factor(full_data5$interval)







## ABUNDANCE ####

m1 <- glmmTMB(abundance ~ period + scaled_precip2 + scaled_temp2 + (1|year2) + (1|interval1) + (1|site_short/trap_ID) , data= full_data5, ziformula = ~ 0, family = nbinom2)

tab_model(m1)



# DHARMa diagnostics
dharma_sim1 <- simulateResiduals(fittedModel = m1, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)
testZeroInflation(dharma_sim1)



# predict mean abundance based on model
ggpreds1 <- ggpredict(m1, terms= c("period [all]"), ci.lvl = .95, type="fixed", condition = c(scaled_precip2 = 0,scaled_temp2 =0 ))

ggpreds1$plotting_year[ggpreds1$x == "past"] <- 2
ggpreds1$plotting_year[ggpreds1$x == "present"] <- 4




## calculate relative change in abundance ##########

# relative change
(ggpreds1[1,]$predicted - ggpreds1[2,]$predicted) / ggpreds1[1,]$predicted

# bootstrap CIs of relative change 

# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model, terms= c("period [all]"), ci.lvl = .95, type="fixed", condition = c(scaled_precip2 = 0,scaled_temp2 =0 ))
return(ggpreds$predicted)
}

# quantile approach for CIs
abundance_boot <- bootMer(m1, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed= 128)

boot_rates <- c()
abundance_boot$t

for (i in c(1:1000)){
  boot_rates[i] <-(abundance_boot$t[i,1] - abundance_boot$t[i,2]) / abundance_boot$t[i,1]
}

#lwr   
quantile(boot_rates, probs=.025, na.rm=TRUE)

#upr
quantile(boot_rates, probs=.975, na.rm=TRUE)




# per site

m2 <- glmmTMB(abundance ~ period * site_short + scaled_temp2 + (1|year2) + (1|interval1) + (1|trap_ID), data= full_data5, ziformula = ~ 0, family = nbinom2)

ggpreds2 <- ggpredict(m2, terms= c("period [all]", "site_short [all"),ci.lvl = F, type="fixed", condition = c(scaled_precip2 = 0,scaled_temp2 =0 ))

ggpreds2$plotting_year[ggpreds2$x == "past"] <- 2
ggpreds2$plotting_year[ggpreds2$x == "present"] <- 4




## plot abundance change ####
# Fig. 1 (part)

plot1 <- ggplot(data= ggpreds1, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "Abundance", x = "")+
  scale_y_continuous( limits = c(0,10), labels= c("0","3","6","     9"), breaks= c(0,3,6,9))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=10, color= colors[9])+
  geom_errorbar(data= ggpreds1, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[9])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "*", y_position = 9),
    textsize = 10, vjust = -0.2,size = 1,
    manual = F) 

plot2 <- ggplot(ggpreds2, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  ylim(c(0,10))+
  labs(y = "", x = "")+
  scale_y_continuous( limits = c(0,10), labels= c("0","3","6","     9"), breaks= c(0,3,6,9))+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds2, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")
  



## BIOMASS #####

# zero-inflated Gamma good solution without data transformation: https://stackoverflow.com/questions/65745148/is-there-a-difference-between-gamma-hurdle-two-part-models-and-zero-inflated-g


m4 <- glmmTMB(biomass ~ period + scaled_precip2 + scaled_temp2 + (1|year2) + (1|interval1) + (1|site_short/trap_ID), data= full_data5, ziformula = ~ period, family = ziGamma(link="log"))

tab_model(m4)

dharma_sim1 <- simulateResiduals(fittedModel = m4, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)
testZeroInflation(dharma_sim1)



# predict mean biomass based on model
ggpreds4 <- ggpredict(m4, terms= c("period [all]"),ci.lvl = .95, type="fe.zi", condition = c(scaled_precip2 = 0,scaled_temp2 =0 ))
ggpreds4$plotting_year <- c(2,4)


## calculate relative change in abundance ##########

# relative change
(ggpreds4[1,]$predicted - ggpreds4[2,]$predicted) / ggpreds4[1,]$predicted


# bootstrap CIs of relative change 

# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model, terms= c("period [all]"), ci.lvl = .95, type="fe.zi", condition = c(scaled_precip2 = 0,scaled_temp2 =0 ))
return(ggpreds$predicted)
}

# quantile approach for CIs
biomass_boot <- bootMer(m4, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed= 129)


boot_rates <- c()
biomass_boot$t

for (i in c(1:1000)){
  boot_rates[i] <-(biomass_boot$t[i,1] - biomass_boot$t[i,2]) / biomass_boot$t[i,1]
}

#lwr   
quantile(boot_rates, probs=.025, na.rm=TRUE)

#upr
quantile(boot_rates, probs=.975, na.rm=TRUE)




#plot


plot3 <- ggplot(data= ggpreds4, aes(x= plotting_year, y= predicted))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  ylim(c(0,2200))+
  labs(y = "Biomass (mg)", x = "")+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_point(size=10, color= colors[9])+
  geom_errorbar(data= ggpreds4, aes( ymin=conf.low, ymax=conf.high), width= 0.3, size=2.5, color= colors[9])+
  geom_signif(    aes(xmin = 2, xmax = 4, annotations = "**", y_position = 2100),
                  textsize = 10, vjust = -0.2,size = 1,
                  manual = F) 



# per site
m5 <- glmmTMB(biomass ~ period*site_short + scaled_precip2 + scaled_temp2 + (1|year2) + (1|interval1) + (1|trap_ID), data= full_data5, ziformula = ~ period, family = ziGamma(link="log"))


ggpreds5 <- ggpredict(m5, terms= c("period [all]", "site_short [all"),ci.lvl = F, type="fe.zi", condition = c(scaled_precip2 = 0,scaled_temp2 =0 ))

ggpreds5$plotting_year[ggpreds5$x == "past"] <- 2
ggpreds5$plotting_year[ggpreds5$x == "present"] <- 4


plot4 <- ggplot(ggpreds5, aes(x= plotting_year, y=predicted, by=group))+
  theme_classic() +
  theme(text = element_text(size = 25)) +
  ylim(c(0,2200))+
  labs(y = "", x = "")+
  scale_x_continuous( limits = c(1,5),labels= c( "1999-2001", "2020-2022"), breaks = c(2,4) )+
  geom_line(data= ggpreds5, aes(x= plotting_year, y=predicted, by=group),lwd=1)+
  geom_point(data= ggpreds5, aes(x= plotting_year, y=predicted, by=group), size = 6, pch=19)+
  theme(legend.title = element_text("site"))+ 
  theme(legend.position="none")





# Fig. 2
grid.arrange(plot1, plot2, plot3, plot4, widths= c(1,1), ncol=2)
# export 10x10




## extract and save site specific trends ######

# abundance
past <- ggpreds2[ggpreds2$x == "past",] 
present <- ggpreds2[ggpreds2$x == "present",] 

site_trends <- data.frame(site= past$group)
site_trends$abundance <- ((past$predicted - present$predicted) / past$predicted) *(-1)

# biomass
past <- ggpreds5[ggpreds5$x == "past",] 
present <- ggpreds5[ggpreds5$x == "present",] 

site_trends$biomass <- ((past$predicted - present$predicted) / past$predicted) *(-1)





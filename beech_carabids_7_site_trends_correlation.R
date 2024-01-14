
# SUPPORTING R CODE FOR
# Evidence for regional-scale drought-induced declines in carabid beetles in old lowland beech forests
# Fabio Weiss, Susanne Winter, Dirk Pflugmacher, Thomas Kolling, Andreas Linde




# TREND AND SITE-VARIABLE CORRELATIONS

# data
setwd("C:/Users/fweiss/Promotion/PhD - own research/projects/Repeated Survey/data/publish data")
site_vars <- read.csv2("site_variables.csv")

# for additional data see previous code file
site_trends <- site_trends1

site_trends <- site_trends[,-1]

site_trends <- site_trends %>% 
  rename(  "Abundance" ="abundance",
           "Biomass"="biomass",
           "Richness"="richness",
           "Simpson"="simpson",
           "Evenness"="evenness",
           "Size.CWM"="cwm_size",
           "Size.CM"="cm_size",
           "Wings.Ind"="ind_wings",
           "Wings.Spec"="spec_wings",
           "Predator.Ind"="ind_pred",
           "Predator.Spec"= "spec_pred")



## colors
colors <- c4a("renoir", 12)
colors[13] <- adjustcolor(colors[11], alpha=0.9)
colors[14] <- adjustcolor(colors[8], alpha=0.9)
colors[15] <- adjustcolor(colors[2], alpha=0.9)

show_col(colors)



### test normality of trends ####

dat <- site_trends[,-1]

shapiro_test_df <- function(df, bonf= T, alpha= 0.05) {
  l <- lapply(df, shapiro.test)
  s <- do.call("c", lapply(l, "[[", 1))
  p <- do.call("c", lapply(l, "[[", 2))
  if (bonf == TRUE) {
    sig <- ifelse(p > alpha / length(l), "H0", "Ha")
  } else {
    sig <- ifelse(p > alpha, "H0", "Ha")
  }
  return(list(statistic= s,
              p.value= p,
              significance= sig,
              method= ifelse(bonf == TRUE, "Shapiro-Wilks test with Bonferroni Correction",
                             "Shapiro-Wilks test without Bonferroni Correction")))
}

test <- shapiro_test_df(dat)



### test correlations among trends ####

dat <- site_trends[,-c(1,8)]

corr.matrix <- corr.test(dat, method = "pearson")
cor_mat <- corr.matrix$r
p_mat <- corr.matrix$p

corrplot(cor_mat, addCoef.col=T, type="lower", p.mat= p_mat , insig = "blank", diag = F, tl.srt=45,	 tl.col="black",  col= colorRampPalette(c(colors[9], "white",  colors[2]))(10))




### test normality of site variables ####

dat <- site_vars[,-1]
shapiro_test_df <- function(df, bonf= T, alpha= 0.05) {
  l <- lapply(df, shapiro.test)
  s <- do.call("c", lapply(l, "[[", 1))
  p <- do.call("c", lapply(l, "[[", 2))
  if (bonf == TRUE) {
    sig <- ifelse(p > alpha / length(l), "H0", "Ha")
  } else {
    sig <- ifelse(p > alpha, "H0", "Ha")
  }
  return(list(statistic= s,
              p.value= p,
              significance= sig,
              method= ifelse(bonf == TRUE, "Shapiro-Wilks test with Bonferroni Correction",
                             "Shapiro-Wilks test without Bonferroni Correction")))
}

test <- shapiro_test_df(dat)



### test correlation among env. vars ####

dat <- site_vars[,-c(1)]



corr.matrix <- corr.test(dat, method = "pearson")

cor_mat <- corr.matrix$r
p_mat <- corr.matrix$p



corrplot(cor_mat, addCoef.col=T, p.mat= p_mat )


corrplot(cor_mat, addCoef.col=T, type="lower", p.mat= p_mat , insig = "blank", diag = F, tl.srt=45,	 tl.col="black",  col= colorRampPalette(c(colors[9], "white",  colors[2]))(10))






### Correlation and t-tests  ########


# correlation 

dat1 <- site_trends[, c(2:12)]

dat2 <- site_vars[,-c(1)]

cor_mat <- cor(dat2,dat1)

corr.matrix <- corr.test(dat2,dat1, method = "pearson")

cor_mat <- corr.matrix$r
p_mat <- corr.matrix$p


corrplot(cor_mat, method= "circle", addCoef.col=T, p.mat= p_mat , insig = "blank", 	tl.srt=45, tl.col="black",  col= colorRampPalette(c(colors[9], "white",  colors[2]))(10))
# replace row with protection with t-test results in PP


# t-tests
ttest_dat <- site_trends
ttest_dat$protection <- site_vars$Protection[match(ttest_dat$site, site_vars$site)]

for ( i in c(2:12)){
  
  dat <- data.frame(value= ttest_dat[,i], group= ttest_dat$protection)

  test<- t_test(data=dat, value ~ group, paired= F)
  
  print(paste(colnames(ttest_dat[i]), ", p=", test$p , sep = " "))
  
}


### Bootstrapping means and correlations ###########


#  Lat and size CWM ###############
dat <-data.frame(lat= site_vars$Latitude, cwm_size = site_trends$Size.CWM)

plot(dat$cwm_size  ~ dat$lat)

# correlation coef
cor1 <- cor.test(dat$cwm_size,  dat$lat, method="pearson", paired=F)
cor1$estimate

# bootstrap CIs
set.seed(128)

n <- 10000
BootstrapSamples <- data.frame(lat = c(), cwm_size = c(), boot= c())

for( i in c(1:n)){
  boot_data <- dat[sample(nrow(dat), 11, replace = T), ]
  boot_data$boot <- i
  
  BootstrapSamples <- rbind(BootstrapSamples, boot_data)
}

cor_boot_samples <- c() 

for(z in c(1:n)){
  loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
  cor_boot <- cor.test(loop_data$cwm_size,  loop_data$lat, method="pearson", paired=F) 
  cor_boot_samples[z] <- cor_boot$estimate
}

#upr  
quantile(cor_boot_samples, probs=.975, na.rm=TRUE)  

#lwr
quantile(cor_boot_samples, probs=.025, na.rm=TRUE) 


cor_text <- paste("r=", round(cor1$estimate, 3), " [", round(quantile(cor_boot_samples, probs=.025, na.rm=TRUE), 3), ",", round(quantile(cor_boot_samples, probs=.975, na.rm=TRUE), 3), "]", sep = ""   ) 

# p value 

set.seed(128)

n <- 10000
BootstrapSamples <- data.frame(lat = c(), cwm_size = c(), boot= c())

for( i in c(1:n)){
  
  boot_data <- data.frame(lat = c(rep(NA,11)), cwm_size = c(rep(NA,11)), boot= c(rep(NA,11)))
  boot_data$lat <- sample(dat$lat, 11, replace = T)
  boot_data$cwm_size <- sample(dat$cwm_size, 11, replace = T) 
  boot_data$boot <- i
  
  BootstrapSamples <- rbind(BootstrapSamples, boot_data)
}

cor_boot_samples <- c() 

for(z in c(1:n)){
  loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
  cor_boot <- cor.test(loop_data$cwm_size,  loop_data$lat, method="pearson", paired=F) 
  cor_boot_samples[z] <- cor_boot$estimate
}

# p value
pval<- paste ( "p=", mean(cor_boot_samples <= cor1$estimate), sep="")
# p= 0.043

# plot 
# export 6x4
ggplot(dat, aes(x= lat, y=cwm_size))+
  theme_classic() +
  theme(text = element_text(size = 20)) +
  ggtitle("")+
  labs(y = "Mean trend in size (CWM)", x = "Latitude")+
  #scale_x_continuous(labels= c("total reserves", "managed"), breaks = c(0,1) )+
  geom_point(size=8, alpha=0.8)+
  geom_smooth(method='lm', se=FALSE, formula= y ~ x, color=colors[9], size=4)+
  annotate(geom="text", x=5890000, y=1.8, label=cor_text, size= 6)+
  annotate(geom="text", x=5890000, y=1, label=pval, size= 6)





#  mgmt and cwm_size ###################
dat <-data.frame(mgmt= site_vars$Protection, cwm_size = site_trends$Size.CWM)

boxplot(dat$cwm_size ~ dat$mgmt)

#means
means <- data.frame (mgmt = c(0,1), means = with(dat, tapply(cwm_size, mgmt, mean)))


# bootstrap CIs

set.seed(128)

n <- length(dat$cwm_size[dat$mgmt=="0"])
B <- 10000
variable <- dat$cwm_size[dat$mgmt=="0"]
BootstrapSamples1 <- matrix(sample(variable, size = n*B, replace =T), nrow=n, ncol=B)

boot_mean1 <- rep(0,B)
for(i in 1:B){
  boot_mean1[i] <- mean(BootstrapSamples1[i])
}


n <- length(dat$cwm_size[dat$mgmt=="1"])
B <- 10000
variable <- dat$cwm_size[dat$mgmt=="1"]
BootstrapSamples2 <- matrix(sample(variable, size = n*B, replace =T), nrow=n, ncol=B)

boot_mean2 <- rep(0,B)
for(i in 1:B){
  boot_mean2[i] <- mean(BootstrapSamples2[i])
}

means$upr <- c(quantile(boot_mean1, probs=.975, na.rm=TRUE), quantile(boot_mean2, probs=.975, na.rm=TRUE))
means$lwr <- c(quantile(boot_mean1, probs=.025, na.rm=TRUE), quantile(boot_mean2, probs=.025, na.rm=TRUE))
                              


# p value

set.seed(128)

diff1 <- abs(diff(with(dat, tapply(cwm_size, mgmt, mean))))
n <- length(dat$mgmt)
B <- 10000
variable <- dat$cwm_size
BootstrapSamples <- matrix(sample(variable, size = n*B, replace =T), nrow=n, ncol=B)

boot_diff <- rep(0,B)
for(i in 1:B){
  boot_diff[i] <- abs(mean(BootstrapSamples[1:5,i]) - mean(BootstrapSamples[6:11, i]))
}

# p value
pval<- paste ( "p=", mean(boot_diff >= diff1), sep="")
# p= 0.043

means$plotcat <- c(2,4)

# plot
ggplot(means, aes(x= plotcat, y=means))+
  theme_classic() +
  theme(text = element_text(size = 20)) +
  ggtitle("")+
  labs(y = "Mean trend in size (CWM)", x = "Protection")+
  scale_x_continuous(limits = c(1,5), labels= c( "managed", "total reserves"), breaks = c(2,4) )+
  geom_point(size=12, color= "black")+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width= 0.3, size=1, color= "black")+ 
 # geom_jitter(data=site_trends, aes(x= mgmt, y=cwm_size), size=4, width= 0.03)+
  annotate(geom="text", x=3, y=-1.5, label=pval, size= 6)

# export 6x4  




#  comm_size and size###############

show_col(colors)

# try mgmt and cwm_size
dat <-data.frame(comm_size= site_vars$Initial.community.size, size = site_trends$Size.CWM)

plot(dat$size  ~ dat$comm_size)

# correlation coef
cor1 <- cor.test(dat$size,  dat$comm_size, method="pearson", paired=F)
cor1$estimate

# bootstrap CIs
set.seed(128)

n <- 10000
BootstrapSamples <- data.frame(comm_size = c(), size = c(), boot= c())

for( i in c(1:n)){
  boot_data <- dat[sample(nrow(dat), 11, replace = T), ]
  boot_data$boot <- i
  
  BootstrapSamples <- rbind(BootstrapSamples, boot_data)
}

cor_boot_samples <- c() 

for(z in c(1:n)){
  loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
  cor_boot <- cor.test(loop_data$size,  loop_data$comm_size, method="pearson", paired=F) 
  cor_boot_samples[z] <- cor_boot$estimate
}

#upr  
quantile(cor_boot_samples, probs=.975, na.rm=TRUE)  

#lwr
quantile(cor_boot_samples, probs=.025, na.rm=TRUE) 


cor_text <- paste("r=", round(cor1$estimate, 3), " [", round(quantile(cor_boot_samples, probs=.025, na.rm=TRUE), 3), ",", round(quantile(cor_boot_samples, probs=.975, na.rm=TRUE), 3), "]", sep = ""   ) 

# p value 

set.seed(128)

n <- 10000
BootstrapSamples <- data.frame(comm_size = c(), size = c(), boot= c())

for( i in c(1:n)){
  
  boot_data <- data.frame(comm_size = c(rep(NA,11)), size = c(rep(NA,11)), boot= c(rep(NA,11)))
  boot_data$comm_size <- sample(dat$comm_size, 11, replace = T)
  boot_data$size <- sample(dat$size, 11, replace = T) 
  boot_data$boot <- i
  
  BootstrapSamples <- rbind(BootstrapSamples, boot_data)
}

cor_boot_samples <- c() 

for(z in c(1:n)){
  loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
  cor_boot <- cor.test(loop_data$size,  loop_data$comm_size, method="pearson", paired=F) 
  cor_boot_samples[z] <- cor_boot$estimate
}

# p value
pval<- paste ( "p=", mean(cor_boot_samples <= cor1$estimate), sep="")



# plot 
# export 6x4
ggplot(dat, aes(x= comm_size, y=size))+
  theme_classic() +
  theme(text = element_text(size = 20)) +
  ggtitle("")+
  labs(y = "Mean trend in size (CWM)", x = "Mean size in 1999-2001 (CWM)")+
  #scale_x_continuous(labels= c("total reserves", "managed"), breaks = c(0,1) )+
  geom_point(size=8, alpha=0.8)+
  geom_smooth(method='lm', se=FALSE, formula= y ~ x, color=colors[9], size=4)+
  annotate(geom="text", x=18.2, y=0.5, label=cor_text, size= 6)+
  annotate(geom="text", x=18.2, y=-0.55, label=pval, size= 6)




#  comm_size and biomass###############

show_col(colors)
  
  # try mgmt and cwm_size
  dat <-data.frame(comm_size= site_vars$Initial.community.size, biomass = site_trends$Biomass)
  
  plot(dat$biomass  ~ dat$comm_size)
  
# correlation coef
cor1 <- cor.test(dat$biomass,  dat$comm_size, method="pearson", paired=F)
cor1$estimate
  
# bootstrap CIs
  set.seed(128)

  n <- 10000
  BootstrapSamples <- data.frame(comm_size = c(), biomass = c(), boot= c())
  
  for( i in c(1:n)){
    boot_data <- dat[sample(nrow(dat), 11, replace = T), ]
    boot_data$boot <- i
    
    BootstrapSamples <- rbind(BootstrapSamples, boot_data)
  }

  cor_boot_samples <- c() 
  
  for(z in c(1:n)){
    loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
    cor_boot <- cor.test(loop_data$biomass,  loop_data$comm_size, method="pearson", paired=F) 
    cor_boot_samples[z] <- cor_boot$estimate
  }

#upr  
 quantile(cor_boot_samples, probs=.975, na.rm=TRUE)  

#lwr
 quantile(cor_boot_samples, probs=.025, na.rm=TRUE) 

  
  cor_text <- paste("r=", round(cor1$estimate, 3), " [", round(quantile(cor_boot_samples, probs=.025, na.rm=TRUE), 3), ",", round(quantile(cor_boot_samples, probs=.975, na.rm=TRUE), 3), "]", sep = ""   ) 
  
# p value 

  set.seed(128)
  
  n <- 10000
  BootstrapSamples <- data.frame(comm_size = c(), biomass = c(), boot= c())
  
  for( i in c(1:n)){
    
    boot_data <- data.frame(comm_size = c(rep(NA,11)), biomass = c(rep(NA,11)), boot= c(rep(NA,11)))
    boot_data$comm_size <- sample(dat$comm_size, 11, replace = T)
    boot_data$biomass <- sample(dat$biomass, 11, replace = T) 
    boot_data$boot <- i
    
    BootstrapSamples <- rbind(BootstrapSamples, boot_data)
  }
  
  cor_boot_samples <- c() 
  
  for(z in c(1:n)){
    loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
    cor_boot <- cor.test(loop_data$biomass,  loop_data$comm_size, method="pearson", paired=F) 
    cor_boot_samples[z] <- cor_boot$estimate
  }
  
    # p value
    pval<- paste ( "p=", mean(cor_boot_samples <= cor1$estimate), sep="")
    # p= 0.043

  
  
# plot 
# export 6x4
    ggplot(dat, aes(x= comm_size, y=biomass))+
      theme_classic() +
      theme(text = element_text(size = 20)) +
      ggtitle("")+
      labs(y = "Mean trend in biomass", x = "Mean size in 1999-2001 (CWM)")+
      #scale_x_continuous(labels= c("total reserves", "managed"), breaks = c(0,1) )+
      geom_point(size=8, alpha=0.8)+
      geom_smooth(method='lm', se=FALSE, formula= y ~ x, color=colors[9], size=4)+
      annotate(geom="text", x=18.2, y=-0.5, label=cor_text, size= 6)+
      annotate(geom="text", x=18.2, y=-0.55, label=pval, size= 6)


#  cover and biomass###############
    
    
    
    dat <-data.frame(cover= site_vars$Landscape.cover, biomass = site_trends$Biomass)
    
    plot(dat$biomass  ~ dat$cover)
    
    # correlation coef
    cor1 <- cor.test(dat$biomass,  dat$cover, method="pearson", paired=F)
    cor1$estimate
    
    # bootstrap CIs
    set.seed(128)
    
    n <- 10000
    BootstrapSamples <- data.frame(cover = c(), biomass = c(), boot= c())
    
    for( i in c(1:n)){
      boot_data <- dat[sample(nrow(dat), 11, replace = T), ]
      boot_data$boot <- i
      
      BootstrapSamples <- rbind(BootstrapSamples, boot_data)
    }
    
    cor_boot_samples <- c() 
    
    for(z in c(1:n)){
      loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
      cor_boot <- cor.test(loop_data$biomass,  loop_data$cover, method="pearson", paired=F) 
      cor_boot_samples[z] <- cor_boot$estimate
    }
    
    #upr  
    quantile(cor_boot_samples, probs=.975, na.rm=TRUE)  
    
    #lwr
    quantile(cor_boot_samples, probs=.025, na.rm=TRUE) 
    
    
    cor_text <- paste("r=", round(cor1$estimate, 3), " [", round(quantile(cor_boot_samples, probs=.025, na.rm=TRUE), 3), ",", round(quantile(cor_boot_samples, probs=.975, na.rm=TRUE), 3), "]", sep = ""   ) 
    
    # p value 
    
    set.seed(128)
    
    n <- 10000
    BootstrapSamples <- data.frame(cover = c(), biomass = c(), boot= c())
    
    for( i in c(1:n)){
      
      boot_data <- data.frame(cover = c(rep(NA,11)), biomass = c(rep(NA,11)), boot= c(rep(NA,11)))
      boot_data$cover <- sample(dat$cover, 11, replace = T)
      boot_data$biomass <- sample(dat$biomass, 11, replace = T) 
      boot_data$boot <- i
      
      BootstrapSamples <- rbind(BootstrapSamples, boot_data)
    }
    
    cor_boot_samples <- c() 
    
    for(z in c(1:n)){
      loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
      cor_boot <- cor.test(loop_data$biomass,  loop_data$cover, method="pearson", paired=F) 
      cor_boot_samples[z] <- cor_boot$estimate
    }
    
    # p value
    pval<- paste ( "p=", mean(cor_boot_samples >= cor1$estimate), sep="")
    # p= 0.043
    
    
    
    # plot 
    # export 6x4
    
    ggplot(dat, aes(x= cover, y=biomass))+
      theme_classic() +
      theme(text = element_text(size = 20)) +
      ggtitle("")+
      labs(y = "Mean trend in biomass", x = "Landscape (forest) cover")+
      #scale_x_continuous(labels= c("total reserves", "managed"), breaks = c(0,1) )+
      geom_point(size=8, alpha=0.8)+
      geom_smooth(method='lm', se=FALSE, formula= y ~ x, color=colors[2], size=4)+
      annotate(geom="text", x=75, y=-0.5, label=cor_text, size= 6)+
      annotate(geom="text", x=75, y=-0.55, label=pval, size= 6)
    
  
    
#  canopy and wings ###############
    
    
    
    dat <-data.frame(canopy= site_vars$Canopy, wings = site_trends$Wings.Spec)
    
    plot(dat$wings  ~ dat$canopy)
    
    # correlation coef
    cor1 <- cor.test(dat$wings,  dat$canopy, method="pearson", paired=F)
    cor1$estimate
    
    # bootstrap CIs
    set.seed(128)
    
    n <- 10000
    BootstrapSamples <- data.frame(canopy = c(), wings = c(), boot= c())
    
    for( i in c(1:n)){
      boot_data <- dat[sample(nrow(dat), 11, replace = T), ]
      boot_data$boot <- i
      
      BootstrapSamples <- rbind(BootstrapSamples, boot_data)
    }
    
    cor_boot_samples <- c() 
    
    for(z in c(1:n)){
      loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
      cor_boot <- cor.test(loop_data$wings,  loop_data$canopy, method="pearson", paired=F) 
      cor_boot_samples[z] <- cor_boot$estimate
    }
    
    #upr  
    quantile(cor_boot_samples, probs=.975, na.rm=TRUE)  
    
    #lwr
    quantile(cor_boot_samples, probs=.025, na.rm=TRUE) 
    
    
    cor_text <- paste("r=", round(cor1$estimate, 3), " [", round(quantile(cor_boot_samples, probs=.025, na.rm=TRUE), 3), ",", round(quantile(cor_boot_samples, probs=.975, na.rm=TRUE), 3), "]", sep = ""   ) 
    
    # p value 
    
    set.seed(128)
    
    n <- 10000
    BootstrapSamples <- data.frame(canopy = c(), wings = c(), boot= c())
    
    for( i in c(1:n)){
      
      boot_data <- data.frame(canopy = c(rep(NA,11)), wings = c(rep(NA,11)), boot= c(rep(NA,11)))
      boot_data$canopy <- sample(dat$canopy, 11, replace = T)
      boot_data$wings <- sample(dat$wings, 11, replace = T) 
      boot_data$boot <- i
      
      BootstrapSamples <- rbind(BootstrapSamples, boot_data)
    }
    
    cor_boot_samples <- c() 
    
    for(z in c(1:n)){
      loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
      cor_boot <- cor.test(loop_data$wings,  loop_data$canopy, method="pearson", paired=F) 
      cor_boot_samples[z] <- cor_boot$estimate
    }
    
    # p value
    pval<- paste ( "p=", mean(cor_boot_samples >= cor1$estimate), sep="")
    # p= 0.043
    
    
    
    # plot 
    # export 6x4
    
    ggplot(dat, aes(x= canopy, y=wings))+
      theme_classic() +
      theme(text = element_text(size = 18)) +
      ggtitle("")+
      labs(y = "Change of prop. for flying species", x = "Canopy cover")+
      #scale_x_continuous(labels= c("total reserves", "managed"), breaks = c(0,1) )+
      geom_point(size=8, alpha=0.8)+
      geom_smooth(method='lm', se=FALSE, formula= y ~ x, color=colors[2], size=4)+
      annotate(geom="text", x=40, y=0.35, label=cor_text, size= 6)+
      annotate(geom="text", x=40, y=0.25, label=pval, size= 6)
    
    
    
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


### visualisation

# https://towardsdatascience.com/how-to-create-a-correlation-matrix-with-too-many-variables-309cc0c0a57

# Mann-Whitney test for management
# spearman or pearson for other variables



### Bootstrapping??
# https://www.uvm.edu/~statdhtx/StatPages/Randomization%20Tests/BootstCorr/bootstrapping_correlations.html
# https://www.youtube.com/watch?v=9STZ7MxkNVg
# https://www.youtube.com/watch?v=Zet-qmEEfCU












    
    
#  canopy and predators ###############
    
    
    
    dat <-data.frame(canopy= site_vars$Canopy, predator = site_trends$Predator.Spec)
    
    plot(dat$predator  ~ dat$canopy)
    
    # correlation coef
    cor1 <- cor.test(dat$predator,  dat$canopy, method="pearson", paired=F)
    cor1$estimate
    
    # bootstrap CIs
    set.seed(128)
    
    n <- 10000
    BootstrapSamples <- data.frame(canopy = c(), predator = c(), boot= c())
    
    for( i in c(1:n)){
      boot_data <- dat[sample(nrow(dat), 11, replace = T), ]
      boot_data$boot <- i
      
      BootstrapSamples <- rbind(BootstrapSamples, boot_data)
    }
    
    cor_boot_samples <- c() 
    
    for(z in c(1:n)){
      loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
      cor_boot <- cor.test(loop_data$predator,  loop_data$canopy, method="pearson", paired=F) 
      cor_boot_samples[z] <- cor_boot$estimate
    }
    
    #upr  
    quantile(cor_boot_samples, probs=.975, na.rm=TRUE)  
    
    #lwr
    quantile(cor_boot_samples, probs=.025, na.rm=TRUE) 
    
    
    cor_text <- paste("r=", round(cor1$estimate, 3), " [", round(quantile(cor_boot_samples, probs=.025, na.rm=TRUE), 3), ",", round(quantile(cor_boot_samples, probs=.975, na.rm=TRUE), 3), "]", sep = ""   ) 
    
    # p value 
    
    set.seed(128)
    
    n <- 10000
    BootstrapSamples <- data.frame(canopy = c(), predator = c(), boot= c())
    
    for( i in c(1:n)){
      
      boot_data <- data.frame(canopy = c(rep(NA,11)), predator = c(rep(NA,11)), boot= c(rep(NA,11)))
      boot_data$canopy <- sample(dat$canopy, 11, replace = T)
      boot_data$predator <- sample(dat$predator, 11, replace = T) 
      boot_data$boot <- i
      
      BootstrapSamples <- rbind(BootstrapSamples, boot_data)
    }
    
    cor_boot_samples <- c() 
    
    for(z in c(1:n)){
      loop_data <- BootstrapSamples[BootstrapSamples$boot == z,]
      cor_boot <- cor.test(loop_data$predator,  loop_data$canopy, method="pearson", paired=F) 
      cor_boot_samples[z] <- cor_boot$estimate
    }
    
    # p value
    pval<- paste ( "p=", mean(cor_boot_samples <= cor1$estimate), sep="")
    # p= 0.043
    
    
    
    # plot 
    # export 6x4
    
    ggplot(dat, aes(x= canopy, y=predator))+
      theme_classic() +
      theme(text = element_text(size = 17)) +
      ggtitle("")+
      labs(y = "Change of prop. for predatory species", x = "Canopy cover")+
      #scale_x_continuous(labels= c("total reserves", "managed"), breaks = c(0,1) )+
      geom_point(size=8, alpha=0.8)+
      geom_smooth(method='lm', se=FALSE, formula= y ~ x, color=colors[9], size=4)+
      annotate(geom="text", x=60, y=0.18, label=cor_text, size= 6)+
      annotate(geom="text", x=60, y=0.12, label=pval, size= 6)
    
    

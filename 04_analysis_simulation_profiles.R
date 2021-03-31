#title: Construction of migratory connectivity profiles from simulation from Knight et al. 2021. Comprehensive estimation of spatial and temporal migratory connectivity across the annual cycle to direct conservation efforts. Ecography ECOG-05111
#author: Elly C. Knight
#date: Aug. 19, 2020

#1. PRELIMINARY####
#load packages
library(tidyverse)
library(mgcv)
library(sf)
library(pracma)
library(meanShiftR)

#limit display of scientific notation
options(scipen = 999)

#Change working directory to subfolder
setwd(paste0(getwd(),"/simulation"))

#2. SPATIAL####

##2a. Find peaks----
files.spat <- list.files(paste0(rd, "MCSpatial"))
peaks.spat.final <- data.frame()
for(i in 1:length(files.spat)){
  
  spat <- read.csv(paste0(rd, "MCSpatial/", files.spat[i])) %>% 
    dplyr::rename(M=mantel) %>% 
    dplyr::filter(!is.na(M))
  
  spat.sum <- spat %>% 
    group_by(lat, l1o) %>% 
    summarize(Mmean=mean(M),
              Mlow83=quantile(M, probs=0.083),
              Mhigh83=quantile(M, probs=0.917),
              Mlow95=quantile(M, probs=0.025),
              Mhigh95=quantile(M, probs=0.975),
              Mlow99=quantile(M, probs=0.005),
              Mhigh99=quantile(M, probs=0.995),
              Mlow=min(M),
              Mhigh=max(M),
              Mlowsd=mean(M)-sd(M),
              Mhighsd=mean(M)+sd(M)) %>% 
    ungroup()
  
  gamlist <- expand.grid(l1o=unique(spat$l1o),
                         k=c(5:15))
  
  pred.spat.all <- data.frame()
  for(j in 1:nrow(gamlist)){
    
    spat.j <- spat %>% 
      filter(l1o==gamlist$l1o[j])
    
    gam.j <- gam(M ~ s(lat, k = gamlist$k[j], bs = 'cs'), data=spat.j)
    
    pred.j <- data.frame(predict(gam.j, newdata=data.frame(lat=min(spat.j$lat):max(spat.j$lat)), se.fit=TRUE)) %>% 
      cbind(data.frame(lat=min(spat.j$lat):max(spat.j$lat))) %>% 
      mutate(upr = fit + (1.96*se.fit),
             lwr = fit - (1.96*se.fit),
             l1o = gamlist$l1o[j],
             k = gamlist$k[j])
    
    pred.spat.all <- rbind(pred.spat.all, pred.j)
    
  }
  
  ndoys <- table(spat.sum$l1o) %>% 
    data.frame() %>% 
    rename(l1o=Var1, lats=Freq) %>% 
    mutate(l1o=as.integer(as.character(l1o))) 
  
  pred.spat.k <- pred.spat.all %>% 
    inner_join(spat.sum, by=c("lat", "l1o")) %>% 
    mutate(use=ifelse(fit < Mhigh99 & fit > Mlow99, 1, 0)) %>% 
    group_by(l1o, k) %>% 
    summarize(sum=sum(use)) %>% 
    left_join(ndoys, by="l1o") %>% 
    dplyr::mutate(diff=lats-sum) %>% 
    dplyr::filter(diff==min(diff)) %>% 
    group_by(k) %>% 
    summarize(pops=n()) %>% 
    ungroup() %>% 
    dplyr::filter(pops==max(pops))
  
  spat.k <- min(pred.spat.k$k)
  
  print(paste0("The selected number of knots is ", spat.k))
  
  pred.spat <- pred.spat.all %>% 
    dplyr::filter(k==spat.k)
  
  ggplot() +
    geom_hex(aes(x=lat, y=M), data=spat) +
    geom_line(aes(x=lat, y=fit), data=pred.spat) +
    facet_wrap(~l1o)
  
  peaks.spat <- data.frame()
  for(j in 1:length(unique(pred.spat$l1o))){
    
    l1o.j <- unique(pred.spat$l1o)[j]
    
    pred.j <- pred.spat %>% 
      dplyr::filter(l1o==l1o.j)
    
    max.ind.j <- findpeaks(pred.j$fit, nups = 1, ndowns=1)
    min.ind.j <- findpeaks(-pred.j$fit, nups = 1, ndowns=1)
    
    max.j <- pred.j[max.ind.j[,2],] %>% 
      mutate(peak="max") %>% 
      data.frame() %>% 
      mutate(latlagdiff = abs(lat-lag(lat)),
             latleaddiff = abs(lat-lead(lat)),
             use = case_when(latlagdiff < 5 & fit < lag(fit) ~ 0,
                             latleaddiff < 5 & fit < lead(fit) ~ 0),
             use = ifelse(is.na(use), 1, use))
    
    min.j <- pred.j[min.ind.j[,2],] %>% 
      mutate(peak="min") %>% 
      data.frame() %>% 
      mutate(latlagdiff = NA,
             latleaddiff = NA,
             use = NA)
    
    peaks.spat <- rbind(peaks.spat, max.j, min.j) %>% 
      arrange(l1o, lat)
    
  }
  
  peaks.spat <- peaks.spat %>% 
    group_by(l1o) %>% 
    mutate(latlagdiff=as.numeric(latlagdiff),
           latleaddiff=as.numeric(latleaddiff),
           use=as.numeric(use)) %>% 
    mutate(use=case_when(peak=="min" & lead(use)==1 & lag(use)==0 ~ 0,
                         peak=="min" & lead(use)==0 & lag(use)==1 ~ 0,
                         !is.na(use) ~ use),
           use = ifelse(is.na(use), 1, use)) %>% 
    dplyr::filter(use!=0) %>% 
    mutate(id=row_number()) %>% 
    ungroup()
  
  mat.spat1 <- matrix((peaks.spat %>% 
                         dplyr::filter(peak=="max"))$lat)
  mat.spat2 <- matrix((peaks.spat %>% 
                         dplyr::filter(peak=="max"))$fit)
  mat.spat <- cbind(mat.spat1, mat.spat2)
  
  spat.ms <- meanShift(mat.spat,
                               algorithm="LINEAR",
                               bandwidth=c(2,2))
  peaks.spat.ms <- peaks.spat %>% 
    dplyr::filter(peak=="max") %>% 
    mutate(cluster = spat.ms$assignment)
  
  ggplot() +
    geom_point(aes(x=lat, y=fit, colour=factor(cluster)), data=peaks.spat.ms)
  
  if(nrow(peaks.spat.ms)>0){
    
    spat.clust <- table(peaks.spat.ms$l1o, peaks.spat.ms$cluster) %>% 
      data.frame() %>% 
      rename(l1o=Var1, cluster=Var2) %>% 
      group_by(cluster) %>% 
      summarize(count=sum(Freq)) %>% 
      ungroup() %>% 
      dplyr::filter(count >= length(unique(spat$l1o)))
    
    peaks.spat.clust <- peaks.spat.ms %>% 
      dplyr::filter(cluster %in% spat.clust$cluster)
    peaks.spat.clust
    
    peaks.spat.ci <- peaks.spat %>%
      left_join(spat.sum, by=c("lat", "l1o")) %>% 
      dplyr::filter(l1o==99) %>% 
      mutate(Mlagdiff = round((Mmean - lag(Mmean)), 1),
             Mleaddiff = round((Mmean - lead(Mmean)), 1),
             lagci = ifelse(Mlow83 < lag(Mhigh83), 0, 1),
             leadci = ifelse(Mhigh83 < lead(Mhigh83), 0, 1),
             lagci = ifelse(is.na(lagci), 99, lagci),
             leadci = ifelse(is.na(leadci), 99, leadci)) %>% 
      mutate(Mlagdiff = ifelse(is.na(Mlagdiff), 99, Mlagdiff),
             Mleaddiff= ifelse(is.na(Mleaddiff), 99, Mleaddiff)) %>% 
      dplyr::filter(peak=="max",
                    lagci!=0,
                    leadci!=0,
                    !(Mlagdiff <= 0.1 & Mleaddiff <= 0.1))
    peaks.spat.ci
    
    peaks.spat.all <- peaks.spat.clust %>% 
      inner_join(peaks.spat.ci, by = c("fit", "se.fit", "lat", "upr", "lwr", "l1o", "k", "peak", "id")) %>% 
      mutate(file=files.spat[i])  %>% 
      separate(file, into=c("metric", "connectivity", "n.pop", "n.ind", "n.gps", "sim", ".filetype"), remove=FALSE)
    
    if(nrow(peaks.spat.all) >0){
      
      raw <- read.csv(paste0("Data/StartingPoints/SimulationStartingPoints_",
                             unique(peaks.spat.all$n.pop),
                             "-",
                             unique(peaks.spat.all$n.ind),
                             "-",
                             unique(peaks.spat.all$n.gps),
                             "_",
                             unique(peaks.spat.all$sim),
                             ".csv"))
      raw <- raw %>% 
        sf::st_as_sf(coords = c("X","Y"), crs=3857) %>% 
        sf::st_transform(4326) %>% 
        st_coordinates() %>% 
        data.frame() %>% 
        dplyr::rename(lat=Y, long=X) %>% 
        cbind(raw)
        
      stop2lat <- raw %>% 
        dplyr::filter(season=="stop2") %>% 
        dplyr::select(pop, lat) %>% 
        unique()
      
      mean(stop2lat$lat)
      peaks.spat.final.i <- peaks.spat.all %>% 
        mutate(true = ifelse(lat <= (mean(stop2lat$lat+5)) & lat >= (mean(stop2lat$lat-5)), 1, 0))
      
      peaks.spat.final <- rbind(peaks.spat.final, peaks.spat.final.i)
      
      
    }
    
    
  }
  
  else{
    next
  }

  print(paste0("Completed processing simulation ", i, " of ", length(files.spat)))
  
}

##2b. Wrangle----
peaks.spat.final.clean <- peaks.spat.final %>% 
  full_join(data.frame(file=files.spat)) %>% 
  separate(file, into=c("metric", "connectivity", "n.pop", "n.ind", "n.gps", "sim", ".filetype"), remove=FALSE) %>% 
  mutate(detect = ifelse(is.na(cluster), 0, 1),
         true = ifelse(is.na(true), 1, true)) %>% 
  dplyr::select(file, n.pop, n.ind, n.gps, sim, lat, Mmean, Mlow83, Mhigh83, Mlow95, Mhigh95, true, detect) %>% 
  mutate(status=case_when(true==1 & detect==1 ~ "tp",
                          true==1 & detect==0 ~ "fn",
                          true==0 & detect==1 ~ "fp",
                          true==0 & detect==0 ~ "tn"),
         sim = as.numeric(sim),
         n.pop = as.numeric(n.pop),
         n.ind = as.numeric(n.ind),
         n.gps = as.numeric(n.gps)) %>% 
  arrange(n.pop, n.ind, n.gps, sim) %>% 
#  dplyr::filter(n.gps <= 50) %>% 
  unique() 

write.csv(peaks.spat.final.clean, "CleanedSpatialPeaks_All.csv", row.names = FALSE)

##2c. Calculate probability of false negative----
peaks.spat.final.clean <- read.csv("CleanedSpatialPeaks_All.csv") %>% 
  dplyr::filter(n.gps <= 50)

#Data
peaks.spat.final.fn <- peaks.spat.final.clean %>% 
  dplyr::filter(true==1)

#Rate
table(peaks.spat.final.fn$detect)

#Visualize
ggplot(peaks.spat.final.fn) +
  geom_jitter(aes(x=n.ind, y=detect)) +
  geom_smooth(aes(x=n.ind, y=detect), method="lm")

ggplot(peaks.spat.final.fn) +
  geom_jitter(aes(x=n.gps, y=detect)) +
  geom_smooth(aes(x=n.gps, y=detect), method="lm")

#Model
fn1 <- glm(detect ~ n.ind+n.gps, data=peaks.spat.final.fn, family="binomial")
summary(fn1)

#Predict
#Individuals
newdat <- data.frame(n.ind=seq(min(peaks.spat.final.fn$n.ind), max(peaks.spat.final.fn$n.ind),1),
                     n.gps=round(mean(peaks.spat.final.fn$n.gps)))

fp1pred <- predict(fn1, newdat, se=TRUE, type="response")

pred.spat.fn.ind <- data.frame(fit=fp1pred$fit,
                         se=fp1pred$se.fit) %>% 
  cbind(newdat) %>% 
  mutate(diff = fit-lag(fit),
         con = "spatial")

#GPS points
newdat <- data.frame(n.gps=seq(min(peaks.spat.final.fn$n.gps), max(peaks.spat.final.fn$n.gps),1),
                     n.ind=round(mean(peaks.spat.final.fn$n.ind)))

fp1pred <- predict(fn1, newdat, se=TRUE, type="response")

pred.spat.fn.gps <- data.frame(fit=fp1pred$fit,
                           se=fp1pred$se.fit) %>% 
  cbind(newdat) %>% 
  mutate(diff = fit-lag(fit),
         con = "spatial")

#Summarize predictions
summary(pred.spat.fn.ind$diff)
summary(pred.spat.fn.gps$diff)

#Plot predictions
ggplot(pred.spat.fn.ind) +
  geom_line(aes(x=n.ind, y=fit)) +
  geom_ribbon(aes(x=n.ind, ymax=fit+2*se, ymin=fit-2*se), alpha=0.5)

ggplot(pred.spat.fn.gps) +
  geom_line(aes(x=n.gps, y=fit)) +
  geom_ribbon(aes(x=n.gps, ymax=fit+2*se, ymin=fit-2*se), alpha=0.5)

##2d. Calcualte probability of false positive----
#Data
peaks.spat.final.fp <- peaks.spat.final.clean %>% 
  dplyr::filter(status=="fp") %>% 
  dplyr::select(file) %>% 
  unique() %>% 
  mutate(detect=1) %>% 
  full_join(data.frame(file=files.spat)) %>% 
  separate(file, into=c("metric", "connectivity", "n.pop", "n.ind", "n.gps", "sim", ".filetype"), remove=FALSE) %>% 
  mutate(sim = as.numeric(sim),
         n.pop = as.numeric(n.pop),
         n.ind = as.numeric(n.ind),
         n.gps = as.numeric(n.gps),
         detect=ifelse(is.na(detect), 0, detect)) %>% 
  arrange(n.pop, n.ind, n.gps, sim) %>% 
  dplyr::filter(n.gps <= 50) %>% 
  unique() 

#Rate
table(peaks.spat.final.fp$detect)

#Visualize
ggplot(peaks.spat.final.fp) +
  geom_jitter(aes(x=n.ind, y=detect)) +
  geom_smooth(aes(x=n.ind, y=detect), method="lm")

ggplot(peaks.spat.final.fp) +
  geom_jitter(aes(x=n.gps, y=detect)) +
  geom_smooth(aes(x=n.gps, y=detect), method="lm")

#Model
fp1 <- glm(detect ~ n.ind+n.gps, data=peaks.spat.final.fp, family="binomial")
summary(fp1)

#Predict
#Individuals
newdat <- data.frame(n.ind=seq(min(peaks.spat.final.fp$n.ind), max(peaks.spat.final.fp$n.ind),1),
                     n.gps=round(mean(peaks.spat.final.fp$n.gps)))

fp1pred <- predict(fp1, newdat, se=TRUE, type="response")

pred.spat.fp.ind <- data.frame(fit=fp1pred$fit,
                               se=fp1pred$se.fit) %>% 
  cbind(newdat) %>% 
  mutate(diff = fit-lag(fit),
         con = "spatial")

#GPS points
newdat <- data.frame(n.gps=seq(min(peaks.spat.final.fp$n.gps), max(peaks.spat.final.fp$n.gps),1),
                     n.ind=round(mean(peaks.spat.final.fp$n.ind)))

fp1pred <- predict(fp1, newdat, se=TRUE, type="response")

pred.spat.fp.gps <- data.frame(fit=fp1pred$fit,
                               se=fp1pred$se.fit) %>% 
  cbind(newdat) %>% 
  mutate(diff = fit-lag(fit),
         con = "spatial") 

#Plot predictions
ggplot(pred.spat.fp.ind) +
  geom_line(aes(x=n.ind, y=fit)) +
  geom_ribbon(aes(x=n.ind, ymax=fit+2*se, ymin=fit-2*se), alpha=0.5)

ggplot(pred.spat.fp.gps) +
  geom_line(aes(x=n.gps, y=fit)) +
  geom_ribbon(aes(x=n.gps, ymax=fit+2*se, ymin=fit-2*se), alpha=0.5)

  
##2e. Calculate mean peak value (true positives only)----
peaks.spat.final.tp <- peaks.spat.final.clean %>% 
  dplyr::filter(true==1,
                !is.na(Mmean)) %>% 
  mutate(CI95 = Mhigh95 - Mlow95)

summary(peaks.spat.final.tp)
sd(peaks.spat.final.tp$Mmean)

m1 <- lm(Mmean ~ n.ind+n.gps, data=peaks.spat.final.tp)
summary(m1)

ggplot(peaks.spat.final.tp) +
  geom_point(aes(x=n.ind, y=Mmean)) +
  geom_smooth(aes(x=n.ind, y=Mmean), method="lm")

ggplot(peaks.spat.final.tp) +
  geom_point(aes(x=n.gps, y=Mmean)) +
  geom_smooth(aes(x=n.gps, y=Mmean), method="lm")

##2f. Calcualte mean CI (true positives only)----
ci1 <- lm(CI95 ~ log(n.ind)+log(n.gps), data=peaks.spat.final.tp)
summary(ci1)

ggplot(peaks.spat.final.tp) +
  geom_point(aes(x=log(n.ind), y=CI95)) +
  geom_smooth(aes(x=log(n.ind), y=CI95))

ggplot(peaks.spat.final.tp) +
  geom_point(aes(x=log(n.gps), y=CI95)) +
  geom_smooth(aes(x=log(n.gps), y=CI95))

#2. TEMPORAL CONNECTIVITY####
##3a. Find peaks----
files.temp <- list.files(paste0(rd, "MCTemporal"))
peaks.temp.final <- data.frame()
for(i in 1:length(files.temp)){
  
  temp <- read.csv(paste0(rd, "MCTemporal/", files.temp[i])) %>% 
    dplyr::rename(M=mantel) %>% 
    dplyr::filter(!is.na(M))
  
  temp.sum <- temp %>% 
    group_by(doy, l1o) %>% 
    summarize(Mmean=mean(M),
              Mlow83=quantile(M, probs=0.083),
              Mhigh83=quantile(M, probs=0.917),
              Mlow95=quantile(M, probs=0.025),
              Mhigh95=quantile(M, probs=0.975),
              Mlow99=quantile(M, probs=0.005),
              Mhigh99=quantile(M, probs=0.995),
              Mlow=min(M),
              Mhigh=max(M),
              Mlowsd=mean(M)-sd(M),
              Mhighsd=mean(M)+sd(M)) %>% 
    ungroup()
  
  gamlist <- expand.grid(l1o=unique(temp$l1o),
                         k=c(5:15))
  
  pred.temp.all <- data.frame()
  for(j in 1:nrow(gamlist)){
    
    temp.j <- temp %>% 
      filter(l1o==gamlist$l1o[j])
    
    gam.j <- gam(M ~ s(doy, k = gamlist$k[j], bs = 'cs'), data=temp.j)
    
    pred.j <- data.frame(predict(gam.j, newdata=data.frame(doy=min(temp.j$doy):max(temp.j$doy)), se.fit=TRUE)) %>% 
      cbind(data.frame(doy=min(temp.j$doy):max(temp.j$doy))) %>% 
      mutate(upr = fit + (1.96*se.fit),
             lwr = fit - (1.96*se.fit),
             l1o = gamlist$l1o[j],
             k = gamlist$k[j])
    
    pred.temp.all <- rbind(pred.temp.all, pred.j)
    
  }
  
  ndoys <- table(temp.sum$l1o) %>% 
    data.frame() %>% 
    rename(l1o=Var1, lats=Freq) %>% 
    mutate(l1o=as.integer(as.character(l1o))) 
  
  pred.temp.k <- pred.temp.all %>% 
    inner_join(temp.sum, by=c("doy", "l1o")) %>% 
    mutate(use=ifelse(fit < Mhigh99 & fit > Mlow99, 1, 0)) %>% 
    group_by(l1o, k) %>% 
    summarize(sum=sum(use)) %>% 
    left_join(ndoys, by="l1o") %>% 
    dplyr::mutate(diff=lats-sum) %>% 
    dplyr::filter(diff==min(diff)) %>% 
    group_by(k) %>% 
    summarize(pops=n()) %>% 
    ungroup() %>% 
    dplyr::filter(pops==max(pops))
  
  temp.k <- min(pred.temp.k$k)
  
  print(paste0("The selected number of knots is ", temp.k))
  
  pred.temp <- pred.temp.all %>% 
    dplyr::filter(k==temp.k)
  
  ggplot() +
    geom_hex(aes(x=doy, y=M), data=temp) +
    geom_line(aes(x=doy, y=fit), data=pred.temp) +
    facet_wrap(~l1o)
  
  peaks.temp <- data.frame()
  for(j in 1:length(unique(pred.temp$l1o))){
    
    l1o.j <- unique(pred.temp$l1o)[j]
    
    pred.j <- pred.temp %>% 
      dplyr::filter(l1o==l1o.j)
    
    max.ind.j <- findpeaks(pred.j$fit, nups = 1, ndowns=1)
    min.ind.j <- findpeaks(-pred.j$fit, nups = 1, ndowns=1)
    
    max.j <- pred.j[max.ind.j[,2],] %>% 
      mutate(peak="max") %>% 
      data.frame() %>% 
      mutate(doylagdiff = abs(doy-lag(doy)),
             doyleaddiff = abs(doy-lead(doy)),
             use = case_when(doylagdiff < 10 & fit < lag(fit) ~ 0,
                             doyleaddiff < 10 & fit < lead(fit) ~ 0),
             use = ifelse(is.na(use), 1, use))
    
    min.j <- pred.j[min.ind.j[,2],] %>% 
      mutate(peak="min") %>% 
      data.frame() %>% 
      mutate(doylagdiff = NA,
             doyleaddiff = NA,
             use = NA)
    
    peaks.temp <- rbind(peaks.temp, max.j, min.j) %>% 
      arrange(l1o, doy)
    
  }
  
  peaks.temp <- peaks.temp %>% 
    group_by(l1o) %>% 
    mutate(doylagdiff=as.numeric(doylagdiff),
           doyleaddiff=as.numeric(doyleaddiff),
           use=as.numeric(use)) %>% 
    mutate(use=case_when(peak=="min" & lead(use)==1 & lag(use)==0 ~ 0,
                         peak=="min" & lead(use)==0 & lag(use)==1 ~ 0,
                         !is.na(use) ~ use),
           use = ifelse(is.na(use), 1, use)) %>% 
    dplyr::filter(use!=0) %>% 
    mutate(id=row_number()) %>% 
    ungroup()
  
  mat.spat1 <- matrix((peaks.temp %>% 
                         dplyr::filter(peak=="max"))$doy)
  mat.spat2 <- matrix((peaks.temp %>% 
                         dplyr::filter(peak=="max"))$fit)
  mat.temp <- cbind(mat.spat1, mat.spat2)
  
  temp.ms <- meanShift(mat.temp,
                       algorithm="LINEAR",
                       bandwidth=c(2,2))
  peaks.temp.ms <- peaks.temp %>% 
    dplyr::filter(peak=="max") %>% 
    mutate(cluster = temp.ms$assignment)
  
  if(nrow(peaks.temp.ms)>0){
    
    temp.clust <- table(peaks.temp.ms$l1o, peaks.temp.ms$cluster) %>% 
      data.frame() %>% 
      rename(l1o=Var1, cluster=Var2) %>% 
      group_by(cluster) %>% 
      summarize(count=sum(Freq)) %>% 
      ungroup() %>% 
      dplyr::filter(count >= length(unique(temp$l1o)))
    
    peaks.temp.clust <- peaks.temp.ms %>% 
      dplyr::filter(cluster %in% temp.clust$cluster)
    peaks.temp.clust
    
    peaks.temp.ci <- peaks.temp %>%
      left_join(temp.sum, by=c("doy", "l1o")) %>% 
      dplyr::filter(l1o==99) %>% 
      mutate(Mlagdiff = round((Mmean - lag(Mmean)), 1),
             Mleaddiff = round((Mmean - lead(Mmean)), 1),
             lagci = ifelse(Mlow83 < lag(Mhigh83), 0, 1),
             leadci = ifelse(Mhigh83 < lead(Mhigh83), 0, 1),
             lagci = ifelse(is.na(lagci), 99, lagci),
             leadci = ifelse(is.na(leadci), 99, leadci)) %>% 
      mutate(Mlagdiff = ifelse(is.na(Mlagdiff), 99, Mlagdiff),
             Mleaddiff= ifelse(is.na(Mleaddiff), 99, Mleaddiff)) %>% 
      dplyr::filter(peak=="max",
                    lagci!=0,
                    leadci!=0,
                    !(Mlagdiff <= 0.1 & Mleaddiff <= 0.1))
    peaks.temp.ci
    
    peaks.temp.all <- peaks.temp.clust %>% 
      inner_join(peaks.temp.ci, by = c("fit", "se.fit", "doy", "upr", "lwr", "l1o", "k", "peak", "id")) %>% 
      mutate(file=files.temp[i])  %>% 
      separate(file, into=c("metric", "connectivity", "n.pop", "n.ind", "n.gps", "sim", ".filetype"), remove=FALSE)
    
    if(nrow(peaks.temp.all) >0){
      
      raw <- read.csv(paste0("Data/StartingPoints/SimulationStartingPoints_",
                             unique(peaks.temp.all$n.pop),
                             "-",
                             unique(peaks.temp.all$n.ind),
                             "-",
                             unique(peaks.temp.all$n.gps),
                             "_",
                             unique(peaks.temp.all$sim),
                             ".csv")) %>% 
        mutate(yday = yday(day))
      
      stop2day <- raw %>% 
        dplyr::filter(season=="stop2") %>% 
        dplyr::select(pop, yday) %>% 
        unique()
      
      mean(stop2day$yday)
      peaks.temp.final.i <- peaks.temp.all %>% 
        mutate(true = ifelse(doy <= (mean(stop2day$yday+10)) & doy >= (mean(stop2day$yday-10)), 1, 0))
      
      peaks.temp.final <- rbind(peaks.temp.final, peaks.temp.final.i)
      
      
    }
      
  }
  
  else{
    next
  }
  
  print(paste0("Completed processing simulation ", i, " of ", length(files.temp)))
  
}

##3b. Wrangle----
peaks.temp.final.clean <- peaks.temp.final %>% 
  full_join(data.frame(file=files.temp)) %>% 
  separate(file, into=c("metric", "connectivity", "n.pop", "n.ind", "n.gps", "sim", ".filetype"), remove=FALSE) %>% 
  mutate(detect = ifelse(is.na(cluster), 0, 1),
         true = ifelse(is.na(true), 1, true)) %>% 
  dplyr::select(file, n.pop, n.ind, n.gps, sim, doy, Mmean, Mlow83, Mhigh83, Mlow95, Mhigh95, true, detect) %>% 
  mutate(status=case_when(true==1 & detect==1 ~ "tp",
                          true==1 & detect==0 ~ "fn",
                          true==0 & detect==1 ~ "fp",
                          true==0 & detect==0 ~ "tn"),
         sim = as.numeric(sim),
         n.pop = as.numeric(n.pop),
         n.ind = as.numeric(n.ind),
         n.gps = as.numeric(n.gps)) %>% 
  arrange(n.pop, n.ind, n.gps, sim)

write.csv(peaks.temp.final.clean, "CleanedTemporalPeaks_All.csv", row.names = FALSE)

##3c. Calculate probability of false negative----
peaks.temp.final.clean <- read.csv("CleanedTemporalPeaks_All.csv") %>% 
  dplyr::filter(n.gps <= 50)

#Data
peaks.temp.final.fn <- peaks.temp.final.clean %>% 
  dplyr::filter(true==1)

#Rate
table(peaks.temp.final.fn$detect)

#Visualize
ggplot(peaks.temp.final.fn) +
  geom_jitter(aes(x=n.ind, y=detect)) +
  geom_smooth(aes(x=n.ind, y=detect), method="lm")

ggplot(peaks.temp.final.fn) +
  geom_jitter(aes(x=n.gps, y=detect)) +
  geom_smooth(aes(x=n.gps, y=detect), method="lm")

#Model
fn1 <- glm(detect ~ n.ind+n.gps, data=peaks.temp.final.fn, family="binomial")
summary(fn1)

#Predict
#Individuals
newdat <- data.frame(n.ind=seq(min(peaks.temp.final.fn$n.ind), max(peaks.temp.final.fn$n.ind),1),
                     n.gps=round(mean(peaks.temp.final.fn$n.gps)))

fp1pred <- predict(fn1, newdat, se=TRUE, type="response")

pred.temp.fn.ind <- data.frame(fit=fp1pred$fit,
                               se=fp1pred$se.fit) %>% 
  cbind(newdat) %>% 
  mutate(diff = fit-lag(fit),
         con = "temporal")

#GPS points
newdat <- data.frame(n.gps=seq(min(peaks.temp.final.fn$n.gps), max(peaks.temp.final.fn$n.gps),1),
                     n.ind=round(mean(peaks.temp.final.fn$n.ind)))

fp1pred <- predict(fn1, newdat, se=TRUE, type="response")

pred.temp.fn.gps <- data.frame(fit=fp1pred$fit,
                               se=fp1pred$se.fit) %>% 
  cbind(newdat) %>% 
  mutate(diff = fit-lag(fit),
         con = "temporal")

#Summarize predictions
summary(pred.temp.fn.ind$diff)
summary(pred.temp.fn.gps$diff)

#Plot predictions
ggplot(pred.temp.fn.ind) +
  geom_line(aes(x=n.ind, y=fit)) +
  geom_ribbon(aes(x=n.ind, ymax=fit+2*se, ymin=fit-2*se), alpha=0.5)

ggplot(pred.temp.fn.gps) +
  geom_line(aes(x=n.gps, y=fit)) +
  geom_ribbon(aes(x=n.gps, ymax=fit+2*se, ymin=fit-2*se), alpha=0.5)

##3d. Calcualte probability of false positive----
#Data
peaks.temp.final.fp <- peaks.temp.final.clean %>% 
  dplyr::filter(status=="fp") %>% 
  dplyr::select(file) %>% 
  unique() %>% 
  mutate(detect=1) %>% 
  full_join(data.frame(file=files.temp)) %>% 
  separate(file, into=c("metric", "connectivity", "n.pop", "n.ind", "n.gps", "sim", ".filetype"), remove=FALSE) %>% 
  mutate(sim = as.numeric(sim),
         n.pop = as.numeric(n.pop),
         n.ind = as.numeric(n.ind),
         n.gps = as.numeric(n.gps),
         detect=ifelse(is.na(detect), 0, detect)) %>% 
  arrange(n.pop, n.ind, n.gps, sim) %>% 
  dplyr::filter(n.gps <= 50) %>% 
  unique() 

#Rate
table(peaks.temp.final.fp$detect)

#Visualize
ggplot(peaks.temp.final.fp) +
  geom_jitter(aes(x=n.ind, y=detect)) +
  geom_smooth(aes(x=n.ind, y=detect), method="lm")

ggplot(peaks.temp.final.fp) +
  geom_jitter(aes(x=n.gps, y=detect)) +
  geom_smooth(aes(x=n.gps, y=detect), method="lm")

#Model
fp1 <- glm(detect ~ n.ind+n.gps, data=peaks.temp.final.fp, family="binomial")
summary(fp1)

#Predict
#Individuals
newdat <- data.frame(n.ind=seq(min(peaks.temp.final.fp$n.ind), max(peaks.temp.final.fp$n.ind),1),
                     n.gps=round(mean(peaks.temp.final.fp$n.gps)))

fp1pred <- predict(fp1, newdat, se=TRUE, type="response")

pred.temp.fp.ind <- data.frame(fit=fp1pred$fit,
                               se=fp1pred$se.fit) %>% 
  cbind(newdat) %>% 
  mutate(diff = fit-lag(fit),
         con = "temporal")

#GPS points
newdat <- data.frame(n.gps=seq(min(peaks.temp.final.fp$n.gps), max(peaks.temp.final.fp$n.gps),1),
                     n.ind=round(mean(peaks.temp.final.fp$n.ind)))

fp1pred <- predict(fp1, newdat, se=TRUE, type="response")

pred.temp.fp.gps <- data.frame(fit=fp1pred$fit,
                               se=fp1pred$se.fit) %>% 
  cbind(newdat) %>% 
  mutate(diff = fit-lag(fit),
         con = "temporal")

#Plot predictions
ggplot(pred.temp.fp.ind) +
  geom_line(aes(x=n.ind, y=fit)) +
  geom_ribbon(aes(x=n.ind, ymax=fit+2*se, ymin=fit-2*se), alpha=0.5)

ggplot(pred.temp.fp.gps) +
  geom_line(aes(x=n.gps, y=fit)) +
  geom_ribbon(aes(x=n.gps, ymax=fit+2*se, ymin=fit-2*se), alpha=0.5)

##3e. Calculate mean peak value (true positives only)----
peaks.temp.final.tp <- peaks.temp.final.clean %>% 
  dplyr::filter(true==1,
                !is.na(Mmean)) %>% 
  mutate(CI95 = Mhigh95 - Mlow95)

summary(peaks.temp.final.tp)
sd(peaks.temp.final.tp$Mmean)

m1 <- lm(Mmean ~ n.ind+n.gps, data=peaks.temp.final.tp)
summary(m1)

ggplot(peaks.temp.final.tp) +
  geom_point(aes(x=n.ind, y=Mmean)) +
  geom_smooth(aes(x=n.ind, y=Mmean), method="lm")

ggplot(peaks.temp.final.tp) +
  geom_point(aes(x=n.gps, y=Mmean)) +
  geom_smooth(aes(x=n.gps, y=Mmean), method="lm")

##3f. Calcualte mean CI (true positives only)----
ci1 <- lm(CI95 ~ log(n.ind)+log(n.gps), data=peaks.temp.final.tp)
summary(ci1)
plot(ci1)

ggplot(peaks.temp.final.tp) +
  geom_point(aes(x=log(n.ind), y=CI95)) +
  geom_smooth(aes(x=log(n.ind), y=CI95))

ggplot(peaks.temp.final.tp) +
  geom_point(aes(x=log(n.gps), y=CI95)) +
  geom_smooth(aes(x=log(n.gps), y=CI95))

#4. SUMMARIZE####
#4a. Predictions for figure----
pred.ind.fn <- rbind(pred.spat.fn.ind, pred.temp.fn.ind)
pred.ind.fp <- rbind(pred.spat.fp.ind, pred.temp.fp.ind)
pred.gps.fn <- rbind(pred.spat.fn.gps, pred.temp.fn.gps)
pred.gps.fp <- rbind(pred.spat.fp.gps, pred.temp.fp.gps)

write.csv(pred.ind.fn, "SimulationResults_FN_Individual.csv", row.names = FALSE)
write.csv(pred.ind.fp, "SimulationResults_FP_Individual.csv", row.names = FALSE)
write.csv(pred.gps.fn, "SimulationResults_FN_GPS.csv", row.names = FALSE)
write.csv(pred.gps.fp, "SimulationResults_FP_GPS.csv", row.names = FALSE)
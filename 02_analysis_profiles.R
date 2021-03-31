#title: Construction of migratory connectivity profiles from Knight et al. 2021. Comprehensive estimation of spatial and temporal migratory connectivity across the annual cycle to direct conservation efforts. Ecography ECOG-05111
#author: Elly C. Knight
#date: Oct. 23, 2020

#1. PRELIMINARY####
#load packages
library(tidyverse)
library(mgcv)
library(sf)
library(pracma)
library(meanShiftR)

#limit display of scientific notation
options(scipen = 999)

#Each analysis is done separately to allow for visual assessment of results through the analysis process

#2. FALL SPATIAL#####

##2a. Read in data----
fall.spatial.count <- read.csv("Fall_Spatial_Latitudes.csv") %>% 
  rename(lat = latr) %>% 
  dplyr::select(lat, count.pop, count.ind)

fall.spatial <- read.csv("MC_Fall_Spatial_LeaveOneOut.csv") %>% 
  dplyr::filter(lat %in% fall.spatial.count$lat)

#Summarize across bootstraps
fall.spatial.sum <- fall.spatial %>% 
  group_by(lat, pop) %>% 
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
            Mhighsd=mean(M)+sd(M))

##2b. GAM----
###2bi. Run GAM across a range of k----
#set number of k to try and number of leave one out iteration
fall.pops <- expand.grid(pop = unique(fall.spatial$pop),
                           k = c(5:15))

#empty dataframe for results
pred.fall.spatial.all <- data.frame()

#loop
for(i in 1:nrow(fall.pops)){
  
  #subset to leave one out iteration
  fall.spatial.i <- fall.spatial %>% 
    filter(pop==fall.pops$pop[i])
  
  #gam
  gam.i <- gam(M ~ s(lat, k = fall.pops$k[i], bs = 'cs'), data=fall.spatial.i)
  
  #predict
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(lat=min(fall.spatial.i$lat):max(fall.spatial.i$lat)), se.fit=TRUE)) %>% 
    cbind(data.frame(lat=min(fall.spatial.i$lat):max(fall.spatial.i$lat))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           pop = fall.pops$pop[i],
           k = fall.pops$k[i])
  
  #save out results
  pred.fall.spatial.all <- rbind(pred.fall.spatial.all, pred.i)
  
  #report status
  print(paste0("Finished number ", i, " of ", nrow(fall.pops), " iterations"))
  
}

###2bii. Select min k for which predicted mean is within 99% CI for all pops----
#Count number of latitudes for each leave one out iteration
nlats <- table(fall.spatial.sum$pop) %>% 
  data.frame() %>% 
  rename(pop=Var1, lats=Freq) %>% 
  mutate(pop=as.integer(as.character(pop)))

#remove gam k values for which mean prediction is outside the 99% CI
pred.fall.spatial.k <- pred.fall.spatial.all %>% 
  inner_join(fall.spatial.sum) %>% 
  mutate(use=ifelse(fit < Mhigh99 & fit > Mlow99, 1, 0)) %>% 
  group_by(pop, k) %>% 
  summarize(sum=sum(use)) %>% 
  left_join(nlats) %>% 
  dplyr::mutate(diff=lats-sum) %>% 
  dplyr::filter(diff==min(diff)) %>% 
  group_by(k) %>% 
  summarize(pops=n()) %>% 
  ungroup() %>% 
  dplyr::filter(pops==max(pops))

#take the minimum of remaining k values
k.fall.spatial <- min(pred.fall.spatial.k$k)
k.fall.spatial

###2biii. Filter to just the selected predictions----
pred.fall.spatial <- pred.fall.spatial.all %>% 
  dplyr::filter(k==k.fall.spatial)

###2biv. Visualize----
plot.fall.sp <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlowsd, ymax=Mhighsd), alpha=0.4, data=fall.spatial.sum) +
  geom_line(aes(x=lat, y=Mmean), data=fall.spatial.sum, colour="red") +
  geom_line(aes(x=lat, y=fit), data=pred.fall.spatial, colour="green") +
  scale_x_reverse() +
  facet_wrap(.~pop)
plot.fall.sp

write.csv(pred.fall.spatial, "GAM_Fall_Spatial.csv", row.names = FALSE)

##2c. Find peaks----
#empty dataframe for results
peaks.fall.spatial <- data.frame()

###2ci. Loop through leave one out iterations----
for(i in 1:length(unique(fall.pops$pop))){
  
  #set leave one out population
  pop.i <- unique(fall.pops$pop)[i]
  
  #filter data
  pred.i <- pred.fall.spatial %>% 
    dplyr::filter(pop==pop.i)
  
  #find maximums
  max.ind.i <- findpeaks(pred.i$fit, nups = 1, ndowns=1)
  max.i <- pred.i[max.ind.i[,2],] %>% 
    mutate(peak="max") %>% 
    data.frame()
  
  #find minimums
  min.ind.i <- findpeaks(-pred.i$fit, nups = 1, ndowns=1)
  min.i <- pred.i[min.ind.i[,2],] %>% 
    mutate(peak="min") %>% 
    data.frame()
  
  #save out results
  peaks.fall.spatial <- rbind(peaks.fall.spatial, max.i, min.i) %>% 
    arrange(pop, lat)
  
}

peaks.fall.spatial <- peaks.fall.spatial %>% 
  mutate(id=row_number())

###2cii. Visualize----
peaks <- ggplot(peaks.fall.spatial) +
  geom_point(aes(x=lat, y=fit, colour=factor(pop))) +
  facet_wrap(~peak, nrow=2)

gams <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=fall.spatial.sum.99) +
  geom_line(aes(x=lat, y=fit, colour=factor(pop)), data=pred.fall.spatial) +
  geom_line(aes(x=lat, y=fit), colour="black", size=1.5, data=pred.fall.spatial.99)

gridExtra::grid.arrange(gams, peaks, nrow=2)

##2d. Remove clusters of peaks that aren't in all leave one out iterations----

###2di. Create matrix of peaks by latitude----
mat.fall1 <- matrix((peaks.fall.spatial %>% 
                         dplyr::filter(peak=="max"))$lat)
mat.fall2 <- matrix((peaks.fall.spatial %>% 
                         dplyr::filter(peak=="max"))$fit)
mat.fall <- cbind(mat.fall1, mat.fall2)

###2dii. Cluster with mean shift classification----
fall.spatial.ms <- meanShift(mat.fall,
                               algorithm="LINEAR",
                               bandwidth=c(2,2))
peaks.fall.spatial.ms <- peaks.fall.spatial %>% 
  dplyr::filter(peak=="max") %>% 
  mutate(cluster = fall.spatial.ms$assignment)

###2diii. Visualize----
ggplot(peaks.fall.spatial.ms) +
  geom_point(aes(x=lat, y=fit, colour=factor(cluster)))

###2div. Filter----
fall.spatial.clust <- table(peaks.fall.spatial.ms$pop, peaks.fall.spatial.ms$cluster) %>% 
  data.frame() %>% 
  rename(pop=Var1, cluster=Var2) %>% 
  group_by(cluster) %>% 
  summarize(count=sum(Freq)) %>% 
  ungroup() %>% 
  dplyr::filter(count >= length(unique(fall.spatial$pop)))

peaks.fall.spatial.clust <- peaks.fall.spatial.ms %>% 
  dplyr::filter(cluster %in% fall.spatial.clust$cluster)
peaks.fall.spatial.clust

##2e. Test remaining peaks with 83.4% CI----
peaks.fall.spatial.ci <- peaks.fall.spatial %>%
  left_join(fall.spatial.sum) %>% 
  dplyr::filter(pop==99) %>% 
  mutate(lagdiff = round((Mmean - lag(Mmean)), 1),
         leaddiff = round((Mmean - lead(Mmean)), 1),
         lagci = ifelse(Mlow83 < lag(Mhigh83), 0, 1),
         leadci = ifelse(Mhigh83 < lead(Mhigh83), 0, 1)) %>% 
  dplyr::filter(peak=="max",
                lagci==1 & leadci==1,
                !(lagdiff <= 0.1 & leaddiff <= 0.1))
peaks.fall.spatial.ci

##2f. Final list of peaks----
peaks.fall.spatial.final <- peaks.fall.spatial.clust %>% 
  inner_join(peaks.fall.spatial.ci)
peaks.fall.spatial.final

#3. SPRING SPATIAL#####

##3a. Read in data----
spring.spatial.count <- read.csv("Spring_Spatial_Latitudes.csv") %>% 
  rename(lat = latr) %>% 
  dplyr::select(lat, count.pop, count.ind)

spring.spatial <- read.csv("MC_Spring_Spatial_LeaveOneOut.csv") %>% 
  dplyr::filter(lat %in% spring.spatial.count$lat)

#Summarize across bootstraps
spring.spatial.sum <- spring.spatial %>% 
  group_by(lat, pop) %>% 
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
            Mhighsd=mean(M)+sd(M))

##3b. GAM----
###3bi. Run GAM across a range of k----
#set number of k to try and number of leave one out iteration
spring.pops <- expand.grid(pop = unique(spring.spatial$pop),
                         k = c(5:15))

#empty dataframe for results
pred.spring.spatial.all <- data.frame()

#loop
for(i in 1:nrow(spring.pops)){
  
  #subset to leave one out iteration
  spring.spatial.i <- spring.spatial %>% 
    filter(pop==spring.pops$pop[i])
  
  #gam
  gam.i <- gam(M ~ s(lat, k = spring.pops$k[i], bs = 'cs'), data=spring.spatial.i)
  
  #predict
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(lat=min(spring.spatial.i$lat):max(spring.spatial.i$lat)), se.fit=TRUE)) %>% 
    cbind(data.frame(lat=min(spring.spatial.i$lat):max(spring.spatial.i$lat))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           pop = spring.pops$pop[i],
           k = spring.pops$k[i])
  
  #save out results
  pred.spring.spatial.all <- rbind(pred.spring.spatial.all, pred.i)
  
  #report status
  print(paste0("Finished number ", i, " of ", nrow(spring.pops), " iterations"))
  
}

###3bii. Select min k for which predicted mean is within 99% CI for all pops----
#Count number of latitudes for each leave one out iteration
nlats <- table(spring.spatial.sum$pop) %>% 
  data.frame() %>% 
  rename(pop=Var1, lats=Freq) %>% 
  mutate(pop=as.integer(as.character(pop)))

#remove gam k values for which mean prediction is outside the 99% CI
pred.spring.spatial.k <- pred.spring.spatial.all %>% 
  inner_join(spring.spatial.sum) %>% 
  mutate(use=ifelse(fit < Mhigh99 & fit > Mlow99, 1, 0)) %>% 
  group_by(pop, k) %>% 
  summarize(sum=sum(use)) %>% 
  left_join(nlats) %>% 
  dplyr::mutate(diff=lats-sum) %>% 
  dplyr::filter(diff==min(diff)) %>% 
  group_by(k) %>% 
  summarize(pops=n()) %>% 
  ungroup() %>% 
  dplyr::filter(pops==max(pops))

#take the minimum of remaining k values
k.spring.spatial <- min(pred.spring.spatial.k$k)
k.spring.spatial

###3biii. Filter to just the selected predictions----
pred.spring.spatial <- pred.spring.spatial.all %>% 
  dplyr::filter(k==k.spring.spatial)

###3biv. Visualize----
plot.spring.sp <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlowsd, ymax=Mhighsd), alpha=0.4, data=spring.spatial.sum) +
  geom_line(aes(x=lat, y=Mmean), data=spring.spatial.sum, colour="red") +
  geom_line(aes(x=lat, y=fit), data=pred.spring.spatial, colour="green") +
  scale_x_reverse() +
  facet_wrap(.~pop)
plot.spring.sp

write.csv(pred.spring.spatial, "GAM_Spring_Spatial.csv", row.names = FALSE)

##3c. Find peaks----
#empty dataframe for results
peaks.spring.spatial <- data.frame()

###3ci. Loop through leave one out iterations----
for(i in 1:length(unique(spring.pops$pop))){
  
  #set leave one out population
  pop.i <- unique(spring.pops$pop)[i]
  
  #filter data
  pred.i <- pred.spring.spatial %>% 
    dplyr::filter(pop==pop.i)
  
  #find maximums
  max.ind.i <- findpeaks(pred.i$fit, nups = 1, ndowns=1)
  max.i <- pred.i[max.ind.i[,2],] %>% 
    mutate(peak="max") %>% 
    data.frame()
  
  #find minimums
  min.ind.i <- findpeaks(-pred.i$fit, nups = 1, ndowns=1)
  min.i <- pred.i[min.ind.i[,2],] %>% 
    mutate(peak="min") %>% 
    data.frame()
  
  #save out results
  peaks.spring.spatial <- rbind(peaks.spring.spatial, max.i, min.i) %>% 
    arrange(pop, lat)
  
}

peaks.spring.spatial <- peaks.spring.spatial %>% 
  mutate(id=row_number())

###3cii. Visualize----
peaks <- ggplot(peaks.spring.spatial) +
  geom_point(aes(x=lat, y=fit, colour=factor(pop))) +
  facet_wrap(~peak, nrow=2)

gams <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=spring.spatial.sum.99) +
  geom_line(aes(x=lat, y=fit, colour=factor(pop)), data=pred.spring.spatial) +
  geom_line(aes(x=lat, y=fit), colour="black", size=1.5, data=pred.spring.spatial.99)

gridExtra::grid.arrange(gams, peaks, nrow=2)

##3d. Remove clusters of peaks that aren't in all leave one out iterations----

###3di. Create matrix of peaks by latitude----
mat.spring1 <- matrix((peaks.spring.spatial %>% 
                       dplyr::filter(peak=="max"))$lat)
mat.spring2 <- matrix((peaks.spring.spatial %>% 
                       dplyr::filter(peak=="max"))$fit)
mat.spring <- cbind(mat.spring1, mat.spring2)

###3dii. Cluster with mean shift classification----
spring.spatial.ms <- meanShift(mat.spring,
                             algorithm="LINEAR",
                             bandwidth=c(2,2))
peaks.spring.spatial.ms <- peaks.spring.spatial %>% 
  dplyr::filter(peak=="max") %>% 
  mutate(cluster = spring.spatial.ms$assignment)

###3diii. Visualize----
ggplot(peaks.spring.spatial.ms) +
  geom_point(aes(x=lat, y=fit, colour=factor(cluster)))

###3div. Filter----
spring.spatial.clust <- table(peaks.spring.spatial.ms$pop, peaks.spring.spatial.ms$cluster) %>% 
  data.frame() %>% 
  rename(pop=Var1, cluster=Var2) %>% 
  group_by(cluster) %>% 
  summarize(count=sum(Freq)) %>% 
  ungroup() %>% 
  dplyr::filter(count >= length(unique(spring.spatial$pop)))

peaks.spring.spatial.clust <- peaks.spring.spatial.ms %>% 
  dplyr::filter(cluster %in% spring.spatial.clust$cluster)
peaks.spring.spatial.clust

##3e. Test remaining peaks with 83.4% CI----
peaks.spring.spatial.ci <- peaks.spring.spatial %>%
  left_join(spring.spatial.sum) %>% 
  dplyr::filter(pop==99) %>% 
  mutate(lagdiff = round((Mmean - lag(Mmean)), 1),
         leaddiff = round((Mmean - lead(Mmean)), 1),
         lagci = ifelse(Mlow83 < lag(Mhigh83), 0, 1),
         leadci = ifelse(Mhigh83 < lead(Mhigh83), 0, 1)) %>% 
  dplyr::filter(peak=="max",
                lagci==1 & leadci==1,
                !(lagdiff <= 0.1 & leaddiff <= 0.1))
peaks.spring.spatial.ci

##3f. Final list of peaks----
peaks.spring.spatial.final <- peaks.spring.spatial.clust %>% 
  inner_join(peaks.spring.spatial.ci)
peaks.spring.spatial.final

#4. FALL TEMPORAL#####

##4a. Read in data----
fall.temporal.count <- read.csv("Fall_Temporal_Days.csv") %>% 
  dplyr::select(doy, count.pop, count.ind)

fall.temporal <- read.csv("MC_Fall_Temporal.csv") %>% 
  dplyr::filter(doy %in% fall.temporal.count$doy)

#Summarize across bootstraps
fall.temporal.sum <- fall.temporal %>% 
  group_by(doy, day) %>% 
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
            Mhighsd=mean(M)+sd(M))

##4b. GAM----
##4bi. Run GAM across a range of k----
#set number of k to try and number of leave one out iteration
fall.pops <- expand.grid(day = unique(fall.temporal$day),
                         k = c(5:15))

#empty dataframe for results
pred.fall.temporal.all <- data.frame()

#loop
for(i in 1:nrow(fall.pops)){
  
  #subset to leave one out iteration
  fall.temporal.i <- fall.temporal %>% 
    filter(day==fall.pops$day[i])
  
  #gam
  gam.i <- gam(M ~ s(doy, k = fall.pops$k[i], bs = 'cs'), data=fall.temporal.i)
  
  #predict
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(doy=min(fall.temporal.i$doy):max(fall.temporal.i$doy)), se.fit=TRUE)) %>% 
    cbind(data.frame(doy=min(fall.temporal.i$doy):max(fall.temporal.i$doy))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           day = fall.pops$day[i],
           k = fall.pops$k[i])
  
  #save out results
  pred.fall.temporal.all <- rbind(pred.fall.temporal.all, pred.i)
  
  #report status
  print(paste0("Finished number ", i, " of ", nrow(fall.pops), " iterations"))
  
}

###4bii. Select min k for which predicted mean is within 99% CI for all pops----
#Count number of days for each leave one out iteration
ndays <- table(fall.temporal.sum$day) %>% 
  data.frame() %>% 
  rename(day=Var1, lats=Freq) %>% 
  mutate(day=as.integer(as.character(day)))

#remove gam k values for which mean prediction is outside the 99% CI
pred.fall.temporal.k <- pred.fall.temporal.all %>% 
  inner_join(fall.temporal.sum) %>% 
  mutate(use=ifelse(fit < Mhigh & fit > Mlow, 1, 0)) %>% 
  group_by(day, k) %>% 
  summarize(sum=sum(use)) %>% 
  left_join(ndays) %>% 
  dplyr::mutate(diff=lats-sum) %>% 
  dplyr::filter(diff==min(diff)) %>% 
  group_by(k) %>% 
  summarize(days=n()) %>% 
  ungroup() %>% 
  dplyr::filter(days==max(days))

#take the minimum of remaining k values
k.fall.temporal <- min(pred.fall.temporal.k$k)
k.fall.temporal

###4biii. Filter to just the selected predictions----
pred.fall.temporal <- pred.fall.temporal.all %>% 
  dplyr::filter(k==k.fall.temporal)

###4biv. Visualize----
plot.fall.sp <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow, ymax=Mhigh), alpha=0.4, data=fall.temporal.sum) +
  geom_line(aes(x=doy, y=fit), data=pred.fall.temporal, colour="green") +
  geom_line(aes(x=doy, y=Mmean), data=fall.temporal.sum, colour="red") +
  scale_x_reverse() +
  facet_wrap(.~day)
plot.fall.sp

write.csv(pred.fall.temporal, "FAM_Fall_Temporal.csv", row.names = FALSE)

##4c. Find peaks----
#empty dataframe for results
peaks.fall.temporal <- data.frame()

###4ci. Loop through leave one out iterations----
for(i in 1:length(unique(fall.pops$day))){
  
  #set leave one out day
  day.i <- unique(fall.pops$day)[i]
  
  #filter data
  pred.i <- pred.fall.temporal %>% 
    dplyr::filter(day==day.i)
  
  #find maximums
  max.ind.i <- findpeaks(pred.i$fit, nups = 1, ndowns=1)
  max.i <- pred.i[max.ind.i[,2],] %>% 
    mutate(peak="max") %>% 
    data.frame()
  
  #find minimums
  min.ind.i <- findpeaks(-pred.i$fit, nups = 1, ndowns=1)
  min.i <- pred.i[min.ind.i[,2],] %>% 
    mutate(peak="min") %>% 
    data.frame()
  
  #save out results
  peaks.fall.temporal <- rbind(peaks.fall.temporal, max.i, min.i) %>% 
    arrange(day, doy)
  
}

peaks.fall.temporal <- peaks.fall.temporal %>% 
  mutate(id=row_number())

###4cii. Visualize----
peaks <- ggplot(peaks.fall.temporal) +
  geom_point(aes(x=doy, y=fit, colour=factor(day))) +
  facet_wrap(~peak, nrow=2)

gams <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=fall.temporal.sum.99) +
  geom_line(aes(x=doy, y=fit, colour=factor(day)), data=pred.fall.temporal) +
  geom_line(aes(x=doy, y=fit), colour="black", size=1.5, data=pred.fall.temporal.99)

gridExtra::grid.arrange(gams, peaks, nrow=2)

##4d. Remove clusters of peaks that aren't in all leave one out iterations----

###4di. Create matrix of peaks by day----
mat.fall1 <- matrix((peaks.fall.temporal %>% 
                         dplyr::filter(peak=="max"))$doy)
mat.fall2 <- matrix((peaks.fall.temporal %>% 
                         dplyr::filter(peak=="max"))$fit)
mat.fall <- cbind(mat.fall1, mat.fall2)

###4dii. Cluster with mean shift classifcation----
fall.temporal.ms <- meanShift(mat.fall,
                             algorithm="LINEAR",
                             bandwidth=c(2,2))
peaks.fall.temporal.ms <- peaks.fall.temporal %>% 
  dplyr::filter(peak=="max") %>% 
  mutate(cluster = fall.temporal.ms$assignment)

###4diii. Visualize----
ggplot(peaks.fall.temporal.ms) +
  geom_point(aes(x=doy, y=fit, colour=factor(cluster)))

###4div. Filter----
fall.temporal.clust <- table(peaks.fall.temporal.ms$day, peaks.fall.temporal.ms$cluster) %>% 
  data.frame() %>% 
  rename(day=Var1, cluster=Var2) %>% 
  group_by(cluster) %>% 
  summarize(count=sum(Freq)) %>% 
  ungroup() %>% 
  dplyr::filter(count >= length(unique(fall.temporal$day)))

peaks.fall.temporal.clust <- peaks.fall.temporal.ms %>% 
  dplyr::filter(cluster %in% fall.temporal.clust$cluster)
peaks.fall.temporal.clust

##4e. Test remaining peaks with 83.4% CI----
peaks.fall.temporal.ci <- peaks.fall.temporal %>%
  left_join(fall.temporal.sum) %>% 
  dplyr::filter(day==99) %>% 
  mutate(lagdiff = round((Mmean - lag(Mmean)), 1),
         leaddiff = round((Mmean - lead(Mmean)), 1),
         lagci = ifelse(Mlow83 < lag(Mhigh83), 0, 1),
         leadci = ifelse(Mhigh83 < lead(Mhigh83), 0, 1)) %>% 
  dplyr::filter(peak=="max",
                lagci==1 & leadci==1,
                !(lagdiff <= 0.1 & leaddiff <= 0.1))
peaks.fall.temporal.ci

#4f. Final list of peaks----
peaks.fall.temporal.final <- peaks.fall.temporal.clust %>% 
  inner_join(peaks.fall.temporal.ci)
peaks.fall.temporal.final

##4g. Save out positions of peaks----
###4gi. Read in predictions from crawl----
pred.peaks.fall.temporal <- read.csv("CRAWL_Fall_Temporal_Predictions.csv") %>% 
  filter(doy %in% peaks.fall.temporal.final$doy) %>% 
  sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
  sf::st_set_crs(3857) %>% 
  sf::st_transform(4326)

###4gii. Save to day of peak----
pred.peaks.fall.temporal.xy <- pred.peaks.fall.temporal %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(data.frame(pred.peaks.fall.temporal)) %>% 
  dplyr::select(-geometry)

write.csv(pred.peaks.fall.temporal.xy, "Peaks_Fall_Temporal.csv", row.names = FALSE)

#5. SPRING TEMPORAL#####

##5a. Read in data----
spring.temporal.count <- read.csv("Spring_Temporal_Days.csv") %>% 
  dplyr::select(doy, count.pop, count.ind)

spring.temporal <- read.csv("MC_Spring_Temporal.csv") %>% 
  dplyr::filter(doy %in% spring.temporal.count$doy)

#Summarize across bootstraps
spring.temporal.sum <- spring.temporal %>% 
  group_by(doy, day) %>% 
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
            Mhighsd=mean(M)+sd(M))

##5b. GAM----
##5bi. Run GAM across a range of k----
#set number of k to try and number of leave one out iteration
spring.pops <- expand.grid(day = unique(spring.temporal$day),
                         k = c(5:15))

#empty dataframe for results
pred.spring.temporal.all <- data.frame()

#loop
for(i in 1:nrow(spring.pops)){
  
  #subset to leave one out iteration
  spring.temporal.i <- spring.temporal %>% 
    filter(day==spring.pops$day[i])
  
  #gam
  gam.i <- gam(M ~ s(doy, k = spring.pops$k[i], bs = 'cs'), data=spring.temporal.i)
  
  #predict
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(doy=min(spring.temporal.i$doy):max(spring.temporal.i$doy)), se.fit=TRUE)) %>% 
    cbind(data.frame(doy=min(spring.temporal.i$doy):max(spring.temporal.i$doy))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           day = spring.pops$day[i],
           k = spring.pops$k[i])
  
  #save out results
  pred.spring.temporal.all <- rbind(pred.spring.temporal.all, pred.i)
  
  #report status
  print(paste0("Finished number ", i, " of ", nrow(spring.pops), " iterations"))
  
}

###5bii. Select min k for which predicted mean is within 99% CI for all pops----
#Count number of days for each leave one out iteration
ndays <- table(spring.temporal.sum$day) %>% 
  data.frame() %>% 
  rename(day=Var1, lats=Freq) %>% 
  mutate(day=as.integer(as.character(day)))

#remove gam k values for which mean prediction is outside the 99% CI
pred.spring.temporal.k <- pred.spring.temporal.all %>% 
  inner_join(spring.temporal.sum) %>% 
  mutate(use=ifelse(fit < Mhigh & fit > Mlow, 1, 0)) %>% 
  group_by(day, k) %>% 
  summarize(sum=sum(use)) %>% 
  left_join(ndays) %>% 
  dplyr::mutate(diff=lats-sum) %>% 
  dplyr::filter(diff==min(diff)) %>% 
  group_by(k) %>% 
  summarize(days=n()) %>% 
  ungroup() %>% 
  dplyr::filter(days==max(days))

#take the minimum of remaining k values
k.spring.temporal <- min(pred.spring.temporal.k$k)
k.spring.temporal

###5biii. Filter to just the selected predictions----
pred.spring.temporal <- pred.spring.temporal.all %>% 
  dplyr::filter(k==k.spring.temporal)

###5biv. Visualize----
plot.spring.sp <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow, ymax=Mhigh), alpha=0.4, data=spring.temporal.sum) +
  geom_line(aes(x=doy, y=fit), data=pred.spring.temporal, colour="green") +
  geom_line(aes(x=doy, y=Mmean), data=spring.temporal.sum, colour="red") +
  scale_x_reverse() +
  facet_wrap(.~day)
plot.spring.sp

write.csv(pred.spring.temporal, "FAM_Spring_Temporal.csv", row.names = FALSE)

##5c. Find peaks----
#empty dataframe for results
peaks.spring.temporal <- data.frame()

###5ci. Loop through leave one out iterations----
for(i in 1:length(unique(spring.pops$day))){
  
  #set leave one out day
  day.i <- unique(spring.pops$day)[i]
  
  #filter data
  pred.i <- pred.spring.temporal %>% 
    dplyr::filter(day==day.i)
  
  #find maximums
  max.ind.i <- findpeaks(pred.i$fit, nups = 1, ndowns=1)
  max.i <- pred.i[max.ind.i[,2],] %>% 
    mutate(peak="max") %>% 
    data.frame()
  
  #find minimums
  min.ind.i <- findpeaks(-pred.i$fit, nups = 1, ndowns=1)
  min.i <- pred.i[min.ind.i[,2],] %>% 
    mutate(peak="min") %>% 
    data.frame()
  
  #save out results
  peaks.spring.temporal <- rbind(peaks.spring.temporal, max.i, min.i) %>% 
    arrange(day, doy)
  
}

peaks.spring.temporal <- peaks.spring.temporal %>% 
  mutate(id=row_number())

###5cii. Visualize----
peaks <- ggplot(peaks.spring.temporal) +
  geom_point(aes(x=doy, y=fit, colour=factor(day))) +
  facet_wrap(~peak, nrow=2)

gams <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=spring.temporal.sum.99) +
  geom_line(aes(x=doy, y=fit, colour=factor(day)), data=pred.spring.temporal) +
  geom_line(aes(x=doy, y=fit), colour="black", size=1.5, data=pred.spring.temporal.99)

gridExtra::grid.arrange(gams, peaks, nrow=2)

##5d. Remove clusters of peaks that aren't in all leave one out iterations----

###5di. Create matrix of peaks by day----
mat.spring1 <- matrix((peaks.spring.temporal %>% 
                       dplyr::filter(peak=="max"))$doy)
mat.spring2 <- matrix((peaks.spring.temporal %>% 
                       dplyr::filter(peak=="max"))$fit)
mat.spring <- cbind(mat.spring1, mat.spring2)

###5dii. Cluster with mean shift classifcation----
spring.temporal.ms <- meanShift(mat.spring,
                              algorithm="LINEAR",
                              bandwidth=c(2,2))
peaks.spring.temporal.ms <- peaks.spring.temporal %>% 
  dplyr::filter(peak=="max") %>% 
  mutate(cluster = spring.temporal.ms$assignment)

###5diii. Visualize----
ggplot(peaks.spring.temporal.ms) +
  geom_point(aes(x=doy, y=fit, colour=factor(cluster)))

###5div. Filter----
spring.temporal.clust <- table(peaks.spring.temporal.ms$day, peaks.spring.temporal.ms$cluster) %>% 
  data.frame() %>% 
  rename(day=Var1, cluster=Var2) %>% 
  group_by(cluster) %>% 
  summarize(count=sum(Freq)) %>% 
  ungroup() %>% 
  dplyr::filter(count >= length(unique(spring.temporal$day)))

peaks.spring.temporal.clust <- peaks.spring.temporal.ms %>% 
  dplyr::filter(cluster %in% spring.temporal.clust$cluster)
peaks.spring.temporal.clust

##5e. Test remaining peaks with 83.4% CI----
peaks.spring.temporal.ci <- peaks.spring.temporal %>%
  left_join(spring.temporal.sum) %>% 
  dplyr::filter(day==99) %>% 
  mutate(lagdiff = round((Mmean - lag(Mmean)), 1),
         leaddiff = round((Mmean - lead(Mmean)), 1),
         lagci = ifelse(Mlow83 < lag(Mhigh83), 0, 1),
         leadci = ifelse(Mhigh83 < lead(Mhigh83), 0, 1)) %>% 
  dplyr::filter(peak=="max",
                lagci==1 & leadci==1,
                !(lagdiff <= 0.1 & leaddiff <= 0.1))
peaks.spring.temporal.ci

#5f. Final list of peaks----
peaks.spring.temporal.final <- peaks.spring.temporal.clust %>% 
  inner_join(peaks.spring.temporal.ci)
peaks.spring.temporal.final

##5g. Save out positions of peaks----
###5gi. Read in predictions from crawl----
pred.peaks.spring.temporal <- read.csv("CRAWL_Spring_Temporal_Predictions.csv") %>% 
  filter(doy %in% peaks.spring.temporal.final$doy) %>% 
  sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
  sf::st_set_crs(3857) %>% 
  sf::st_transform(4326)

###5gii. Save to day of peak----
pred.peaks.spring.temporal.xy <- pred.peaks.spring.temporal %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(data.frame(pred.peaks.spring.temporal)) %>% 
  dplyr::select(-geometry)

write.csv(pred.peaks.spring.temporal.xy, "Peaks_Spring_Temporal.csv", row.names = FALSE)
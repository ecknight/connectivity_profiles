library(tidyverse)
library(mgcv)
library(sf)
library(pracma)
library(meanShiftR)

#Version 8: Selecting k for GAM automatically and using 83.4% CI on vs of mantel to select peaks

options(scipen = 999999)

count <- 3

#A. FALL SPATIAL#####

#1. Read in data----
fall.spatial.count <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_FallSpatial.csv") %>% 
  rename(lat = latr) %>% 
  dplyr::select(lat, count.pop, count.ind) %>% 
  dplyr::filter(count.ind >= count)

fall.spatial <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_Spatial_LeaveOneOut_MantelV4.csv") %>% 
  dplyr::filter(lat %in% fall.spatial.count$lat)

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

#2. GAM----

#2a. Run GAM across a range of k----
fall.pops <- expand.grid(pop = unique(fall.spatial$pop),
                           k = c(5:15))

pred.fall.spatial.all <- data.frame()
for(i in 1:nrow(fall.pops)){
  
  fall.spatial.i <- fall.spatial %>% 
    filter(pop==fall.pops$pop[i])
  
  gam.i <- gam(M ~ s(lat, k = fall.pops$k[i], bs = 'cs'), data=fall.spatial.i)
  
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(lat=min(fall.spatial.i$lat):max(fall.spatial.i$lat)), se.fit=TRUE)) %>% 
    cbind(data.frame(lat=min(fall.spatial.i$lat):max(fall.spatial.i$lat))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           pop = fall.pops$pop[i],
           k = fall.pops$k[i])
  
  pred.fall.spatial.all <- rbind(pred.fall.spatial.all, pred.i)
  
  print(paste0("Finished number ", i, " of ", nrow(fall.pops), " iterations"))
  
}

#2b. Select min k for which predicted mean is within 99% CI for all pops----
nlats <- table(fall.spatial.sum$pop) %>% 
  data.frame() %>% 
  rename(pop=Var1, lats=Freq) %>% 
  mutate(pop=as.integer(as.character(pop)))

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

k.fall.spatial <- min(pred.fall.spatial.k$k)
k.fall.spatial

#2c. Filter to just the selected predictions----

pred.fall.spatial <- pred.fall.spatial.all %>% 
  dplyr::filter(k==k.fall.spatial)

#2d. Visualize----
plot.fall.sp <- ggplot() +
  #  geom_hex(aes(x=lat, y=M), data=fall.spatial.ci) +
  geom_ribbon(aes(x=lat, ymin=Mlowsd, ymax=Mhighsd), alpha=0.4, data=fall.spatial.sum) +
  geom_line(aes(x=lat, y=Mmean), data=fall.spatial.sum, colour="red") +
  geom_line(aes(x=lat, y=fit), data=pred.fall.spatial, colour="green") +
  scale_x_reverse() +
  facet_wrap(.~pop)
plot.fall.sp

write.csv(pred.fall.spatial, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_Spatial_GAM_Mantel.csv", row.names = FALSE)


#3. Find peaks----
peaks.fall.spatial <- data.frame()
for(i in 1:length(unique(fall.pops$pop))){
  
  pop.i <- unique(fall.pops$pop)[i]
  
  pred.i <- pred.fall.spatial %>% 
    dplyr::filter(pop==pop.i)
  
  max.ind.i <- findpeaks(pred.i$fit, nups = 1, ndowns=1)
  min.ind.i <- findpeaks(-pred.i$fit, nups = 1, ndowns=1)
  
  max.i <- pred.i[max.ind.i[,2],] %>% 
    mutate(peak="max") %>% 
    data.frame()
  
  min.i <- pred.i[min.ind.i[,2],] %>% 
    mutate(peak="min") %>% 
    data.frame()
  
  peaks.fall.spatial <- rbind(peaks.fall.spatial, max.i, min.i) %>% 
    arrange(pop, lat)
  
}

peaks.fall.spatial <- peaks.fall.spatial %>% 
  mutate(id=row_number())

#4. Visualize----
pred.fall.spatial.99 <- pred.fall.spatial %>% 
  dplyr::filter(pop==99)

fall.spatial.99 <- fall.spatial %>% 
  dplyr::filter(pop==99)

fall.spatial.sum.99 <- fall.spatial.sum %>% 
  dplyr::filter(pop==99)

peaks <- ggplot(peaks.fall.spatial) +
  geom_point(aes(x=lat, y=fit, colour=factor(pop))) +
  facet_wrap(~peak, nrow=2)

gams <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=fall.spatial.sum.99) +
  geom_line(aes(x=lat, y=fit, colour=factor(pop)), data=pred.fall.spatial) +
  geom_line(aes(x=lat, y=fit), colour="black", size=1.5, data=pred.fall.spatial.99)

gridExtra::grid.arrange(gams, peaks, nrow=2)


#5. Remove clusters that aren't in all pops----

mat.fall1 <- matrix((peaks.fall.spatial %>% 
                         dplyr::filter(peak=="max"))$lat)
mat.fall2 <- matrix((peaks.fall.spatial %>% 
                         dplyr::filter(peak=="max"))$fit)
mat.fall <- cbind(mat.fall1, mat.fall2)

fall.spatial.ms <- meanShift(mat.fall,
                               algorithm="LINEAR",
                               bandwidth=c(2,2))
peaks.fall.spatial.ms <- peaks.fall.spatial %>% 
  dplyr::filter(peak=="max") %>% 
  mutate(cluster = fall.spatial.ms$assignment)

ggplot(peaks.fall.spatial.ms) +
  geom_point(aes(x=lat, y=fit, colour=factor(cluster)))

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

#6. Test remaining peaks with 83.4% CI----
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

#7. Final list of peaks----
peaks.fall.spatial.final <- peaks.fall.spatial.clust %>% 
  inner_join(peaks.fall.spatial.ci)
peaks.fall.spatial.final

#8. Fit final model----
fall.spatial.99 <- fall.spatial %>% 
  dplyr::filter(pop==99)

gam.fall.spatial <- gam(M ~ s(lat, bs = 'cs', k=k.fall.spatial), data=fall.spatial.99)
summary(gam.fall.spatial)

#B. SPRING SPATIAL#####

#1. Read in data----
spring.spatial.count <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_SpringSpatial.csv") %>% 
  rename(lat = latr) %>% 
  dplyr::select(lat, count.pop, count.ind) %>% 
  dplyr::filter(count.ind >= count)

spring.spatial <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Spatial_LeaveOneOut_MantelV4.csv") %>% 
  dplyr::filter(lat %in% spring.spatial.count$lat)

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

#2. GAM----

#2a. Run GAM across a range of k----
spring.pops <- expand.grid(pop = unique(spring.spatial$pop),
                           k = c(5:15))

pred.spring.spatial.all <- data.frame()
for(i in 1:nrow(spring.pops)){
  
  spring.spatial.i <- spring.spatial %>% 
    filter(pop==spring.pops$pop[i])
  
  gam.i <- gam(M ~ s(lat, k = spring.pops$k[i], bs = 'cs'), data=spring.spatial.i)
  
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(lat=min(spring.spatial.i$lat):max(spring.spatial.i$lat)), se.fit=TRUE)) %>% 
    cbind(data.frame(lat=min(spring.spatial.i$lat):max(spring.spatial.i$lat))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           pop = spring.pops$pop[i],
           k = spring.pops$k[i])
  
  pred.spring.spatial.all <- rbind(pred.spring.spatial.all, pred.i)
  
  print(paste0("Finished number ", i, " of ", nrow(spring.pops), " iterations"))
  
}

#2b. Select min k for which predicted mean is within 95% CI for all pops----
nlats <- table(spring.spatial.sum$pop) %>% 
  data.frame() %>% 
  rename(pop=Var1, lats=Freq) %>% 
  mutate(pop=as.integer(as.character(pop)))

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

k.spring.spatial <- min(pred.spring.spatial.k$k)
k.spring.spatial

#2c. Filter to just the selected predictions----

pred.spring.spatial <- pred.spring.spatial.all %>% 
  dplyr::filter(k==k.spring.spatial)

#2d. Visualize----
plot.spring.sp <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlow, ymax=Mhigh), alpha=0.4, data=spring.spatial.sum) +
  geom_line(aes(x=lat, y=Mmean), data=spring.spatial.sum, colour="red") +
  geom_line(aes(x=lat, y=fit), data=pred.spring.spatial, colour="green") +
  scale_x_reverse() +
  facet_wrap(.~pop)
plot.spring.sp

write.csv(pred.spring.spatial, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Spatial_GAM_Mantel.csv", row.names = FALSE)


#3. Find peaks----
peaks.spring.spatial <- data.frame()
for(i in 1:length(unique(spring.pops$pop))){
  
  pop.i <- unique(spring.pops$pop)[i]
  
  pred.i <- pred.spring.spatial %>% 
    dplyr::filter(pop==pop.i)
  
  max.ind.i <- findpeaks(pred.i$fit, nups = 1, ndowns=1)
  min.ind.i <- findpeaks(-pred.i$fit, nups = 1, ndowns=1)
  
  max.i <- pred.i[max.ind.i[,2],] %>% 
    mutate(peak="max") %>% 
    data.frame()
  
  min.i <- pred.i[min.ind.i[,2],] %>% 
    mutate(peak="min") %>% 
    data.frame()
  
  peaks.spring.spatial <- rbind(peaks.spring.spatial, max.i, min.i) %>% 
    arrange(pop, lat)
  
}

peaks.spring.spatial <- peaks.spring.spatial %>% 
  mutate(id=row_number())

#4. Visualize----
pred.spring.spatial.99 <- pred.spring.spatial %>% 
  dplyr::filter(pop==99)

spring.spatial.99 <- spring.spatial %>% 
  dplyr::filter(pop==99)

spring.spatial.sum.99 <- spring.spatial.sum %>% 
  dplyr::filter(pop==99)

peaks <- ggplot(peaks.spring.spatial) +
  geom_point(aes(x=lat, y=fit, colour=factor(pop))) +
  facet_wrap(~peak, nrow=2)

gams <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=spring.spatial.sum.99) +
  geom_line(aes(x=lat, y=fit, colour=factor(pop)), data=pred.spring.spatial) +
  geom_line(aes(x=lat, y=fit), colour="black", size=1.5, data=pred.spring.spatial.99)

gridExtra::grid.arrange(gams, peaks, nrow=2)


#5. Remove clusters that aren't in all pops----

mat.spring1 <- matrix((peaks.spring.spatial %>% 
                         dplyr::filter(peak=="max"))$lat)
mat.spring2 <- matrix((peaks.spring.spatial %>% 
                         dplyr::filter(peak=="max"))$fit)
mat.spring <- cbind(mat.spring1, mat.spring2)

spring.spatial.ms <- meanShift(mat.spring,
                               algorithm="LINEAR",
                               bandwidth=c(2,2))
peaks.spring.spatial.ms <- peaks.spring.spatial %>% 
  dplyr::filter(peak=="max") %>% 
  mutate(cluster = spring.spatial.ms$assignment)

ggplot(peaks.spring.spatial.ms) +
  geom_point(aes(x=lat, y=fit, colour=factor(cluster)))

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

#6. Test remaining peaks with 83.4% CI----
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

#7. Final list of peaks----
peaks.spring.spatial.final <- peaks.spring.spatial.clust %>% 
  inner_join(peaks.spring.spatial.ci)
peaks.spring.spatial.final

#8. Fit final model----
spring.spatial.99 <- spring.spatial %>% 
  dplyr::filter(pop==99)

gam.spring.spatial <- gam(M ~ s(lat, bs = 'cs', k=k.spring.spatial), data=spring.spatial.99)
summary(gam.spring.spatial)

#C. FALL TEMPORAL#####

#1. Read in data----
fall.temporal.count <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_FallTemporal.csv") %>% 
  dplyr::select(doy, count.pop, count.ind) %>% 
  dplyr::filter(count.ind >= count)

fall.temporal <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_Temporal_LeaveOneOut_MantelV4.csv") %>% 
  dplyr::filter(doy %in% fall.temporal.count$doy)

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

#2. GAM----

#2a. Run GAM across a range of k----
fall.pops <- expand.grid(day = unique(fall.temporal$day),
                         k = c(5:15))

pred.fall.temporal.all <- data.frame()
for(i in 1:nrow(fall.pops)){
  
  fall.temporal.i <- fall.temporal %>% 
    filter(day==fall.pops$day[i])
  
  gam.i <- gam(M ~ s(doy, k = fall.pops$k[i], bs = 'cs'), data=fall.temporal.i)
  
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(doy=min(fall.temporal.i$doy):max(fall.temporal.i$doy)), se.fit=TRUE)) %>% 
    cbind(data.frame(doy=min(fall.temporal.i$doy):max(fall.temporal.i$doy))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           day = fall.pops$day[i],
           k = fall.pops$k[i])
  
  pred.fall.temporal.all <- rbind(pred.fall.temporal.all, pred.i)
  
  print(paste0("Finished number ", i, " of ", nrow(fall.pops), " iterations"))
  
}

#2b. Select min k for which predicted mean is within 95% CI for all pops----
nlats <- table(fall.temporal.sum$day) %>% 
  data.frame() %>% 
  rename(day=Var1, lats=Freq) %>% 
  mutate(day=as.integer(as.character(day)))

pred.fall.temporal.k <- pred.fall.temporal.all %>% 
  inner_join(fall.temporal.sum) %>% 
  mutate(use=ifelse(fit < Mhigh & fit > Mlow, 1, 0)) %>% 
  group_by(day, k) %>% 
  summarize(sum=sum(use)) %>% 
  left_join(nlats) %>% 
  dplyr::mutate(diff=lats-sum) %>% 
  dplyr::filter(diff==min(diff)) %>% 
  group_by(k) %>% 
  summarize(days=n()) %>% 
  ungroup() %>% 
  dplyr::filter(days==max(days))

k.fall.temporal <- min(pred.fall.temporal.k$k)
k.fall.temporal

#2c. Filter to just the selected predictions----

pred.fall.temporal <- pred.fall.temporal.all %>% 
  dplyr::filter(k==k.fall.temporal)

#2d. Visualize----
plot.fall.sp <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow, ymax=Mhigh), alpha=0.4, data=fall.temporal.sum) +
  geom_line(aes(x=doy, y=fit), data=pred.fall.temporal, colour="green") +
  geom_line(aes(x=doy, y=Mmean), data=fall.temporal.sum, colour="red") +
  scale_x_reverse() +
  facet_wrap(.~day)
plot.fall.sp

write.csv(pred.fall.temporal, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_Temporal_GAM_Mantel.csv", row.names = FALSE)


#3. Find peaks----
peaks.fall.temporal <- data.frame()
for(i in 1:length(unique(fall.pops$day))){
  
  day.i <- unique(fall.pops$day)[i]
  
  pred.i <- pred.fall.temporal %>% 
    dplyr::filter(day==day.i)
  
  max.ind.i <- findpeaks(pred.i$fit, nups = 1, ndowns=1)
  min.ind.i <- findpeaks(-pred.i$fit, nups = 1, ndowns=1)
  
  max.i <- pred.i[max.ind.i[,2],] %>% 
    mutate(peak="max") %>% 
    data.frame()
  
  min.i <- pred.i[min.ind.i[,2],] %>% 
    mutate(peak="min") %>% 
    data.frame()
  
  peaks.fall.temporal <- rbind(peaks.fall.temporal, max.i, min.i) %>% 
    arrange(day, doy)
  
}

peaks.fall.temporal <- peaks.fall.temporal %>% 
  mutate(id=row_number())

#4. Visualize----
pred.fall.temporal.99 <- pred.fall.temporal %>% 
  dplyr::filter(day==99)

fall.temporal.99 <- fall.temporal %>% 
  dplyr::filter(day==99)

fall.temporal.sum.99 <- fall.temporal.sum %>% 
  dplyr::filter(day==99)

peaks <- ggplot(peaks.fall.temporal) +
  geom_point(aes(x=doy, y=fit, colour=factor(day))) +
  facet_wrap(~peak, nrow=2)

gams <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=fall.temporal.sum.99) +
  geom_line(aes(x=doy, y=fit, colour=factor(day)), data=pred.fall.temporal) +
  geom_line(aes(x=doy, y=fit), colour="black", size=1.5, data=pred.fall.temporal.99)

gridExtra::grid.arrange(gams, peaks, nrow=2)


#5. Remove clusters that aren't in all pops----

mat.fall1 <- matrix((peaks.fall.temporal %>% 
                         dplyr::filter(peak=="max"))$doy)
mat.fall2 <- matrix((peaks.fall.temporal %>% 
                         dplyr::filter(peak=="max"))$fit)
mat.fall <- cbind(mat.fall1, mat.fall2)

fall.temporal.ms <- meanShift(mat.fall,
                             algorithm="LINEAR",
                             bandwidth=c(2,2))
peaks.fall.temporal.ms <- peaks.fall.temporal %>% 
  dplyr::filter(peak=="max") %>% 
  mutate(cluster = fall.temporal.ms$assignment)

ggplot(peaks.fall.temporal.ms) +
  geom_point(aes(x=doy, y=fit, colour=factor(cluster)))

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

#6. Test remaining peaks with 83.4% CI----
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

#7. Final list of peaks----
peaks.fall.temporal.final <- peaks.fall.temporal.clust %>% 
  inner_join(peaks.fall.temporal.ci)
peaks.fall.temporal.final

#8. Fit final model----
fall.temporal.99 <- fall.temporal %>% 
  dplyr::filter(day==99)

gam.fall.temporal <- gam(M ~ s(doy, bs = 'cs', k=k.fall.temporal), data=fall.temporal.99)
summary(gam.fall.temporal)

#9. Save out positions of peaks----
pred.peaks.fall.temporal <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MovementModelPredictions_FallTemporal.csv") %>% 
  filter(doy %in% peaks.fall.temporal.final$doy) %>% 
  sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
  sf::st_set_crs(3857) %>% 
  sf::st_transform(4326)

pred.peaks.fall.temporal.xy <- pred.peaks.fall.temporal %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(data.frame(pred.peaks.fall.temporal)) %>% 
  dplyr::select(-geometry)

write.csv(pred.peaks.fall.temporal.xy, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/FallTemporalConnectivityPeaks_Mantel.csv", row.names = FALSE)


#D. SPRING TEMPORAL#####

#1. Read in data----
spring.temporal.count <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_SpringTemporal.csv") %>% 
  dplyr::select(doy, count.pop, count.ind) %>% 
  dplyr::filter(count.ind >= count)

spring.temporal <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Temporal_LeaveOneOut_MantelV3.csv") %>% 
  dplyr::filter(doy %in% spring.temporal.count$doy)

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

#2. GAM----

#2a. Run GAM across a range of k----
spring.pops <- expand.grid(day = unique(spring.temporal$day),
                         k = c(5:15))

pred.spring.temporal.all <- data.frame()
for(i in 1:nrow(spring.pops)){
  
  spring.temporal.i <- spring.temporal %>% 
    filter(day==spring.pops$day[i])
  
  gam.i <- gam(M ~ s(doy, k = spring.pops$k[i], bs = 'cs'), data=spring.temporal.i)
  
  pred.i <- data.frame(predict(gam.i, newdata=data.frame(doy=min(spring.temporal.i$doy):max(spring.temporal.i$doy)), se.fit=TRUE)) %>% 
    cbind(data.frame(doy=min(spring.temporal.i$doy):max(spring.temporal.i$doy))) %>% 
    mutate(upr = fit + (1.96*se.fit),
           lwr = fit - (1.96*se.fit),
           day = spring.pops$day[i],
           k = spring.pops$k[i])
  
  pred.spring.temporal.all <- rbind(pred.spring.temporal.all, pred.i)
  
  print(paste0("Finished number ", i, " of ", nrow(spring.pops), " iterations"))
  
}

#2b. Select min k for which predicted mean is within 95% CI for all pops----
nlats <- table(spring.temporal.sum$day) %>% 
  data.frame() %>% 
  rename(day=Var1, lats=Freq) %>% 
  mutate(day=as.integer(as.character(day)))

pred.spring.temporal.k <- pred.spring.temporal.all %>% 
  inner_join(spring.temporal.sum) %>% 
  mutate(use=ifelse(fit < Mhigh & fit > Mlow, 1, 0)) %>% 
  group_by(day, k) %>% 
  summarize(sum=sum(use)) %>% 
  left_join(nlats) %>% 
  dplyr::mutate(diff=lats-sum) %>% 
  dplyr::filter(diff==min(diff)) %>% 
  group_by(k) %>% 
  summarize(days=n()) %>% 
  ungroup() %>% 
  dplyr::filter(days==max(days))

k.spring.temporal <- min(pred.spring.temporal.k$k)
k.spring.temporal

#2c. Filter to just the selected predictions----

pred.spring.temporal <- pred.spring.temporal.all %>% 
  dplyr::filter(k==k.spring.temporal)

#2d. Visualize----
plot.spring.sp <- ggplot() +
#  geom_hex(aes(x=doy, y=M), data=spring.temporal.ci) +
  geom_ribbon(aes(x=doy, ymin=Mlow, ymax=Mhigh), alpha=0.4, data=spring.temporal.sum) +
  geom_line(aes(x=doy, y=fit), data=pred.spring.temporal, colour="green") +
  geom_line(aes(x=doy, y=Mmean), data=spring.temporal.sum, colour="red") +
  scale_x_reverse() +
  facet_wrap(.~day)
plot.spring.sp

write.csv(pred.spring.temporal, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Temporal_GAM_Mantel.csv", row.names = FALSE)

#3. Find peaks----
peaks.spring.temporal <- data.frame()
for(i in 1:length(unique(spring.pops$day))){
  
  day.i <- unique(spring.pops$day)[i]
  
  pred.i <- pred.spring.temporal %>% 
    dplyr::filter(day==day.i)
  
  max.ind.i <- findpeaks(pred.i$fit, nups = 1, ndowns=1)
  min.ind.i <- findpeaks(-pred.i$fit, nups = 1, ndowns=1)
  
  max.i <- pred.i[max.ind.i[,2],] %>% 
    mutate(peak="max") %>% 
    data.frame()
  
  min.i <- pred.i[min.ind.i[,2],] %>% 
    mutate(peak="min") %>% 
    data.frame()
  
  peaks.spring.temporal <- rbind(peaks.spring.temporal, max.i, min.i) %>% 
    arrange(day, doy)
  
}

peaks.spring.temporal <- peaks.spring.temporal %>% 
  mutate(id=row_number())

#4. Visualize----
pred.spring.temporal.99 <- pred.spring.temporal %>% 
  dplyr::filter(day==99)

spring.temporal.99 <- spring.temporal %>% 
  dplyr::filter(day==99)

spring.temporal.sum.99 <- spring.temporal.sum %>% 
  dplyr::filter(day==99)

peaks <- ggplot(peaks.spring.temporal) +
  geom_point(aes(x=doy, y=fit, colour=factor(day))) +
  facet_wrap(~peak, nrow=2)

gams <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=spring.temporal.sum.99) +
  geom_line(aes(x=doy, y=fit, colour=factor(day)), data=pred.spring.temporal) +
  geom_line(aes(x=doy, y=fit), colour="black", size=1.5, data=pred.spring.temporal.99)

gridExtra::grid.arrange(gams, peaks, nrow=2)

#5. Remove clusters that aren't in all pops----

mat.spring1 <- matrix((peaks.spring.temporal %>% 
                         dplyr::filter(peak=="max"))$doy)
mat.spring2 <- matrix((peaks.spring.temporal %>% 
                         dplyr::filter(peak=="max"))$fit)
mat.spring <- cbind(mat.spring1, mat.spring2)

spring.temporal.ms <- meanShift(mat.spring,
                              algorithm="LINEAR",
                              bandwidth=c(2,2))
peaks.spring.temporal.ms <- peaks.spring.temporal %>% 
  dplyr::filter(peak=="max") %>% 
  mutate(cluster = spring.temporal.ms$assignment)

ggplot(peaks.spring.temporal.ms) +
  geom_point(aes(x=doy, y=fit, colour=factor(cluster)))

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

#6. Test remaining peaks with 83.4% CI----
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

#7. Final list of peaks----
peaks.spring.temporal.final <- peaks.spring.temporal.clust %>% 
  inner_join(peaks.spring.temporal.ci)
peaks.spring.temporal.final

#8. Fit final model----
spring.temporal.99 <- spring.temporal %>% 
  dplyr::filter(day==99)

gam.spring.temporal <- gam(M ~ s(doy, bs = 'cs', k=k.spring.temporal), data=spring.temporal.99)
summary(gam.spring.temporal)


#8. Save out positions of peaks----
pred.peaks.spring.temporal <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MovementModelPredictions_SpringTemporal.csv") %>% 
  filter(doy %in% peaks.spring.temporal.final$doy) %>% 
  sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
  sf::st_set_crs(3857) %>% 
  sf::st_transform(4326)

pred.peaks.spring.temporal.xy <- pred.peaks.spring.temporal %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(data.frame(pred.peaks.spring.temporal)) %>% 
  dplyr::select(-geometry) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326)

write.csv(pred.peaks.spring.temporal.xy, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/SpringTemporalConnectivityPeaks_Mantel.csv", row.names = FALSE)


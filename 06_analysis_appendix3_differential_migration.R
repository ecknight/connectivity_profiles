#title: Differential migration data screening for Appendix 3 from Knight et al. 2021. Comprehensive estimation of spatial and temporal migratory connectivity across the annual cycle to direct conservation efforts. Ecography ECOG-05111
#author: Elly C. Knight
#date: Nov. 6, 2020

library(tidyverse)
library(lubridate)
library(lme4)
library(AICcmodavg)
library(MuMIn)
library(gridExtra)
library(mgcv)

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        plot.title=element_text(size=14))

make_names <- function(x) {
  new_names <- make.names(colnames(x))
  new_names <- gsub("\\.", "_", new_names)
  new_names <- tolower(new_names)
  colnames(x) <- new_names
  x
}

std <- function(x) sd(x)/sqrt(length(x))

world <- ggplot() +
  borders("world", colour = "gray85", fill = "gray80") +
  theme_classic() +
  xlim(-150, -30) +
  ylim(-55, 72)

#Minimum # of points per individual
min <- 1

#1. Load data----
dat <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  mutate(DateTime = ymd_hms(DateTime),
         BandDate = ymd(BandDate),
         Year = year(DateTime),
         YearDep = year(BandDate),
         Season=Season2) %>% 
  mutate(Winter2 = case_when(PinpointID %in% c(81, 439, 443, 490, 825, 826, 828) & Season=="Winter2" ~ 2,
                             PinpointID %in% c(81, 439, 443, 490, 825, 826, 828) & Season=="Winter" ~ 3,
                             !PinpointID %in% c(81, 439, 443, 490, 825, 826, 828) & Season=="Winter" ~ 1,
                             is.na(Season) ~ 0))

pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/tbl_population_abundance.csv")

#2. Wrangle----
n <- table(dat$PinpointID, dat$Season) %>% 
  data.frame() %>% 
  dplyr::filter(Freq >= min) %>% 
  dplyr::rename(PinpointID=Var1, Season=Var2) %>% 
  mutate(PinpointID=as.numeric(as.character(PinpointID)),
         Season=as.character(Season))

dat.wint <- dat %>% 
  dplyr::filter(Season=="Winter")
#  group_by(Population/PinpointID, Sex) %>% 
#  summarize(Lat=mean(Lat),
#            Long=mean(Long)) %>% 
#  ungroup()

dat.wint2 <- dat %>% 
  dplyr::filter(Winter2 %in% c(1,2),
                Season2 %in% c("Winter", "Winter2"))
#  group_by(Population/PinpointID, Sex) %>% 
#  summarize(Lat=mean(Lat),
#            Long=mean(Long)) %>% 
#  ungroup()

dat.mig <- dat %>% 
  dplyr::filter(Season %in% c("SpringMig", "FallMig")) %>% 
  dplyr::mutate(LatR=round(Lat, -1)) %>% 
  inner_join(n)

dat.mig.fall <- dat.mig %>% 
  dplyr::filter(Season=="FallMig")
table(dat.mig.fall$PinpointID)
nrow(dat.mig.fall)

dat.mig.spring <- dat.mig %>% 
  dplyr::filter(Season=="SpringMig")
table(dat.mig.spring$PinpointID)
nrow(dat.mig.spring)


#3. Wintering ground analyses----
ggplot(dat.wint) +
  geom_violin(aes(x=Sex, y=Lat))
ggplot(dat.wint2) +
  geom_violin(aes(x=Sex, y=Lat))

ggplot(dat.wint) +
  geom_violin(aes(x=Sex, y=Long))
ggplot(dat.wint2) +
  geom_violin(aes(x=Sex, y=Long))


lm.wint.lat <- lmer(Lat ~ Sex + (1|Population/PinpointID), data=dat.wint, na.action="na.fail", REML=F)
dredge(lm.wint.lat)

lm.wint.long <- lmer(Long ~ Sex + (1|Population/PinpointID), data=dat.wint, na.action="na.fail", REML=F)
dredge(lm.wint.long)

lm.wint2.lat <- lmer(Lat ~ Sex + (1|Population/PinpointID), data=dat.wint2, na.action="na.fail", REML=F)
dredge(lm.wint2.lat)

lm.wint2.long <- lmer(Long ~ Sex + (1|Population/PinpointID), data=dat.wint2, na.action="na.fail", REML=F)
dredge(lm.wint2.long)

#4. Migration analyses: spatial----
plot1 <- ggplot(dat.mig) +
  geom_point(aes(x=Long, y=Lat, colour=Sex)) +
  geom_smooth(aes(x=Long, y=Lat, colour=Sex)) +
  facet_wrap(~Season)

plot2 <- ggplot(dat.mig) +
  geom_point(aes(x=Long, y=Lat, colour=factor(YearDep))) +
  geom_smooth(aes(x=Long, y=Lat, colour=factor(YearDep))) +
  facet_wrap(~Season)

grid.arrange(plot1, plot2)

#Fall
lm.long.mig.fall.shape1 <- lmer(Lat ~ poly(Long,1) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.long.mig.fall.shape2 <- lmer(Lat ~ poly(Long,2) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.long.mig.fall.shape3 <- lmer(Lat ~ poly(Long,3) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
aictab(list(lm.long.mig.fall.shape1, lm.long.mig.fall.shape2, lm.long.mig.fall.shape3), sort=F)

lm.long.mig.fall1 <- lmer(Lat ~ poly(Long,3)*Sex + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.long.mig.fall2 <- lmer(Lat ~ poly(Long,3) + Sex + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.long.mig.fall3 <- lmer(Lat ~ poly(Long,3) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
aictab(list(lm.long.mig.fall1, lm.long.mig.fall2, lm.long.mig.fall3), sort=FALSE)

#Spring
lm.long.mig.spring.shape1 <- lmer(Lat ~ Long + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.long.mig.spring.shape2 <- lmer(Lat ~ poly(Long,2) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.long.mig.spring.shape3 <- lmer(Lat ~ poly(Long,3) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
aictab(list(lm.long.mig.spring.shape1, lm.long.mig.spring.shape2, lm.long.mig.spring.shape3), sort=F)

lm.long.mig.spring1 <- lmer(Lat ~ poly(Long,3)*Sex + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.long.mig.spring2 <- lmer(Lat ~ poly(Long,3) + Sex + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.long.mig.spring3 <- lmer(Lat ~ poly(Long,3) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
aictab(list(lm.long.mig.spring1, lm.long.mig.spring2, lm.long.mig.spring3), sort=F)

#5. Migration analyses: temporal----
plot3 <- ggplot(dat.mig) +
  geom_point(aes(x=doy, y=Lat, colour=Sex)) +
  geom_smooth(aes(x=doy, y=Lat, colour=Sex)) +
  facet_wrap(~Season, scales="free_x")

plot4 <- ggplot(dat.mig) +
  geom_point(aes(x=doy, y=Lat, colour=factor(YearDep))) +
  geom_smooth(aes(x=doy, y=Lat, colour=factor(YearDep))) +
  facet_wrap(~Season, scales="free_x")

grid.arrange(plot3, plot4)

#Fall
lm.doy.mig.fall.shape1 <- lmer(Lat ~ poly(doy,1) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall.shape2 <- lmer(Lat ~ poly(doy,2) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall.shape3 <- lmer(Lat ~ poly(doy,3) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
aictab(list(lm.doy.mig.fall.shape1, lm.doy.mig.fall.shape2, lm.doy.mig.fall.shape3), sort=F)

lm.doy.mig.fall1 <- lmer(Lat ~ poly(doy,3)*Sex + poly(doy,3)*factor(YearDep) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall2 <- lmer(Lat ~ poly(doy,3)*Sex + factor(YearDep) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall3 <- lmer(Lat ~ Sex + poly(doy,3)*factor(YearDep) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall4 <- lmer(Lat ~ poly(doy,3)*Sex + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall5 <- lmer(Lat ~ poly(doy,3)*factor(YearDep) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall6 <- lmer(Lat ~ poly(doy,3) + Sex + factor(YearDep) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall7 <- lmer(Lat ~ poly(doy,3) + Sex + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall8 <- lmer(Lat ~ poly(doy,3) + factor(YearDep) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
lm.doy.mig.fall9 <- lmer(Lat ~ poly(doy,3) + (1|Population/PinpointID), data=dat.mig.fall, REML=F)
aictab(list(lm.doy.mig.fall1, lm.doy.mig.fall2, lm.doy.mig.fall3, lm.doy.mig.fall4, lm.doy.mig.fall5, lm.doy.mig.fall6, lm.doy.mig.fall7, lm.doy.mig.fall8, lm.doy.mig.fall9), sort=F)

plot(lm.doy.mig.fall9)

#Spring
lm.doy.mig.spring.shape1 <- lmer(doy ~ poly(Lat,1) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring.shape2 <- lmer(doy ~ poly(Lat,2) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring.shape3 <- lmer(doy ~ poly(Lat,3) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
aictab(list(lm.doy.mig.spring.shape1, lm.doy.mig.spring.shape2, lm.doy.mig.spring.shape3), sort=F)

lm.doy.mig.spring1 <- lmer(Lat ~ poly(doy,3)*Sex + poly(doy,3)*factor(YearDep) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring2 <- lmer(Lat ~ poly(doy,3)*Sex + factor(YearDep) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring3 <- lmer(Lat ~ Sex + poly(doy,3)*factor(YearDep) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring4 <- lmer(Lat ~ poly(doy,3)*Sex + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring5 <- lmer(Lat ~ poly(doy,3)*factor(YearDep) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring6 <- lmer(Lat ~ poly(doy,3) + Sex + factor(YearDep) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring7 <- lmer(Lat ~ poly(doy,3) + Sex + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring8 <- lmer(Lat ~ poly(doy,3) + factor(YearDep) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
lm.doy.mig.spring9 <- lmer(Lat ~ poly(doy,3) + (1|Population/PinpointID), data=dat.mig.spring, REML=F)
aictab(list(lm.doy.mig.spring1, lm.doy.mig.spring2, lm.doy.mig.spring3, lm.doy.mig.spring4, lm.doy.mig.spring5, lm.doy.mig.spring6, lm.doy.mig.spring7, lm.doy.mig.spring8, lm.doy.mig.spring9), sort=F)

plot(lm.doy.mig.spring9)

pred <- predict(lm.doy.mig.spring4, re.form=~0, se.fit=TRUE) %>% 
  cbind(dat.mig.spring) %>% 
  mutate(upr=fit+1.96*se.fit,
         lwr=fit-1.96*se.fit)

plot5 <- ggplot(pred) +
  geom_line(aes(x=doy, y=fit, colour=Sex)) +
  geom_ribbon(aes(x=doy, ymax=upr, ymin=lwr, colour=Sex), alpha=0.5) +
  geom_point(aes(x=doy, y=Lat, colour=factor(Population)))

plot6 <- ggplot(pred) +
  geom_line(aes(x=doy, y=fit, colour=Sex)) +
  geom_ribbon(aes(x=doy, ymax=upr, ymin=lwr, colour=Sex), alpha=0.5) +
  geom_point(aes(x=doy, y=Lat, colour=factor(Sex)))

grid.arrange(plot5, plot6)

write.csv(pred, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/springmigratepred.csv", row.names=FALSE)


#6. Number of females----
dat.wint %>% 
  dplyr::filter(Sex=="F") %>% 
  dplyr::select(PinpointID) %>% 
  unique() %>% 
  nrow()

dat.wint %>% 
  dplyr::filter(PinpointID %in% c(81, 439, 443, 490, 825, 826, 828)) %>% 
  dplyr::select(PinpointID, Sex) %>% 
  unique()
table(dat.wint$Sex)

dat.wint %>% 
  dplyr::filter(Sex=="F") %>% 
  summarize(mean=mean(Long), 
            sd=sd(Long),
            max=max(Long),
            min=min(Long))

dat.wint %>% 
  dplyr::filter(Sex=="M") %>% 
  summarize(mean=mean(Long), 
            sd=sd(Long),
            max=max(Long),
            min=min(Long))

dat.wint2 %>% 
  dplyr::filter(Sex=="F") %>% 
  summarize(mean=mean(Long), 
            sd=sd(Long),
            max=max(Long),
            min=min(Long))

dat.wint2 %>% 
  dplyr::filter(Sex=="M") %>% 
  summarize(mean=mean(Long), 
            sd=sd(Long),
            max=max(Long),
            min=min(Long))

dat.mig.fall %>% 
  dplyr::filter(Sex=="F") %>% 
  dplyr::select(PinpointID) %>% 
  unique() %>% 
  nrow()

dat.mig.spring %>% 
  dplyr::filter(Sex=="F") %>% 
  dplyr::select(Population, PinpointID) %>% 
  unique()

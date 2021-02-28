#title: Wintering ground migratory connectivity analysis from Knight et al. 2021. Comprehensive estimation of spatial and temporal migratory connectivity across the annual cycle to direct conservation efforts. Ecography ECOG-05111
#author: Elly C. Knight
#date: Oct. 20, 2020

library(sp)
library(sf)
library(rgdal)
library(raster)
library(MigConnectivity)
library(tidyverse)
library(ncf)
library(gridExtra)
library(crawl)
library(ggspatial)
library(lubridate)

options(scipen = 999)

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

#PRELIM: DEFINE NEW FUNCTIONS####

estMantel <- function(targetPoints,
                      originPoints,
                      targetErrorX,
                      targetErrorY,
                      nBoot = 1000,
                      alpha = 0.05,
                      resampleProjection=MigConnectivity::projections$EquidistConic){
  
  #Prepare data
  targetPoints <- sp::spTransform(targetPoints, sp::CRS(MigConnectivity::projections$EquidistConic))
  targetDist <- dist(targetPoints@coords)
  originPoints <- sp::spTransform(originPoints, sp::CRS(MigConnectivity::projections$EquidistConic))
  originDist <- dist(originPoints@coords)
  
  #Calculate mantel without error
  pointCorr <- ade4::mantel.rtest(targetDist, originDist, nrepet = 1)
  
  #Bootstrap with error
  corr <- rep(NA, nBoot)
  
  for(boot in 1:nBoot){
    
    #Randomly sample spatial error
    sample.error.x <- apply(array(targetErrorX), 1, function(x) rnorm(1, mean=0, sd=x))
    sample.error.y <- apply(array(targetErrorY), 1, function(x) rnorm(1, mean=0, sd=x))
    
    #Add to coordinates
    point.sample <- data.frame(targetPoints@coords) %>% 
      cbind(sample.error.x, sample.error.y)  %>% 
      mutate(newx1 = coords.x1 + sample.error.x,
             newx2 = coords.x2 + sample.error.y) %>% 
      dplyr::select(newx1, newx2) %>% 
      sp::SpatialPoints(proj4string = sp::CRS(MigConnectivity::projections$EquidistConic))
    #Calculate distance matrix
    sampleDist <- dist(point.sample@coords)
    
    #Calculate mantel coefficient
    corr[boot] <- ade4::mantel.rtest(sampleDist, originDist, nrepet = 1)$obs
    
  }
  
  meanCorr <- mean(corr, na.rm = TRUE)
  medianCorr <- median(corr, na.rm = TRUE)
  seCorr <- sd(corr, na.rm = TRUE)
  simpleCICorr <- quantile(corr, c(alpha/2, 1 - alpha/2), na.rm = TRUE, 
                           type = 8, names = F)
  corr.z0 <- qnorm(sum((corr) < meanCorr)/nBoot)
  bcCICorr <- quantile(corr, pnorm(2 * corr.z0 + qnorm(c(alpha/2, 
                                                         1 - alpha/2))), na.rm = TRUE, type = 8, names = F)
  
  return(structure(list(sampleCorr = corr, pointCorr = pointCorr, 
                        meanCorr = meanCorr, medianCorr = medianCorr, seCorr = seCorr, 
                        simpleCICorr = simpleCICorr, bcCICorr = bcCICorr, alpha = alpha), 
                   class = c("estMantel")))
  
  
}

#1. LOAD DATA----
dat <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv")  %>% 
  mutate(Winter2 = case_when(PinpointID %in% c(81, 439, 443, 490, 825, 826, 828) & Season2=="Winter2" ~ 2,
                             PinpointID %in% c(81, 439, 443, 490, 825, 826, 828) & Season2=="Winter" ~ 3,
                             !PinpointID %in% c(81, 439, 443, 490, 825, 826, 828) & Season2=="Winter" ~ 1,
                             is.na(Season2) ~ 0)) %>% 
  dplyr::filter(PinpointID!=829)

pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/tbl_population_abundance.csv")

originids <- dat %>% 
  dplyr::filter(Winter==1) %>% 
  dplyr::select(PinpointID) %>% 
  unique()

targetids <- dat %>% 
  dplyr::filter(Winter==1,
                Season2=="Winter") %>% 
  dplyr::select(PinpointID) %>% 
  unique()


#2. Wrangle----
originPoints <- dat  %>% 
  dplyr::filter(Winter==1,
                Season2=="Breed1") %>% 
  dplyr::filter(!(Season2=="Breed1" & Type!="Band")) %>% 
  group_by(PinpointID) %>% 
  summarize(Lat=mean(Lat),
            Long=mean(Long)) %>% 
  ungroup() %>% 
  dplyr::select(Long, Lat) %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_geometry() %>% 
  as_Spatial()

#First wintering grounds----
targetPoints1 <- dat %>%
  dplyr::filter(Winter==1,
                Season2=="Winter") %>% 
  group_by(PinpointID) %>% 
  summarize(Lat=mean(Lat),
            Long=mean(Long)) %>% 
  ungroup() %>% 
  dplyr::select(Lat, Long) %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_geometry() %>% 
  as_Spatial()

targetError1 <- dat %>% 
  dplyr::filter(Winter==1,
                Season2=="Winter") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(dat %>% 
          dplyr::filter(Winter==1,
                        Season2=="Winter")) %>% 
  group_by(PinpointID) %>% 
  summarize(X=sd(X),
            Y=sd(Y)) %>% 
  ungroup()

targetError1X <- targetError1$X
targetError1Y <- targetError1$Y

#Second wintering grounds----
targetPoints2 <- dat %>%
  dplyr::filter(Winter2 %in% c(1,2),
                Season2 %in% c("Winter", "Winter2")) %>% 
  group_by(PinpointID) %>% 
  summarize(Lat=mean(Lat),
            Long=mean(Long)) %>% 
  ungroup() %>% 
  dplyr::select(Lat, Long) %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_geometry() %>% 
  as_Spatial()

targetError2 <- dat %>%
  dplyr::filter(Winter2 %in% c(1,2),
                Season2 %in% c("Winter", "Winter2")) %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(dat %>% 
          dplyr::filter(Winter2 %in% c(1,2),
                        Season2 %in% c("Winter", "Winter2"))) %>% 
#  dplyr::select(PinpointID, DateTime, Lat, Long, Season2, X, Y, WintDistRound) %>% 
  group_by(PinpointID) %>% 
  summarize(X=sd(X),
            Y=sd(Y)) %>% 
  ungroup()

targetError2X <- targetError2$X
targetError2Y <- targetError2$Y

#3. Calculate connectivity
Man1 <- estMantel(originPoints = originPoints,
                targetPoints = targetPoints1,
                targetErrorX = targetError1X,
                targetErrorY = targetError1Y,
                nBoot = 1000)

Man2 <- estMantel(originPoints = originPoints,
                  targetPoints = targetPoints2,
                  targetErrorX = targetError2X,
                  targetErrorY = targetError2Y,
                  nBoot = 1000)


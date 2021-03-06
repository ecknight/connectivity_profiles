#title: Spring migration spatial migratory connectivity analysis from Knight et al. 2021. Comprehensive estimation of spatial and temporal migratory connectivity across the annual cycle to direct conservation efforts. Ecography ECOG-05111
#author: Elly C. Knight
#date: Oct. 26, 2020

#1. PRELIMINARY####
#load packages
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
library(data.table)
library(ade4)

#limit display of scientific notation
options(scipen = 999)

#function for tidying dataframe names
make_names <- function(x) {
  new_names <- make.names(colnames(x))
  new_names <- gsub("\\.", "_", new_names)
  new_names <- tolower(new_names)
  colnames(x) <- new_names
  x
}

#write function for estimating connectivity
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

#Get map base data for visualization
world <- ggplot() +
  borders("world", colour = "gray85", fill = "gray80") +
  theme_classic() +
  xlim(-150, -30) +
  ylim(-55, 72)

#2. WRANGLE####
##2a. Load data----
tbl_locs_all <- read.csv("CONIMCP_CleanDataAll.csv") %>% 
  make_names()

pop <- read.csv("tbl_population_abundance.csv") %>% 
  make_names()

##2b. Select individuals for analysis----
###2bi. Try threshold for days & distance between points----
tbl_locs_spring <- tbl_locs_all %>% 
  dplyr::filter(season %in% c("SpringMig"))

ggplot(tbl_locs_spring, aes(x=duration, y=travdist)) +
  geom_point() +
  geom_smooth() 
#Relationship linear below huge values. Threshold won't work

###2bii. Use rules----
#1) a known wintering site, 2) at least one point in each of North, Central, and South America

tbl_IDs <- tbl_locs_all %>% 
  dplyr::filter(season=="SpringMig" |
                  (dplyr::lead(season)=="SpringMig" & season=="Winter") |
                  (dplyr::lag(season)=="SpringMig" & season=="Breed2") |
                  (type=="Breed1" & nbreed2==0 & nspring > 1)) %>% 
  dplyr::select(population, pinpointid, winter, contpoint) %>% 
  unique() %>% 
  dplyr::group_by(population, pinpointid, winter) %>% 
  summarize(cont=n()) %>% 
  mutate(cont=ifelse(pinpointid %in% c(61, 430, 433, 821, 826, 827), 3, cont)) %>% 
  dplyr::filter(winter==1,
                cont==3)

write.csv(tbl_IDs, "SpringMigration_Individuals.csv", row.names = FALSE)

#Number of birds selected per population
table(tbl_IDs$population)

IDs <- tbl_IDs %>% 
  dplyr::select(population, pinpointid) %>% 
  unique()

#3. PREPARE FOR MOVEMENT MODEL####
##3a: Tidy----
#filter to 1) selected individuals, 2) spring migration plus breeding and wintering sites, 3) points less than 1 day apart
#filter out banding site for the one that went 160 km NE from capture location
#Fix some date stuff
tbl_locs <- tbl_locs_all %>% 
  mutate(datetime = ymd_hms(datetime),
         fix=ifelse(type%in%c("Retrieval", "Band"), "3D", fix)) %>% 
  dplyr::arrange(pinpointid, datetime) %>% 
  dplyr::select(-date, -time) %>% 
  dplyr::filter(pinpointid %in% tbl_IDs$pinpointid,
                season=="SpringMig" |
                  (dplyr::lead(season)=="SpringMig" & season=="Winter") |
                  (dplyr::lag(season)=="SpringMig" & season=="Breed2") |
                  (type=="Band" & nbreed2==0 & nspring > 1) |
                  (type=="Breed1" & nbreed2==0 & nspring > 1),
                !(season=="SpringMig" & duration < 1),
                duration >= 1,
                !(pinpointid==828 & type=="Band")) %>% 
  dplyr::mutate(datetime=ifelse(season=="Breed1", ymd_hms("2018-06-12 12:00:00"), datetime),
                datetime=ifelse((pinpointid==441 & season=="Winter"), ymd_hms("2018-03-20 12:00:00"), datetime))

#Check seasons & transmission type for each ID
table(tbl_locs$pinpointid, tbl_locs$season)
table(tbl_locs$type)
table(tbl_locs$crc)

##3b. Define error parameters for GPS locations-----
tbl_locs <- tbl_locs %>% 
  dplyr::mutate(error_semi_major_axis=case_when(fix=="3D" ~ 20,
                                                fix=="2D" ~ 100),
                error_semi_minor_axis=case_when(fix=="3D" ~ 20,
                                                fix=="2D" ~ 100),
                error_ellipse_orientation=0) %>% 
  dplyr::rename(quality=fix)

tbl_locs <- crawl::argosDiag2Cov(tbl_locs$error_semi_major_axis, tbl_locs$error_semi_minor_axis, tbl_locs$error_ellipse_orientation) %>% 
  cbind(tbl_locs)

#save out for later
write.csv(tbl_locs, "CRAWL_Spring_Spatial_Points.csv", row.names = FALSE)

#save out as spatial object
sf_locs_wgs <- sf::st_as_sf(tbl_locs, coords = c("long","lat")) %>% 
  sf::st_set_crs(4326)
st_write(sf_locs_wgs, "CRAWL_Spring_Spatial_Points.shp", delete_dsn=TRUE)

##3c. Split data into a list with each individual as a separate item----
list_locs <- split(tbl_locs, tbl_locs$pinpointid)

##3d. Create a spatial object----
sf_locs_wgs <- lapply(list_locs,
                      function(x){SpatialPoints(
                        cbind(x$long,
                              x$lat), 
                        proj4string = CRS("+init=epsg:4326"))})

#transform to spherical mercator
sf_locs <- mapply(x=sf_locs_wgs,
                  y=list_locs,
                  function(x,y){
                    SpatialPointsDataFrame(spTransform(x,CRS("+init=epsg:3857")),y)
                  })

##3e. Create starting values for model----
inits <- lapply(sf_locs,
                function(x){
                  list(a = c(sp::coordinates(x)[1,1],0,
                             sp::coordinates(x)[1,2],0),
                       P = diag(c(10 ^ 2, 10 ^ 2, 
                                  10 ^ 2, 10 ^ 2)))})

##3f. Define fixed parameters for model----
fixpar <- c(1,1,NA,NA)

#4. FIT MOVEMENT MODEL####
#set seed for reproducability
set.seed(123)

##4a. Fit model----
fit <- mapply(x = sf_locs,
              y = inits,
              function(x,y){
                crawl::crwMLE(mov.model =  ~ 1,
                              err.model = list(
                                x =  ~ ln.sd.x - 1, 
                                y =  ~ ln.sd.y - 1, 
                                rho =  ~ error.corr),
                              data = x,
                              Time.name = "datetime",
                              time.scale = "days",
                              initial.state = y,
                              fixPar = fixpar,
                              initialSANN = list(maxit = 2500),
                              need.hess = TRUE,
                              theta=c(14,-1),
                              control = list(REPORT = 1000, trace = 6),
                              attempts = 8)},
              SIMPLIFY = FALSE) 

##4b. extract model parameters----
params <- lapply(fit,
                 function(x){
                   crawl::tidy_crwFit(x)
                 })

params.summary <- bind_rows(params, .id = "column_label")

#look at sigma
params.summary %>% 
  dplyr::filter(term=="ln sigma (Intercept)") %>% 
  summarize(est=mean(estimate),
            lwr=mean(conf.low),
            upr=mean(conf.high))

#look at beta
params.summary %>% 
  dplyr::filter(term=="ln beta (Intercept)") %>% 
  summarize(est=mean(estimate),
            lwr=mean(conf.low),
            upr=mean(conf.high))

#5. PREDICT FROM MOVEMENT MODEL####
##5a. times to predict for each individual----
predictTimes <- lapply(fit,
                       function(x){
                         seq(lubridate::ceiling_date(min(as.POSIXct(x$data$datetime,tz ="GMT")),"hour"),
                             lubridate::floor_date(max(as.POSIXct(x$data$datetime,tz = "GMT")),"hour"),
                             "1 hour")})

##5b. predict----
PredData <- mapply(x = fit,
                   y = predictTimes,
                   function(x,y){
                     crawl::crwPredict(x,predTime = y,return.type = "flat")},
                   SIMPLIFY = FALSE)

#tidy and save
PredDataAll <- data.frame(do.call(rbind, PredData)) %>% 
  dplyr::mutate(up.x = mu.x+se.mu.x,
                lw.x = mu.x-se.mu.x,
                up.y = mu.y+se.mu.y,
                lw.y = mu.y-se.mu.y)

write.csv(PredDataAll, "CRAWL_Spring_Spatial_Predictions.csv", row.names = FALSE)

#look at mean error
PredDataAll %>% 
  summarize(mean.x = mean(se.mu.x)/1000,
            mean.y = mean(se.mu.y)/1000,
            med.x = median(se.mu.x)/1000,
            med.y = median(se.mu.y)/1000,
            se.x = sd(se.mu.x)/1000,
            se.y = sd(se.mu.y)/1000,
            sd.x = sd(se.mu.x)/1000,
            sd.y = sd(se.mu.y)/1000) %>% 
  mutate(upr.x = mean.x + 1.96*se.x,
         lwr.x = mean.x - 1.96*se.x)

##5c. Convert predictions to spatial features----
#mean predictions
sf_locs_pred <- sf::st_as_sf(PredDataAll, coords = c("mu.x","mu.y")) %>% 
  sf::st_set_crs(3857)

sf_lines_mean <- sf_locs_pred %>% 
  dplyr::arrange(pinpointid, datetime) %>% 
  sf::st_geometry() %>% 
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs_pred$pinpointid))) %>% 
  sf::st_cast("MULTILINESTRING") %>% 
  sf::st_sf(pinpointid = as.factor(unique(sf_locs_pred$pinpointid)),
            pop=as.factor(IDs$population)) %>% 
  sf::st_transform(4326)

st_write(sf_lines_mean, "CRAWL_Spring_Spatial_Mean.shp", delete_dsn=TRUE)

#upper se interval
sf_lines_up <- sf_locs_up %>% 
  dplyr::arrange(pinpointid, datetime) %>% 
  sf::st_geometry() %>% 
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs_pred$pinpointid))) %>% 
  sf::st_cast("MULTILINESTRING") %>% 
  sf::st_sf(pinpointid = as.factor(unique(sf_locs_pred$pinpointid)),
            pop=as.factor(IDs$population)) %>% 
  sf::st_transform(4326)

sf_locs_up <- sf::st_as_sf(PredDataAll, coords=c("up.x", "up.y")) %>% 
  sf::st_set_crs(3857)

st_write(sf_lines_up, "CRAWL_Spring_Spatial_Upper.shp", delete_dsn=TRUE)

#lower se interval
sf_lines_lw <- sf_locs_lw %>% 
  dplyr::arrange(pinpointid, datetime) %>% 
  sf::st_geometry() %>% 
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs_pred$pinpointid))) %>% 
  sf::st_cast("MULTILINESTRING") %>% 
  sf::st_sf(pinpointid = as.factor(unique(sf_locs_pred$pinpointid)),
            pop=as.factor(IDs$population)) %>% 
  sf::st_transform(4326)

sf_locs_lw <- sf::st_as_sf(PredDataAll, coords=c("lw.x", "lw.y")) %>% 
  sf::st_set_crs(3857)

st_write(sf_lines_lw, "CRAWL_Spring_Spatial_Lower.shp", delete_dsn=TRUE)

##5d. Plot predictions----

plot2 <- world + 
  layer_spatial(sf_lines_mean, size = 0.75, aes(color=factor(pop))) +
  layer_spatial(sf_lines_up, size = 0.5, colour="grey30") +
  layer_spatial(sf_lines_lw, size = 0.5, colour="grey30") +
  layer_spatial(sf_locs_wgs, size = 0.75, col="red") +
  facet_wrap(.~pinpointid)

ggsave("CRAWL_Spring_Spatial.jpeg", device="jpeg", width=20, height=20, units="in")

#6. WRANGLE FOR MIGRATORY CONNECTIVITY ESTIMATION####
##6a. Select prediction for each latitude----
PredData_Spring_Spat <- PredDataAll %>% 
  st_as_sf(coords=c("mu.x", "mu.y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  cbind(PredDataAll) %>% 
  dplyr::select(population, pinpointid, datetime, mu.x, mu.y, se.mu.x, se.mu.y, up.x, up.y, lw.x, lw.y, locType, bandlat, bandlong, X, Y) %>% 
  mutate(latr = round(Y),
         latdiff=abs(latr-Y)) %>% 
  group_by(pinpointid, latr) %>% 
  mutate(n=n(),
         mindiff=min(latdiff)) %>% 
  dplyr::filter(mindiff==latdiff) %>% 
  mutate(n=n()) %>% 
  ungroup()

#inspect data
table(PredData_Spring_Spat$latr, PredData_Spring_Spat$population)

##6b. Summarize number of individuals and populations included at each latitude----
#only latitudes with more than 2 individuals and more than 1 population
spring.lats <- PredData_Spring_Spat %>% 
  group_by(latr, population) %>% 
  summarize(count = n()) %>% 
  group_by(latr) %>% 
  summarize(count.pop=n(),
            count.ind=sum(count)) %>% 
  mutate(percent.pop=count.pop/max(count.pop),
         percent.ind=count.ind/max(count.ind)) %>% 
  filter(count.pop > 1,
         count.ind > 2)

write.csv(spring.lats, "Spring_Spatial_Latitudes.csv", row.names = FALSE)

#7. LOOP FOR LEAVE ONE OUT

##7a. List of populations for leave one out----
pops <- c(unique(PredData_Spring_Spat$population), 99)

#list to save out results
mantel.flat <- list()

for(h in 1:length(pops)){
  
  #population to exclude for this iteration  
  pop.h <- pops[h]
  
  ##7b. Exclude population from data----
  dat.h <- PredData_Spring_Spat %>% 
    dplyr::filter(population != pop.h)
  
  ##7c. Determine latitudes to run for this iteration (same rules as above)----
  spring.lats.h <- dat.h %>% 
    group_by(latr, population) %>% 
    summarize(count = n()) %>% 
    group_by(latr) %>% 
    summarize(count.pop=n(),
              count.ind=sum(count)) %>% 
    mutate(percent.pop=count.pop/max(count.pop),
           percent.ind=count.ind/max(count.ind)) %>% 
    filter(count.pop > 1,
           count.ind > 2)
  
  #8. LOOP FOR CONNECTIVITY PROFILE####
  
  #list to save out results
  mantel <- list()
  
  for(i in 1:nrow(spring.lats.h)){
    
    #latitude for this iteration
    lat.i <- spring.lats.h$latr[i]
    
    ##8a. Subset data to latitude----
    dat.i <- dat.h %>% 
      dplyr::filter(latr==lat.i)
    
    ##8b. Create spatial objects----
    #Convert origin points to spatial object
    originPoints <- dat.i %>% 
      dplyr::select(bandlong, bandlat) %>% 
      rename(Long=bandlong, Lat=bandlat) %>% 
      st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
      st_transform(crs=3857) %>% 
      st_geometry() %>% 
      as_Spatial()
    
    #Convert target points to spatial object
    targetPoints <- dat.i %>% 
      dplyr::select(mu.x, mu.y) %>% 
      rename(Long=mu.x, Lat=mu.y) %>% 
      st_as_sf(coords=c("Long", "Lat"), crs=3857) %>% 
      st_geometry() %>% 
      as_Spatial()
    
    ##8c. Define error for each point----
    targetError = dat.i$se.mu.x
    
    ##8d. Calculate connectivity----
    Man <- try(estMantel(originPoints = originPoints,
                         targetPoints = targetPoints,
                         targetError = targetError,
                         nBoot = 1000))
    
    #Save single latitude results to list
    if (class(Man)=="estMantel"){
      man <- Man[["sampleCorr"]]
      mantel.i <- data.frame(mantel=man) %>% 
        mutate(lat=lat.i,
               pop=pop.h)
      mantel[[i]] <- mantel.i
    }
    
    #Report status
    print(paste0("Completed ", i, " of ", nrow(spring.lats), " lats"))
    
  }
  
  #Save out leave one out iteration to list
  mantel.flat[[h]] <- rbindlist(mantel)
  
  #Report status
  print(paste0("COMPLETED ", h, " of ", length(pops), " POPULATIONS"))
  
}

#Compress list to dataframe
spring.spatial<- rbindlist(mantel.flat) %>%  
  dplyr::rename(M = mantel) %>% 
  dplyr::filter(!is.na(M))

write.csv(spring.spatial, "MC_Spring_Spatial_LeaveOneOut.csv", row.names = FALSE)
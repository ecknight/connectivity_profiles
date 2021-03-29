#title: Spring migration spatial migratory connectivity analysis from Knight et al. 2021. Comprehensive estimation of spatial and temporal migratory connectivity across the annual cycle to direct conservation efforts. Ecography ECOG-05111
#author: Elly C. Knight
#date: Oct. 26, 2020

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

options(scipen = 9999)

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

#PRELIM: DEFINE NEW FUNCTIONS####

estMantel <- function(targetPoints,
                      originPoints,
                      targetError,
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
    sample.error <- apply(array(targetError), 1, function(x) rnorm(1, mean=0, sd=x))
    
    #Add to coordinates
    point.sample <- data.frame(targetPoints@coords) %>% 
      cbind(sample.error)  %>% 
      mutate(newx1 = coords.x1 + sample.error,
             newx2 = coords.x2 + sample.error) %>% 
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




#1. Load data----
dat <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  mutate(Season=ifelse(Season2 %in% c("Winter2", "WinterMig"), "Winter", Season2))

pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/tbl_population_abundance.csv")

my_cols <- cols(
  Population = col_integer(),
  PinpointID = col_integer(),
  Type = col_character(),
  CRC = col_character(),
  Fix = col_character(),
  Year = col_integer(),
  Date = col_datetime(format = ""),
  Time = col_time(format = ""),
  DateTime = col_datetime(format = ""),
  Season = col_character(),
  Lat = col_double(),
  Long = col_double()
)

tbl_locs_all <- readr::read_csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv", col_types = my_cols) %>% 
  mutate(Season=Season2)

ggplot(tbl_locs_all) +
  geom_point(aes(x=doy.order, y=factor(PinpointID), colour=Season)) +
  ylab("Pinpoint ID") +
  xlab("Day of schedule from August 10") +
  facet_wrap(.~Population, scales="free_y", ncol=2, strip.position = "right")

#2. Select threshold for days & distance between points----
tbl_locs_rate <- tbl_locs_all %>% 
  dplyr::filter(Season %in% c("SpringMig"))

ggplot(tbl_locs_rate, aes(x=Duration, y=TravDist)) +
  geom_point() +
  geom_smooth() #Relationship linear

#3. Plot points for each individual to visualize----
tbl_locs_plot <- tbl_locs_all %>% 
  dplyr::filter(Winter==1) %>% 
  dplyr::filter(Season=="SpringMig" |
                  (dplyr::lead(Season)=="SpringMig" & Season=="Winter") |
                  (dplyr::lag(Season)=="SpringMig" & Season=="Breed2") |
#                  (Type=="Band" & nbreed2==0 & nspring > 1) |
                  (Type=="Breed1" & nbreed2==0 & nspring > 1)) %>% 
  dplyr::mutate(Duration=ifelse(Type=="Band",0,Duration),
                Duration=ifelse(lead(Season)=="Winter",0,Duration)) %>% 
  arrange(Population, PinpointID, Date)

plot <- world +
  geom_point(aes(x=Long, y=Lat), data=tbl_locs_plot) +
  geom_line(aes(x=Long, y=Lat, group=PinpointID, colour=factor(Population)), data=tbl_locs_plot) +
  facet_wrap(~PinpointID)

#ggsave(plot=plot, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs/SpringSpatialVisualization.jpeg", device="jpeg", width=20, height=20, units="in")

#4. Filter to selected individuals----
tbl_IDs <- tbl_locs_all %>% 
  dplyr::filter(Winter==1) %>% 
  dplyr::filter(Season=="SpringMig" |
                  (dplyr::lead(Season)=="SpringMig" & Season=="Winter") |
                  (dplyr::lag(Season)=="SpringMig" & Season=="Breed2") |
#                  (Type=="Band" & nbreed2==0 & nspring > 1) |
                  (Type=="Breed1" & nbreed2==0 & nspring > 1)) %>% 
  dplyr::mutate(Duration=ifelse(Type=="Band",0,Duration),
                Duration=ifelse(lead(Season)=="Winter",0,Duration)) %>% 
  dplyr::group_by(PinpointID, ContPoint) %>% 
  dplyr::summarize(maxdur=max(Duration),
                   maxdist=max(TravDist),
                   pop=mean(Population)) %>% 
  dplyr::group_by(PinpointID) %>% 
  dplyr::summarize(cont=n(),
                   maxdur=max(maxdur),
                   maxdist=max(maxdist),
                   pop=mean(pop)) %>% 
  dplyr::mutate(selectdur=ifelse(round(maxdur) <= 21, 1, 0),
                selectdist=ifelse(round(maxdist, -3) <= 4000, 1, 0),
                selectcont=ifelse(cont==3, 1, 0))


IDs <- tbl_IDs %>% 
  dplyr::filter((selectcont >= 1) |
                  (PinpointID %in% c(61, 430, 433, 490, 821, 826, 827))) %>% 
  arrange(PinpointID) %>% 
  data.frame()
nbirds <- nrow(IDs)
nbirds
table(IDs$pop)

#3. Tidy data----
tbl_locs <- tbl_locs_all %>% 
  make_names() %>% 
  dplyr::rename(date_time = datetime,
                deployid = pinpointid,
                latitude = lat,
                longitude = long) %>% 
  dplyr::mutate(fix=ifelse(type%in%c("Retrieval", "Band"), "3D", fix)) %>% 
  dplyr::select(-date, -time) %>% 
  dplyr::filter(deployid %in% IDs$PinpointID) %>% 
  dplyr::filter(season=="SpringMig" |
                  (dplyr::lead(season)=="SpringMig" & season=="Winter") |
                  (dplyr::lag(season)=="SpringMig" & season=="Breed2") |
                  (type=="Band" & nbreed2==0 & nspring > 1) |
                  (type=="Breed1" & nbreed2==0 & nspring > 1),
                !(season=="SpringMig" & duration < 1)) %>% 
  dplyr::mutate(date_time=ymd_hms(ifelse(season=="Breed1", "2018-06-12 12:00:00", as.character(date_time))),
                date_time=ymd_hms(ifelse((deployid==441 & season=="Winter"), "2018-03-20 12:00:00", as.character(date_time))))  %>% 
  dplyr::filter(!(deployid==828 & type=="Band")) %>% #filter out banding site for the one that went 160 km NE from capture location %>% 
  group_by(deployid) %>% 
  mutate(travdist = abs(wintdist-lag(wintdist))) %>% 
  ungroup() %>% 
  dplyr::arrange(deployid, date_time)

IDs <- tbl_locs %>% 
  dplyr::select(population, deployid) %>% 
  unique()

#4. Check seasons & transmission type for each ID----
table(tbl_locs$deployid, tbl_locs$season)
table(tbl_locs$type)
table(tbl_locs$crc)

#5. Error parameters for GPS locations-----
tbl_locs <- tbl_locs %>% 
  dplyr::mutate(error_semi_major_axis=case_when(fix=="3D" ~ 20,
                                                fix=="2D" ~ 100),
                error_semi_minor_axis=case_when(fix=="3D" ~ 20,
                                                fix=="2D" ~ 100),
                error_ellipse_orientation=0) %>% 
  dplyr::rename(quality=fix)

tbl_locs <- crawl::argosDiag2Cov(tbl_locs$error_semi_major_axis, tbl_locs$error_semi_minor_axis, tbl_locs$error_ellipse_orientation) %>% 
  cbind(tbl_locs)

#6. Split out into each individual bird----
list_locs <- split(tbl_locs, tbl_locs$deployid)

#7. Create a spatial object----
sf_locs_wgs <- tbl_locs %>% 
  dplyr::select(deployid, longitude, latitude) %>% 
  sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(4326)
sf_locs <- sf::st_transform(sf_locs_wgs, 3857)

sf_locs_wgs <- lapply(list_locs,
                      function(x){SpatialPoints(
                        cbind(x$longitude,
                              x$latitude), 
                        proj4string = CRS("+init=epsg:4326"))})

sf_locs <- mapply(x=sf_locs_wgs,
                  y=list_locs,
                  function(x,y){
                    SpatialPointsDataFrame(spTransform(x,CRS("+init=epsg:3857")),y)
                  })

#8. Create starting values----
inits <- lapply(sf_locs,
                function(x){
                  list(a = c(sp::coordinates(x)[1,1],0,
                             sp::coordinates(x)[1,2],0),
                       P = diag(c(10 ^ 2, 10 ^ 2, 
                                  10 ^ 2, 10 ^ 2)))})

#9. Define fixed parameters----
fixpar <- c(1,1,NA,NA)

#10. Fit with crawl----
set.seed(123)

fit <- mapply(x = sf_locs,
              y = inits,
              function(x,y){
                crawl::crwMLE(mov.model =  ~ 1,
                              err.model = list(
                                x =  ~ ln.sd.x - 1, 
                                y =  ~ ln.sd.y - 1, 
                                rho =  ~ error.corr),
                              data = x,
                              Time.name = "date_time",
                              time.scale = "days",
                              initial.state = y,
                              fixPar = fixpar,
                              initialSANN = list(maxit = 2500),
                              need.hess = TRUE,
                              theta=c(14,-1),
                              control = list(REPORT = 1000, trace = 6),
                              attempts = 8)},
              SIMPLIFY = FALSE)

params <- lapply(fit,
                 function(x){
                   crawl::tidy_crwFit(x)
                 })
params

#Parameter summaries
params.summary <- bind_rows(params, .id = "column_label")

sigma <- params.summary %>% 
  dplyr::filter(term=="ln sigma (Intercept)") %>% 
  summarize(est=mean(estimate),
            lwr=mean(conf.low),
            upr=mean(conf.high))
sigma

beta <- params.summary %>% 
  dplyr::filter(term=="ln beta (Intercept)") %>% 
  summarize(est=mean(estimate),
            lwr=mean(conf.low),
            upr=mean(conf.high))
beta

#11. Predict----
predictTimes <- lapply(fit,
                       function(x){
                         seq(lubridate::ceiling_date(min(as.POSIXct(x$data$date_time,tz ="GMT")),"hour"),
                             lubridate::floor_date(max(as.POSIXct(x$data$date_time,tz = "GMT")),"hour"),
                             "1 hour")})

PredData <- mapply(x = fit,
                   y = predictTimes,
                   function(x,y){
                     crawl::crwPredict(x,predTime = y,return.type = "flat")},
                   SIMPLIFY = FALSE)

PredDataAll <- data.frame(do.call(rbind, PredData)) %>% 
  dplyr::mutate(up.x = mu.x+se.mu.x,
                lw.x = mu.x-se.mu.x,
                up.y = mu.y+se.mu.y,
                lw.y = mu.y-se.mu.y)

#Mean error
error <- PredDataAll %>% 
  summarize(error.x = mean(se.mu.x)/1000,
            error.y = mean(se.mu.y)/1000,
            se.x = std(se.mu.x)/1000,
            se.y = std(se.mu.y)/1000,
            sd.x = sd(se.mu.x)/1000,
            sd.y = sd(se.mu.y)/1000) %>% 
  mutate(upr.x = error.x + 1.96*se.x,
         lwr.x = error.x - 1.96*se.x)
error

write.csv(PredDataAll, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectvitiy_Spring_Spatial_1s-NA.csv", row.names = FALSE)

#12. Convert to spatial features----
sf_locs_pred <- sf::st_as_sf(PredDataAll, coords = c("mu.x","mu.y")) %>% 
  sf::st_set_crs(3857)

sf_locs_up <- sf::st_as_sf(PredDataAll, coords=c("up.x", "up.y")) %>% 
  sf::st_set_crs(3857)

sf_locs_lw <- sf::st_as_sf(PredDataAll, coords=c("lw.x", "lw.y")) %>% 
  sf::st_set_crs(3857)

sf_lines_mean <- sf_locs_pred %>% 
  dplyr::arrange(deployid, date_time) %>% 
  sf::st_geometry() %>% 
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs_pred$deployid))) %>% 
  sf::st_cast("MULTILINESTRING") %>% 
  sf::st_sf(deployid = as.factor(unique(sf_locs_pred$deployid)),
            pop=as.factor(IDs$population)) %>% 
  sf::st_transform(4326)

sf_lines_up <- sf_locs_up %>% 
  dplyr::arrange(deployid, date_time) %>% 
  sf::st_geometry() %>% 
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs_pred$deployid))) %>% 
  sf::st_cast("MULTILINESTRING") %>% 
  sf::st_sf(deployid = as.factor(unique(sf_locs_pred$deployid)),
            pop=as.factor(IDs$population)) %>% 
  sf::st_transform(4326)

sf_lines_lw <- sf_locs_lw %>% 
  dplyr::arrange(deployid, date_time) %>% 
  sf::st_geometry() %>% 
  sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sf_locs_pred$deployid))) %>% 
  sf::st_cast("MULTILINESTRING") %>% 
  sf::st_sf(deployid = as.factor(unique(sf_locs_pred$deployid)),
            pop=as.factor(IDs$population)) %>% 
  sf::st_transform(4326)

sf_locs_wgs <- sf::st_as_sf(tbl_locs, coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(4326)

st_write(sf_lines_mean, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Mean.shp", delete_dsn=TRUE)
st_write(sf_lines_up, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Upper.shp", delete_dsn=TRUE)
st_write(sf_lines_lw, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Lower.shp", delete_dsn=TRUE)
st_write(sf_locs_wgs, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Points.shp", delete_dsn=TRUE)
write.csv(tbl_locs, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Points.csv", row.names = FALSE)

#13. Plot predicted means----
world <- ggplot() +
  borders("world", colour = "gray85", fill = "gray80") +
  theme_classic()

plot2 <- world + 
  layer_spatial(sf_lines_mean, size = 0.75, aes(color=factor(pop))) +
  layer_spatial(sf_lines_up, size = 0.5, colour="grey30") +
  layer_spatial(sf_lines_lw, size = 0.5, colour="grey30") +
  layer_spatial(sf_locs_wgs, size = 0.75, col="red") +
  facet_wrap(.~deployid) +
  xlim(-150, -30) +
  ylim(-55, 72)

ggsave(plot=plot2, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs/CRAWL_MigConnectvitiy_Spring_Short_Spatial_1s-NA.jpeg", device="jpeg", width=20, height=20, units="in")

#14. Wrangle data for mig connectivity----
PredData_Spring_Spat <- PredDataAll %>% 
  st_as_sf(coords=c("mu.x", "mu.y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  cbind(PredDataAll) %>% 
  dplyr::select(population, deployid, date_time, mu.x, mu.y, se.mu.x, se.mu.y, up.x, up.y, lw.x, lw.y, locType, bandlat, bandlong, X, Y) %>% 
  mutate(latr = round(Y),
         latdiff=abs(latr-Y)) %>% 
  group_by(deployid, latr) %>% 
  mutate(n=n(),
         mindiff=min(latdiff)) %>% 
  dplyr::filter(mindiff==latdiff) %>% 
  mutate(n=n()) %>% 
  ungroup()

table(PredData_Spring_Spat$latr, PredData_Spring_Spat$population)

spring.lats <- PredData_Spring_Spat %>% 
  group_by(latr, population) %>% 
  summarize(count = n()) %>% 
  group_by(latr) %>% 
  summarize(count.pop=n(),
            count.ind=sum(count)) %>% 
  mutate(percent.pop=count.pop/max(count.pop),
         percent.ind=count.ind/max(count.ind)) %>% 
  filter(count.pop > 1)

write.csv(spring.lats, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_SpringSpatial.csv", row.names = FALSE)


#15. Set up leave one out----

pops <- c(unique(PredData_Spring_Spat$population), 99)

mantel.flat <- list()
for(h in 1:length(pops)){
  
  pop.h <- pops[h]
  
  dat.h <- PredData_Spring_Spat %>% 
    dplyr::filter(population != pop.h)
  
  spring.lats <- dat.h %>% 
    group_by(latr, population) %>% 
    summarize(count = n()) %>% 
    group_by(latr) %>% 
    summarize(count.pop=n(),
              count.ind=sum(count)) %>% 
    mutate(percent.pop=count.pop/max(count.pop),
           percent.ind=count.ind/max(count.ind)) %>% 
    filter(count.ind > 2)

  #16. Full model connectivity loop----
  mantel <- list()
  for(i in 1:nrow(spring.lats)){
    
    #Set doy & subset data
    lat.i <- spring.lats$latr[i]
    dat.i <- dat.h %>% 
      dplyr::filter(latr==lat.i)
    
      #Convert origin points to sp
      originPoints <- dat.i %>% 
        dplyr::select(bandlong, bandlat) %>% 
        rename(Long=bandlong, Lat=bandlat) %>% 
        st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
        st_transform(crs=3857) %>% 
        st_geometry() %>% 
        as_Spatial()
      
      #Convert target points to sp
      targetPoints <- dat.i %>% 
        dplyr::select(mu.x, mu.y) %>% 
        rename(Long=mu.x, Lat=mu.y) %>% 
        st_as_sf(coords=c("Long", "Lat"), crs=3857) %>% 
        st_geometry() %>% 
        as_Spatial()
      
      #Define error
      targetError <- dat.i$se.mu.x
      
      #Calculate connectivity
      Man <- try(estMantel(originPoints = originPoints,
                           targetPoints = targetPoints,
                           targetError = targetError,
                           nBoot = 1000))
        
        if (class(Man)=="estMantel"){
          
          #Save results
          man <- Man[["sampleCorr"]]
          mantel.i <- data.frame(mantel=man) %>% 
            mutate(lat=lat.i,
                   pop=pop.h)
          mantel[[i]] <- mantel.i
        
      }
      
      #Report status
      print(paste0("Completed ", i, " of ", nrow(spring.lats), " latitudes"))
      
  }
  
  mantel.flat[[h]] <- rbindlist(mantel)
  
  print(paste0("COMPLETED ", h, " of ", length(pops), " POPULATIONS"))
  
}

spring.spatial <- rbindlist(mantel.flat) %>%  
  dplyr::rename(M = mantel) %>% 
  dplyr::filter(!is.na(M))

write.csv(spring.spatial, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Spatial_LeaveOneOut_MantelV4.csv", row.names = FALSE)

spring.spatial <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Spatial_LeaveOneOut_MantelV4.csv")

spring.spatial.sum <- spring.spatial %>% 
  group_by(lat, pop) %>% 
  summarize(Mmean=mean(M),
            Mlowq=quantile(M, probs=0.087),
            Mhighq=quantile(M, probs=0.913))

#18. Visualize final profile----
plot1 <- ggplot() +
  geom_hex(aes(x=lat, y=M), data=spring.spatial) +
  geom_smooth(aes(x=lat, y=M), data=spring.spatial) +
  scale_fill_viridis_c() +
  facet_wrap(~pop) +
  theme(legend.position= "bottom")
plot1

plot2 <- ggplot(spring.spatial.sum) +
  geom_ribbon(aes(x=lat, ymin=Mlowq, ymax=Mhighq), alpha=0.4) +
  geom_line(aes(x=lat, y=Mmean)) +
#  geom_smooth(aes(x=lat, y=M), data=spring.spatial) +
  facet_wrap(~pop) + 
  ylim(c(-0.25, 1)) +
  xlim(c(-20, 60))
plot2

grid.arrange(plot1, plot2)

save.image("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Spatial_MantelV3.R")
#load("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Spatial_MantelV3.R")
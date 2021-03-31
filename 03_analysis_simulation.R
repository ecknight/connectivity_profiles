#title: Simulation from Knight et al. 2021. Comprehensive estimation of spatial and temporal migratory connectivity across the annual cycle to direct conservation efforts. Ecography ECOG-05111
#author: Elly C. Knight
#date: Aug. 14, 2020

#1. Preliminary####
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
library(mgcv)
library(pracma)
library(meanShiftR)

#Change working directory to subfolder
setwd(paste0(getwd(),"/simulation"))

#limit display of scientific notation
options(scipen = 999)

#set seed for reproducability
set.seed(123)

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

#2. Spatial settings####

#Breeding - connectivity = 1
breed.min.lat <- 50
breed.max.lat <- 65
breed.min.long <- -123
breed.max.long <- -96

#Stopover1 - connectivity = 0
stop1.min.lat <- 40
stop1.max.lat <- 45
stop1.min.long <- -120
stop1.max.long <- -96

#Stopover2 - connectivity = 1
stop2.min.lat <- 30
stop2.max.lat <- 35
stop2.min.long <- -115
stop2.max.long <- -96

#Winter - connectivity = 0
winter.min.lat <- 20
winter.max.lat <- 25
winter.min.long <- -105
winter.max.long <- -98

#Buffer distance between breeding populations (metres)
buff.breed <- 500000
buff.stop2 <- 300000

#Size of populations (metres)
pop.size <- 50000

#3. Temporal settings####

mig.length.min <- 20
mig.length.max <- 30

mig.start.min <- 1
mig.start.max <- 20

buff.date <- 4

#4. Other settings####

#Relative abundance
abun.max <- 100
abun.min <- 10

#Number of MC & Mantel samples
set.mc <- 100

#Number of simulation runs
set.sim <- 2000

#5. Sample size settings####

#Number of populations
set.pop <- c(3) #Must be at least 3 to do leave one out

#Number of individuals per population
set.ind.min <- 2
set.ind.max <- 20

#Number of samples per individual
set.gps.min <- 50
set.gps.max <- 1000

#List of settings to run
settings <- expand.grid(n.pop=set.pop, n.ind=set.ind, n.gps=set.gps)

#6. Set up simulation loop----
times <- data.frame()
for(f in 1:set.sim){
  
  n.files <- 0
  
  while(n.files < set.sim){
    
    start.i <- Sys.time()
    
    files <- data.frame(file=list.files(paste0(rd,"/Figs/Profile")))
    
    n.files <- nrow(files)
    
    n.sim <- n.files + 1
    
    n.pop <- set.pop
    n.ind <- round(runif(min=set.ind.min, max=set.ind.max, n=1))
    n.gps <- round(runif(min=set.gps.min, max=set.gps.max, n=1))
    
    #7. Select breeding departure dates with high temporal connectivity----
    dates.breed <- data.frame()
    while (nrow(dates.breed) < n.pop) {
      
      dates.breed.i <- data.frame(breed = round(runif(min=mig.start.min, max=mig.start.max, 1))) %>% 
        rbind(dates.breed)
      
      if(nrow(dates.breed.i)==1){
        dates.breed <- dates.breed.i
        next
      }
      else{
        dates.breed.mat <- dates.breed.i %>% 
          dist() %>% 
          as.matrix %>% 
          reshape2::melt(varnames = c("pop1", "pop2")) %>% 
          dplyr::filter(pop1!=pop2)
        
        if(any(dates.breed.mat$value < buff.date)==TRUE){
          dates.breed <- dates.breed
        }
        else{
          dates.breed <- dates.breed.i
        }
      }
    }
    
    #8. Select  other dates----
    dates <- data.frame(expand.grid(pop=(1:n.pop), ind=1:n.ind)) %>% 
      mutate(breed = rep(dates.breed$breed, n.ind)) %>% 
      mutate(stop1 = breed + round(runif(min=mig.length.min, max=mig.length.max, n.pop*n.ind)),
             stop2 = breed + 2*rep(round(runif(min=mig.length.min, max=mig.length.max, 1)), n.ind*n.pop),
             winter = stop2 + round(runif(min=mig.length.min, max=mig.length.max, n.pop*n.ind)),
             breed = as.Date(breed, origin="2020-01-01"),
             stop1 = as.Date(stop1, origin="2020-01-01"),
             stop2 = as.Date(stop2, origin="2020-01-01"),
             winter = as.Date(winter, origin="2020-01-01")) %>% 
      tidyr::gather("season", "day", breed:winter)

    #9. Select breeding population centroids----
    breed.pop <- data.frame()
    while (nrow(breed.pop) < n.pop) {
      
      breed.pop.i <- data.frame(lat = runif(min=breed.min.lat, max=breed.max.lat, n=1),
                                long = runif(min=breed.min.long, max=breed.max.long, n=1)) %>% 
        rbind(breed.pop)
      
      if(nrow(breed.pop.i)==1){
        breed.pop <- rbind(breed.pop, breed.pop.i) %>% 
          unique()
        next
      }
      else{
        breed.pop.mat <- breed.pop.i %>% 
          st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
          st_transform(3857) %>% 
          st_coordinates() %>% 
          data.frame() %>% 
          distFromPos("plane")
        
        if(any(breed.pop.mat < buff.breed & breed.pop.mat > 0)==TRUE){
          breed.pop <- breed.pop
        }
        else{
          breed.pop <- rbind(breed.pop, breed.pop.i) %>% 
            unique()
        }
      }
    }
    
    breed.pop <- breed.pop %>% 
      st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
      st_transform(3857) %>% 
      st_coordinates() %>% 
      cbind(breed.pop) %>% 
      mutate(pop = row_number())
    
    #10. Select individual breeding locations----
    locs.breed <- data.frame()
    for(i in 1:nrow(breed.pop)){
      pop = rep(breed.pop$pop[i], n.ind)
      ind = c(1:n.ind)
      season = rep("breed", n.ind)
      Y = runif(min = breed.pop$Y[i] - pop.size, max = breed.pop$Y[i] + pop.size, n=n.ind)
      X = runif(min = breed.pop$X[i] - pop.size, max = breed.pop$X[i] + pop.size, n=n.ind)
      
      locs.breed <- rbind(locs.breed, data.frame(X, Y, pop, ind, season))
      
    }

    
    #11. Select individual stopover #1 locations ----
    locs.stop1 <- data.frame(expand.grid(pop=c(1:n.pop), ind=c(1:n.ind)) %>% 
                          mutate(lat = runif(min=stop1.min.lat, max=stop1.max.lat, n=n.pop*n.ind),
                                 long = runif(min=stop1.min.long, max=stop1.max.long, n=n.pop*n.ind))) %>% 
      st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
      st_transform(3857) %>% 
      st_coordinates() %>% 
      cbind(expand.grid(pop=c(1:n.pop), ind=c(1:n.ind))) %>% 
      mutate(season = "stop1")
    
    #12. Select individual stopover #2 centroids ----
    stop2.pop <- data.frame()
    while (nrow(stop2.pop) < n.pop) {
      
      stop2.pop.i <- data.frame(lat = runif(min=stop2.min.lat, max=stop2.max.lat, n=1),
                                long = runif(min=stop2.min.long, max=stop2.max.long, n=1)) %>% 
        rbind(stop2.pop)
      
      if(nrow(stop2.pop.i)==1){
        stop2.pop <- rbind(stop2.pop, stop2.pop.i) %>% 
          unique()
        next
      }
      else{
        stop2.pop.mat <- stop2.pop.i %>% 
          st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
          st_transform(3857) %>% 
          st_coordinates() %>% 
          data.frame() %>% 
          distFromPos("plane")
        
        if(any(stop2.pop.mat < buff.stop2 & stop2.pop.mat > 0)==TRUE){
          stop2.pop <- stop2.pop
        }
        else{
          stop2.pop <- rbind(stop2.pop, stop2.pop.i) %>% 
            unique()
        }
      }
    }
    
    stop2.pop <- stop2.pop %>% 
      st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
      st_transform(3857) %>% 
      st_coordinates() %>% 
      cbind(stop2.pop) %>% 
      mutate(pop = row_number())
    
    #13. Select individual stopover 2 locations----
    locs.stop2 <- data.frame()
    for(i in 1:nrow(stop2.pop)){
      pop = rep(stop2.pop$pop[i], n.ind)
      ind = c(1:n.ind)
      season = rep("stop2", n.ind)
      Y = runif(min = stop2.pop$Y[i] - pop.size, max = stop2.pop$Y[i] + pop.size, n=n.ind)
      X = runif(min = stop2.pop$X[i] - pop.size, max = stop2.pop$X[i] + pop.size, n=n.ind)
      
      locs.stop2 <- rbind(locs.stop2, data.frame(X, Y, pop, ind, season))
      
    }
    
    #14. Select individual wintering locations----
    locs.winter <- data.frame(expand.grid(pop=c(1:n.pop), ind=c(1:n.ind)) %>% 
                               mutate(lat = runif(min=winter.min.lat, max=winter.max.lat, n=n.pop*n.ind),
                                      long = runif(min=winter.min.long, max=winter.max.long, n=n.pop*n.ind))) %>% 
      st_as_sf(coords=c("long", "lat"), crs=4326) %>% 
      st_transform(3857) %>% 
      st_coordinates() %>% 
      cbind(expand.grid(pop=c(1:n.pop), ind=c(1:n.ind))) %>% 
      mutate(season = "winter")
    
    
    #15. Merge dates and locations----
    locs.start <- rbind(locs.breed, locs.stop1, locs.stop2, locs.winter) %>% 
      full_join(dates)
      
    #16. Error parameters for GPS locations-----
    locs.start <- locs.start %>% 
      mutate(ID = paste(pop, ind, sep="-"),
             day = as.POSIXct(day)) %>% 
      dplyr::mutate(error_semi_major_axis=20,
                    error_semi_minor_axis=20,
                    error_ellipse_orientation=0)
    
    locs.start <- crawl::argosDiag2Cov(locs.start$error_semi_major_axis, locs.start$error_semi_minor_axis, locs.start$error_ellipse_orientation) %>% 
      cbind(locs.start)
    
    #17. Save initial simulation locations----
    write.csv(locs.start, file=paste0(rd,"Data/StartingPoints/SimulationStartingPoints_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
    
    #18. Wrangle for crawl----
    locs.start.sf1 <- locs.start %>% 
      dplyr::select(ID, X, Y) %>% 
      sf::st_as_sf(coords = c("X","Y")) %>% 
      sf::st_set_crs(3857) 
    
    locs.list1 <- split(locs.start, locs.start$ID)
    
    locs.list.sf1 <- lapply(locs.list1,
                                function(x){SpatialPoints(
                                  cbind(x$X,
                                        x$Y), 
                                  proj4string = CRS("+init=epsg:3857"))})
    
    locs.list.sf1 <- mapply(x=locs.list.sf1,
                                 y=locs.list1,
                                 function(x,y){
                                   SpatialPointsDataFrame(x,y)
                                 })
    
    #19. Create starting values----
    inits1 <- lapply(locs.list.sf1,
                     function(x){
                       list(a = c(sp::coordinates(x)[1,1],0,
                                  sp::coordinates(x)[1,2],0),
                            P = diag(c(10 ^ 2, 10 ^ 2, 
                                       10 ^ 2, 10 ^ 2)))})
    
    fixpar <- c(1,1,NA,NA)
    
    #20. Fit with crawl----
    fit1 <- mapply(x = locs.list.sf1,
                   y = inits1,
                   function(x,y){
                     try(crawl::crwMLE(mov.model =  ~ 1,
                                   err.model = list(
                                     x =  ~ ln.sd.x - 1, 
                                     y =  ~ ln.sd.y - 1, 
                                     rho =  ~ error.corr),
                                   data = x,
                                   Time.name = "day",
                                   time.scale = "days",
                                   initial.state = y,
                                   fixPar = fixpar,
                                   initialSANN = list(maxit = 2500),
                                   need.hess = TRUE,
                                   theta=c(14,-1),
                                   control = list(REPORT = 1000, trace = 6),
                                   attempts = 8))},
                   SIMPLIFY = FALSE)
    
    #21. Save crawl parameters----
    params1 <- lapply(fit1,
                      function(x){
                        crawl::tidy_crwFit(x)
                      })
    params.summary1 <- bind_rows(params1, .id = "column_label") #SAVE THIS
    
    write.csv(params.summary1, file=paste0(rd,"CrawlParameters/PredictedPaths1/CrawlParameters1_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
    
    #22. Predict----
    times1 <- lapply(fit1,
                     function(x){
                       seq(lubridate::ceiling_date(min(as.POSIXct(x$data$day,tz ="GMT")),"day"),
                           lubridate::floor_date(max(as.POSIXct(x$data$day,tz = "GMT")),"day"),
                           "1 day")})
    
    pred.list1 <- mapply(x = fit1,
                         y = times1,
                         function(x,y){
                           crawl::crwPredict(x,predTime = y,return.type = "flat")},
                         SIMPLIFY = FALSE)
    
    pred1 <- data.frame(do.call(rbind, pred.list1)) %>% 
      dplyr::mutate(up.x = mu.x+se.mu.x,
                    lw.x = mu.x-se.mu.x,
                    up.y = mu.y+se.mu.y,
                    lw.y = mu.y-se.mu.y)
    
    #23. Check if any individuals didn't work----
    fit.na1 <- pred1 %>% 
      dplyr::filter(is.na(se.mu.x))
    
    param.na1 <- params.summary1 %>% 
      dplyr::filter((term=="ln beta (Intercept)" &
                    std.error=="NaN") |
                    (term=="ln sigma (Intercept)" &
                       std.error=="NaN") |
                    (term=="ln beta (Intercept)" &
                         is.na(std.error)) |
                    (term=="ln sigma (Intercept)" &
                       is.na(std.error)))
    
    if(nrow(fit.na1) > 0 | nrow(param.na1) > 0) {
      
      next
    }
    
    else{
      
      
      #24. Convert predictions spatial features----
      pred.locs.sf1 <- sf::st_as_sf(pred1, coords = c("mu.x","mu.y")) %>% 
        sf::st_set_crs(3857)
      
      pred.lines.sf1 <- pred.locs.sf1 %>% 
        dplyr::arrange(ID, day) %>% 
        sf::st_geometry() %>% 
        sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(pred.locs.sf1$ID))) %>% 
        sf::st_cast("MULTILINESTRING") %>% 
        sf::st_sf(ID = as.factor(unique(pred.locs.sf1$ID))) %>% 
        separate(ID, into=c("pop", "ind"), remove=FALSE) %>% 
        sf::st_transform(4326)
      
      locs.start.sf <- locs.start.sf1 %>% 
        st_transform(4326)
      
      #25. Plot predictions----
      world <- ggplot() +
        borders("world", colour = "gray85", fill = "gray80") +
        theme_classic() +
        xlim(-150, -30) +
        ylim(-55, 72)
      
      plot1 <- world + 
        layer_spatial(pred.lines.sf1, size = 0.75, aes(color=factor(pop))) +
        layer_spatial(locs.start.sf, size = 0.75, col="red") +
        ggtitle("1. Predicted from randomly chosen points")
      
      ggsave(plot1, file=paste0(rd,"Figs/PredictedPaths1/PredictedPaths1_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".jpeg"), device="jpeg", width=8, height=8, units="in")
      
      #26. Simulate----
      times1.sim <- lapply(fit1,
                       function(x){
                         seq(lubridate::ceiling_date(min(as.POSIXct(x$data$day,tz ="GMT")),"hour"),
                             lubridate::floor_date(max(as.POSIXct(x$data$day,tz = "GMT")),"hour"),
                             "1 hour")})
      
      sim <- mapply(x = fit1,
                    y = times1.sim,
                    function(x,y){
                      crawl::crwSimulator(x,predTime = y)},
                    SIMPLIFY = FALSE)
      
      sim.list1 <- mapply(x = sim,
                          function(x){
                            crawl::crwPostIS(x, fullPost = FALSE)},
                          SIMPLIFY = FALSE)
      
      #27. Wrangle simulated points---
      sim.list2 <- purrr::flatten(sim.list1)
      sim.coords <- sim.list2[grepl("alpha.sim", names(sim.list2))] 
      names(sim.coords) <- names(sim.list1)
      sim.coords <- mapply(x=sim.coords,
                           function(x){
                             data.frame(x)},
                           SIMPLIFY = FALSE)
      sim.coords <- rbindlist(sim.coords, use.names=TRUE, idcol="ID")
      sim.days <- sim.list2[grepl("TimeNum", names(sim.list2))] 
      sim.days <- mapply(x=sim.days,
                         function(x){
                           data.frame(x)},
                         SIMPLIFY = FALSE)
      sim.days <- rbindlist(sim.days, use.names=TRUE, idcol=NULL) %>% 
        rename(day = x)
      
      sim.locs <- cbind(sim.coords, sim.days) %>% 
        mutate(day = as.Date(day, origin=as.Date("1970-01-01")))
      
      #28. Save simulated points----
#      write.csv(sim.locs, file=paste0(rd,"Data/SimulatedPaths/SimulatedPaths_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
      
      #29. Convert simulations to spatial features----
      sim.locs.sf <- sim.locs %>% sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
        sf::st_set_crs(3857)
      
      sim.lines.sf <- sim.locs.sf %>% 
        dplyr::arrange(ID, day) %>% 
        sf::st_geometry() %>% 
        sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(sim.locs.sf$ID))) %>% 
        sf::st_cast("MULTILINESTRING") %>% 
        sf::st_sf(ID = as.factor(unique(sim.locs.sf$ID))) %>% 
        separate(ID, into=c("pop", "ind"), remove=FALSE) %>% 
        sf::st_transform(4326)
      
      #30. Plot simulations----
      plot2 <- world + 
        layer_spatial(sim.lines.sf, size = 0.75, aes(color=factor(pop))) +
        ggtitle("2. Simulated from posterior of predictions")
      
      ggsave(plot2, file=paste0(rd,"Figs/SimulatedPaths/SimulatedPaths_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".jpeg"), device="jpeg", width=8, height=8, units="in")
      
      #31. Sample GPS points from simulation for crawl----
      period <- dates %>% 
        group_by(season) %>% 
        summarize(mind=min(day),
                  maxd=max(day))
      
      period.start.max <- period %>% 
        dplyr::filter(season=="breed") %>% 
        dplyr::select(maxd) %>% 
        as.integer
      
      period.start.min <- period %>% 
        dplyr::filter(season=="breed") %>% 
        dplyr::select(mind) %>% 
        as.integer
      
      period.end.min <- period %>% 
        dplyr::filter(season=="winter") %>% 
        dplyr::select(mind) %>% 
        as.integer
      
      period.end.max <- period %>% 
        dplyr::filter(season=="winter") %>% 
        dplyr::select(maxd) %>% 
        as.integer
      
      n.daysfix <- ifelse(floor(as.numeric(period.end.min - period.start.max, units="days")/n.gps)==0, 1,
                          floor(as.numeric(period.end.min - period.start.max, units="days")/n.gps))
      
      dates.sample <- data.frame(day = seq(period.start.min, period.end.max, n.daysfix)) %>% 
        mutate(date = as.Date(day, origin=as.Date("1970-01-01")),
               doy = yday(date))
      
      #This is where the spatial and temporal simulations separate out----
      #Spatial
      locs.gps.spat <- sim.locs %>% 
        mutate(doy = yday(day)) %>% 
        dplyr::filter(doy %in% dates.sample$doy) %>% 
        group_by(ID) %>% 
        sample_n(size=n.gps, replace=FALSE) %>% 
        mutate(season="migration") %>% 
        ungroup() %>% 
        separate(ID, into=c("pop", "ind"), remove=FALSE) %>% 
        dplyr::select(ID, pop, ind, mu.x, mu.y, day, season) %>% 
        dplyr::rename(X=mu.x, Y=mu.y) %>% 
        rbind(locs.start %>% 
                dplyr::filter(season=="breed") %>% 
                dplyr::select(ID, pop, ind, X, Y, day, season)) %>% 
        rbind(locs.start %>% 
                dplyr::filter(season=="winter") %>% 
                dplyr::select(ID, pop, ind, X, Y, day, season))
      
      #Temporal
      locs.gps.temp <- sim.locs %>% 
        mutate(doy = yday(day)) %>% 
        dplyr::filter(doy %in% dates.sample$doy,
                      hour(day)==0) %>% 
        group_by(ID) %>% 
        sample_n(size=n.gps, replace=FALSE) %>% 
        mutate(season="migration") %>% 
        ungroup() %>% 
        separate(ID, into=c("pop", "ind"), remove=FALSE) %>% 
        dplyr::select(ID, pop, ind, mu.x, mu.y, day, season) %>% 
        dplyr::rename(X=mu.x, Y=mu.y)
      
      #32. Error parameters for GPS locations-----
      #Spatial
      locs.gps.spat <- locs.gps.spat %>% 
        mutate(day = as.POSIXct(day)) %>% 
        dplyr::mutate(error_semi_major_axis=20,
                      error_semi_minor_axis=20,
                      error_ellipse_orientation=0)
      
      locs.gps.temp <- locs.gps.temp %>% 
        mutate(day = as.POSIXct(day)) %>% 
        dplyr::mutate(error_semi_major_axis=20,
                      error_semi_minor_axis=20,
                      error_ellipse_orientation=0)
      
      locs.gps.spat <- crawl::argosDiag2Cov(locs.gps.spat$error_semi_major_axis, locs.gps.spat$error_semi_minor_axis, locs.gps.spat$error_ellipse_orientation) %>% 
        cbind(locs.gps.spat) %>% 
        arrange(ID, day)
      
      #Temporal
      locs.gps.temp <- crawl::argosDiag2Cov(locs.gps.temp$error_semi_major_axis, locs.gps.temp$error_semi_minor_axis, locs.gps.temp$error_ellipse_orientation) %>% 
        cbind(locs.gps.temp) %>% 
        arrange(ID, day)
      
      write.csv(locs.gps.spat, file=paste0(rd,"Data/SampledPoints/SampledPoints_Spatial_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
      write.csv(locs.gps.temp, file=paste0(rd,"Data/SampledPoints/SampledPoints_Temporal_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
      
      #33. Wrangle for crawl----
      #Spatial
      locs.gps.spat.sf1 <- locs.gps.spat %>% 
        dplyr::select(ID, X, Y) %>% 
        sf::st_as_sf(coords = c("X","Y")) %>% 
        sf::st_set_crs(3857) 
      
      locs.spat.list2 <- split(locs.gps.spat, locs.gps.spat$ID)
      
      locs.spat.list.sf2 <- lapply(locs.spat.list2,
                                   function(x){SpatialPoints(
                                     cbind(x$X,
                                           x$Y), 
                                     proj4string = CRS("+init=epsg:3857"))})
      
      locs.spat.list.sf2 <- mapply(x=locs.spat.list.sf2,
                                   y=locs.spat.list2,
                                   function(x,y){
                                     SpatialPointsDataFrame(x,y)
                                   })
      
      #Temporal
      locs.gps.temp.sf1 <- locs.gps.temp %>% 
        dplyr::select(ID, X, Y) %>% 
        sf::st_as_sf(coords = c("X","Y")) %>% 
        sf::st_set_crs(3857) 
      
      locs.temp.list2 <- split(locs.gps.temp, locs.gps.temp$ID)
      
      locs.temp.list.sf2 <- lapply(locs.temp.list2,
                                   function(x){SpatialPoints(
                                     cbind(x$X,
                                           x$Y), 
                                     proj4string = CRS("+init=epsg:3857"))})
      
      locs.temp.list.sf2 <- mapply(x=locs.temp.list.sf2,
                                   y=locs.temp.list2,
                                   function(x,y){
                                     SpatialPointsDataFrame(x,y)
                                   })
      
      #34. Create starting values----
      #Spatial
      inits.spat2 <- lapply(locs.spat.list.sf2,
                            function(x){
                              list(a = c(sp::coordinates(x)[1,1],0,
                                         sp::coordinates(x)[1,2],0),
                                   P = diag(c(10 ^ 2, 10 ^ 2, 
                                              10 ^ 2, 10 ^ 2)))})
      
      #Temporal
      inits.temp2 <- lapply(locs.temp.list.sf2,
                            function(x){
                              list(a = c(sp::coordinates(x)[1,1],0,
                                         sp::coordinates(x)[1,2],0),
                                   P = diag(c(10 ^ 2, 10 ^ 2, 
                                              10 ^ 2, 10 ^ 2)))})
      
      fixpar <- c(1,1,NA,NA)
      
      #35. Fit with crawl----
      #Spatial
      fit.spat2 <- mapply(x = locs.spat.list.sf2,
                          y = inits.spat2,
                          function(x,y){
                            try(crawl::crwMLE(mov.model =  ~ 1,
                                              err.model = list(
                                                x =  ~ ln.sd.x - 1, 
                                                y =  ~ ln.sd.y - 1, 
                                                rho =  ~ error.corr),
                                              data = x,
                                              Time.name = "day",
                                              time.scale = "days",
                                              initial.state = y,
                                              fixPar = fixpar,
                                              initialSANN = list(maxit = 2500),
                                              need.hess = TRUE,
                                              theta=c(14,-1),
                                              control = list(REPORT = 1000, trace = 6),
                                              attempts = 8))},
                          SIMPLIFY = FALSE)
      
      #Temporal
      fit.temp2 <- mapply(x = locs.temp.list.sf2,
                          y = inits.temp2,
                          function(x,y){
                            try(crawl::crwMLE(mov.model =  ~ 1,
                                              err.model = list(
                                                x =  ~ ln.sd.x - 1, 
                                                y =  ~ ln.sd.y - 1, 
                                                rho =  ~ error.corr),
                                              data = x,
                                              Time.name = "day",
                                              time.scale = "days",
                                              initial.state = y,
                                              fixPar = fixpar,
                                              initialSANN = list(maxit = 2500),
                                              need.hess = TRUE,
                                              theta=c(14,-1),
                                              control = list(REPORT = 1000, trace = 6),
                                              attempts = 8))},
                          SIMPLIFY = FALSE)
      
      
      
      #36. Save crawl parameters----
      #Spatial
      params.spat2 <- lapply(fit.spat2,
                             function(x){
                               crawl::tidy_crwFit(x)
                             })
      params.spat.summary2 <- bind_rows(params.spat2, .id = "column_label")
      
      write.csv(params.spat.summary2, file=paste0(rd,"CrawlParameters/PredictedPaths2/CrawlParameters2_Spatial_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
      
      params.temp2 <- lapply(fit.temp2,
                             function(x){
                               crawl::tidy_crwFit(x)
                             })
      params.temp.summary2 <- bind_rows(params.temp2, .id = "column_label")
      
      write.csv(params.temp.summary2, file=paste0(rd,"CrawlParameters/PredictedPaths2/CrawlParameters2_Temporal_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
      
      #37. Predict----
      #Spatial
      times.spat2 <- lapply(fit.spat2,
                            function(x){
                              seq(lubridate::ceiling_date(min(as.POSIXct(x$data$day,tz ="GMT")),"hour"),
                                  lubridate::floor_date(max(as.POSIXct(x$data$day,tz = "GMT")),"hour"),
                                  "1 hour")})
      
      pred.spat.list2 <- mapply(x = fit.spat2,
                                y = times.spat2,
                                function(x,y){
                                  crawl::crwPredict(x,predTime = y,return.type = "flat")},
                                SIMPLIFY = FALSE)
      
      pred.spat2 <- data.frame(do.call(rbind, pred.spat.list2)) %>% 
        dplyr::mutate(up.x = mu.x+se.mu.x,
                      lw.x = mu.x-se.mu.x,
                      up.y = mu.y+se.mu.y,
                      lw.y = mu.y-se.mu.y)
      
      #Temporal
      times.temp2 <- lapply(fit.temp2,
                            function(x){
                              seq(lubridate::ceiling_date(min(as.POSIXct(x$data$day,tz ="GMT")),"hour"),
                                  lubridate::floor_date(max(as.POSIXct(x$data$day,tz = "GMT")),"hour"),
                                  "1 hour")})
      
      pred.temp.list2 <- mapply(x = fit.temp2,
                                y = times.temp2,
                                function(x,y){
                                  crawl::crwPredict(x,predTime = y,return.type = "flat")},
                                SIMPLIFY = FALSE)
      
      pred.temp2 <- data.frame(do.call(rbind, pred.temp.list2)) %>% 
        dplyr::mutate(up.x = mu.x+se.mu.x,
                      lw.x = mu.x-se.mu.x,
                      up.y = mu.y+se.mu.y,
                      lw.y = mu.y-se.mu.y)
      
      #38. Check if any individuals didn't work----
      #Spatial
      fit.spat.na2 <- pred.spat2 %>% 
        dplyr::filter(is.na(se.mu.x))
      
      param.spat.na2 <- params.spat.summary2 %>% 
        dplyr::filter((term=="ln beta (Intercept)" &
                         std.error=="NaN") |
                        (term=="ln sigma (Intercept)" &
                           std.error=="NaN") |
                        (term=="ln beta (Intercept)" &
                           is.na(std.error)) |
                        (term=="ln sigma (Intercept)" &
                           is.na(std.error)))
      
      #Temporal
      fit.temp.na2 <- pred.temp2 %>% 
        dplyr::filter(is.na(se.mu.x))
      
      param.temp.na2 <- params.temp.summary2 %>% 
        dplyr::filter((term=="ln beta (Intercept)" &
                         std.error=="NaN") |
                        (term=="ln sigma (Intercept)" &
                           std.error=="NaN") |
                        (term=="ln beta (Intercept)" &
                           is.na(std.error)) |
                        (term=="ln sigma (Intercept)" &
                           is.na(std.error)))
      
      if(nrow(fit.spat.na2) > 0 | nrow(param.spat.na2) > 0 | nrow(fit.temp.na2) > 0 | nrow(param.temp.na2) > 0) {
        
        next
      }
      
      else{
      
      #39. Convert predictions to spatial features----
      #Spatial
      pred.locs.spat.sf2 <- sf::st_as_sf(pred.spat2, coords = c("mu.x","mu.y")) %>% 
        sf::st_set_crs(3857)
      
      pred.lines.spat.sf2 <- pred.locs.spat.sf2 %>% 
        dplyr::arrange(ID, day) %>% 
        sf::st_geometry() %>% 
        sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(pred.locs.spat.sf2$ID))) %>% 
        sf::st_cast("MULTILINESTRING") %>% 
        sf::st_sf(ID = as.factor(unique(pred.locs.spat.sf2$ID))) %>% 
        separate(ID, into=c("pop", "ind"), remove=FALSE) %>% 
        sf::st_transform(4326)
      
      locs.gps.spat.sf <- locs.gps.spat %>% 
        dplyr::select(ID, X, Y) %>% 
        sf::st_as_sf(coords = c("X","Y")) %>% 
        separate(ID, into=c("pop", "ind"), remove=FALSE) %>% 
        sf::st_set_crs(3857) %>% 
        sf::st_transform(4326)
      
      #Temporal
      pred.locs.temp.sf2 <- sf::st_as_sf(pred.temp2, coords = c("mu.x","mu.y")) %>% 
        sf::st_set_crs(3857)
      
      pred.lines.temp.sf2 <- pred.locs.temp.sf2 %>% 
        dplyr::arrange(ID, day) %>% 
        sf::st_geometry() %>% 
        sf::st_cast("MULTIPOINT",ids = as.integer(as.factor(pred.locs.temp.sf2$ID))) %>% 
        sf::st_cast("MULTILINESTRING") %>% 
        sf::st_sf(ID = as.factor(unique(pred.locs.temp.sf2$ID))) %>% 
        separate(ID, into=c("pop", "ind"), remove=FALSE) %>% 
        sf::st_transform(4326)
      
      locs.gps.temp.sf <- locs.gps.temp %>% 
        dplyr::select(ID, X, Y) %>% 
        sf::st_as_sf(coords = c("X","Y")) %>% 
        separate(ID, into=c("pop", "ind"), remove=FALSE) %>% 
        sf::st_set_crs(3857) %>% 
        sf::st_transform(4326)
      
      #40. Plot predictions----
      world <- ggplot() +
        borders("world", colour = "gray85", fill = "gray80") +
        theme_classic() +
        xlim(-150, -30) +
        ylim(-55, 72)
      
      #Spatial
      plot.spat3 <- world + 
        layer_spatial(pred.lines.spat.sf2, size = 0.75, aes(color=factor(pop))) +
        layer_spatial(locs.gps.spat.sf, size = 0.75, col="red") +
        ggtitle("3. Predicted from points sampled from simulated tracks - Spatial")
      
      ggsave(plot.spat3, file=paste0(rd,"Figs/PredictedPaths2/PredictedPaths2_Spatial_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".jpeg"), device="jpeg", width=8, height=8, units="in")
      
      #Temporal
      plot.temp3 <- world + 
        layer_spatial(pred.lines.temp.sf2, size = 0.75, aes(color=factor(pop))) +
        layer_spatial(locs.gps.temp.sf, size = 0.75, col="red") +
        ggtitle("3. Predicted from points sampled from simulated tracks - temporal")
      
      ggsave(plot.temp3, file=paste0(rd,"Figs/PredictedPaths2/PredictedPaths2_Temporal_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".jpeg"), device="jpeg", width=8, height=8, units="in")
      
      #41. Wrangle data for connectivity----
      locs.start.breed <- locs.start %>% 
        dplyr::filter(season=="breed") %>% 
        dplyr::select(ID, X, Y) %>% 
        dplyr::rename(bandlong = X, bandlat = Y)
      
      #Spatial
      pred.mc.spat <- pred.spat2 %>% 
        st_as_sf(coords=c("mu.x", "mu.y"), crs=3857) %>% 
        st_transform(crs=4326) %>% 
        st_coordinates() %>% 
        data.frame() %>% 
        dplyr::rename(long=X, lat=Y) %>% 
        cbind(pred.spat2) %>% 
        mutate(season=ifelse(is.na(season), "migration", season)) %>% 
        dplyr::select(ID, pop, ind, day, mu.x, mu.y, se.mu.x, se.mu.y, up.x, up.y, lw.x, lw.y, locType, X, Y, long, lat, season) %>% 
        mutate(latr = round(lat),
               latdiff=abs(latr-lat)) %>% 
        group_by(ID, latr) %>% 
        mutate(n=n(),
               mindiff=min(latdiff)) %>% 
        dplyr::filter(mindiff==latdiff) %>% 
        mutate(n=n()) %>% 
        ungroup() %>% 
        left_join(locs.start.breed)
      
      write.csv(pred.mc.spat, file=paste0(rd,"Data/SpatialConnectivityPoints/SpatialConnectivityPoints_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
      
      #Temporal
      pred.mc.temp <- pred.temp2 %>% 
        st_as_sf(coords=c("mu.x", "mu.y"), crs=3857) %>% 
        st_transform(crs=4326) %>% 
        st_coordinates() %>% 
        data.frame() %>% 
        dplyr::rename(long=X, lat=Y) %>% 
        cbind(pred.temp2) %>% 
        mutate(season=ifelse(is.na(season), "migration", season)) %>% 
        dplyr::select(ID, pop, ind, day, mu.x, mu.y, se.mu.x, se.mu.y, up.x, up.y, lw.x, lw.y, locType, X, Y, long, lat, season) %>%
        mutate(hour=hour(day),
               doy=yday(day)) %>% 
        dplyr::filter(hour==0) %>% 
        group_by(ID) %>% 
        mutate(firstdoy = min(doy),
               lastdoy = max(doy)) %>% 
        ungroup() %>% 
        left_join(locs.start.breed)
      
      write.csv(pred.mc.temp, file=paste0(rd,"Data/TemporalConnectivityPoints/TemporalConnectivityPoints_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
      
      #42. Spatial connectivity loop----
      pops <- c(unique(pred.mc.spat$pop), 99)
      
      mantel.spat.flat <- list()
      for(h in 1:length(pops)){
        
        pop.h <- pops[h]
        
        dat.spat.h <- pred.mc.spat %>% 
          dplyr::filter(pop != pop.h)
        
        #46. Filter out latitudes with only one population----
        lats <- dat.spat.h %>% 
          group_by(latr, pop) %>% 
          summarize(count = n()) %>% 
          group_by(latr) %>% 
          summarize(count.pop=n(),
                    count.ind=sum(count),
                    min.ind=min(count)) %>% 
          mutate(percent.pop=count.pop/max(count.pop),
                 percent.ind=count.ind/max(count.ind)) %>% 
          filter(count.pop > 1,
                 count.ind > 2)
        
        #43. Set up inner loop----
        mantel.spat <- list()
        for(i in 1:nrow(lats)){
          
          #Set latitude & subset data
          lat.i <- lats$latr[i]
          dat.spat.i <- dat.spat.h %>% 
            dplyr::filter(latr==lat.i)
          
          #Convert origin points to sp
          originPoints <- dat.spat.i %>% 
            dplyr::select(bandlong, bandlat) %>% 
            rename(Long=bandlong, Lat=bandlat) %>% 
            st_as_sf(coords=c("Long", "Lat"), crs=3857) %>% 
            st_geometry() %>% 
            as_Spatial()
          
          #Convert target points to sp
          targetPoints <- dat.spat.i %>% 
            dplyr::select(mu.x, mu.y) %>% 
            rename(Long=mu.x, Lat=mu.y) %>% 
            st_as_sf(coords=c("Long", "Lat"), crs=3857) %>% 
            st_geometry() %>% 
            as_Spatial()
          
          #Define tag type
          isGLDF <- dat.spat.i %>% 
            mutate(GL=case_when(locType=="p" ~ TRUE,
                                locType=="o" ~ FALSE))
          isGL <- isGLDF$GL
          
          #Define error
          targetError = dat.spat.i$se.mu.x
          
          #Calculate connectivity
          Man <- try(estMantel(originPoints = originPoints,
                               targetPoints = targetPoints,
                               targetError = targetError,
                               nBoot = 100))
          
          if (class(Man)=="estMantel"){
            
            #Save results
            man <- Man[["sampleCorr"]]
            mantel.i <- data.frame(mantel=man) %>% 
              mutate(lat=lat.i,
                     l1o=pop.h)
            mantel.spat[[i]] <- mantel.i
            
          }
          
          #Report status
          print(paste0("Completed ", i, " of ", nrow(lats), " latitudes: ", lat.i))
          
        }
        
        #44. Flatten spatial connectivity results----
        mantel.spat.flat[[h]] <- rbindlist(mantel.spat)
        
        
        print(paste0("COMPLETED ", h, " of ", length(pops), " LEAVE ONE OUT: POPULATION ", pop.h))
        
      }
      
      
      #45. Temporal connectivity loop----
      days <- pred.mc.temp %>% 
        group_by(doy, pop) %>% 
        summarize(count = n()) %>% 
        group_by(doy) %>% 
        summarize(count.pop=n(),
                  count.ind=sum(count)) %>% 
        mutate(percent.pop=count.pop/max(count.pop),
               percent.ind=count.ind/max(count.ind)) %>% 
        filter(count.pop > 1,
               count.ind > 2) %>% 
        mutate(use=ifelse(percent.pop==lag(percent.pop), 0, 1)) %>% 
        filter(use==1)
      
      days <- c(days$doy, 99)
      
      mantel.temp.flat <- list()
      for(h in 1:length(days)){
        
        day.h <- days[h]
        
        dat.temp.h <- pred.mc.temp %>% 
          dplyr::filter(firstdoy != day.h,
                        lastdoy != day.h)
        
        #46. Filter out days with only one population----
        doys <- dat.temp.h %>% 
          group_by(doy, pop) %>% 
          summarize(count = n()) %>% 
          group_by(doy) %>% 
          summarize(count.pop=n(),
                    count.ind=sum(count),
                    min.ind=min(count)) %>% 
          mutate(percent.pop=count.pop/max(count.pop),
                 percent.ind=count.ind/max(count.ind)) %>% 
          filter(count.pop > 1,
                 count.ind > 2)  
        
        #47. Set up inner loop----
        mantel.temp <- list()
        for(i in 1:nrow(doys)){
          
          #Set doy & subset data
          doy.i <- doys$doy[i]
          dat.temp.i <- dat.temp.h %>% 
            dplyr::filter(doy==doy.i)
          
          #48. Create temporal connectivity inputs----
          
          #Convert origin points to sp
          originPoints <- dat.temp.i %>% 
            dplyr::select(bandlong, bandlat) %>% 
            rename(Long=bandlong, Lat=bandlat) %>% 
            st_as_sf(coords=c("Long", "Lat"), crs=3857) %>% 
            st_geometry() %>% 
            as_Spatial()
          
          #Convert target points to sp
          targetPoints <- dat.temp.i %>% 
            dplyr::select(mu.x, mu.y) %>% 
            rename(Long=mu.x, Lat=mu.y) %>% 
            st_as_sf(coords=c("Long", "Lat"), crs=3857) %>% 
            st_geometry() %>% 
            as_Spatial()
          
          #Define tag type
          isGLDF <- dat.temp.i %>% 
            mutate(GL=case_when(locType=="p" ~ TRUE,
                                locType=="o" ~ FALSE))
          isGL <- isGLDF$GL
          
          #Define error
          targetError = dat.temp.i$se.mu.x
          
          #Calculate connectivity
          Man <- try(estMantel(originPoints = originPoints,
                               targetPoints = targetPoints,
                               targetError = targetError,
                               nBoot = 100))
          
          if (class(Man)=="estMantel"){
            
            #Save results
            man <- Man[["sampleCorr"]]
            mantel.i <- data.frame(mantel=man) %>% 
              mutate(doy=doy.i,
                     l1o=day.h)
            mantel.temp[[i]] <- mantel.i
            
          } 
          
          #Report status
          print(paste0("Completed ", i, " of ", nrow(doys), " migration days: ", doy.i))
          
        }
        
        #49. Flatten temporal connectivity results----
        mantel.temp.flat[[h]] <- rbindlist(mantel.temp)
        
        
        print(paste0("COMPLETED ", h, " of ", length(days), " LEAVE ONE OUT: DATE ", day.h))
        
      }
      
      #50. Flatten and save results for this iteration----
      mantel.spat.flat2 <- rbindlist(mantel.spat.flat)
      write.csv(mantel.spat.flat2, file=paste0(rd,"Results/MCSpatial/MC_Spatial_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
      
      mantel.temp.flat2 <- rbindlist(mantel.temp.flat)
      write.csv(mantel.temp.flat2, file=paste0(rd,"Results/MCTemporal/MC_Temporal_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".csv"), row.names=FALSE)
      
      #51. Summarize results----
      mantel.spat.sum <- data.frame(mantel.spat.flat2) %>% 
        dplyr::filter(!is.na(mantel)) %>% 
        group_by(lat, l1o) %>% 
        summarize(Mmean=mean(mantel),
                  Mlowq=quantile(mantel, probs=0.085),
                  Mhighq=quantile(mantel, probs=0.915)) %>% 
        ungroup()
      
      mantel.temp.sum <- data.frame(mantel.temp.flat2) %>% 
        dplyr::filter(!is.na(mantel)) %>% 
        group_by(doy, l1o) %>% 
        summarize(Mmean=mean(mantel),
                  Mlowq=quantile(mantel, probs=0.085),
                  Mhighq=quantile(mantel, probs=0.915)) %>% 
        ungroup()
      
      #52. Visualize profiles----
      
      mantelplot1 <- ggplot() +
        geom_hex(aes(x=lat, y=mantel), data=mantel.spat.flat2) +
        geom_smooth(aes(x=lat, y=mantel), data=mantel.spat.flat2) +
        scale_fill_viridis_c() +
        facet_wrap(~l1o) +
        theme(legend.position = "bottom") +
        ggtitle("Spatial")
      
      mantelplot2 <- ggplot(mantel.spat.sum) +
        geom_ribbon(aes(x=lat, ymin=Mlowq, ymax=Mhighq), alpha=0.4) +
        geom_line(aes(x=lat, y=Mmean)) +
        geom_smooth(aes(x=lat, y=mantel), data=mantel.spat.flat2) +
        facet_wrap(~l1o) +
        ggtitle("Spatial")
      
      mantelplot3 <- ggplot() +
        geom_hex(aes(x=doy, y=mantel), data=mantel.temp.flat2) +
        geom_smooth(aes(x=doy, y=mantel), data=mantel.temp.flat2) +
        scale_fill_viridis_c() +
        facet_wrap(~l1o) +
        theme(legend.position = "bottom") +
        ggtitle("Temporal")
      
      mantelplot4 <- ggplot(mantel.temp.sum) +
        geom_ribbon(aes(x=doy, ymin=Mlowq, ymax=Mhighq), alpha=0.4) +
        geom_line(aes(x=doy, y=Mmean)) +
        geom_smooth(aes(x=doy, y=mantel), data=mantel.temp.flat2) +
        facet_wrap(~l1o) +
        ggtitle("Temporal")
      
      #53. Save profiles----
      grid <- grid.arrange(mantelplot1, mantelplot2, mantelplot3, mantelplot4, ncol=2, nrow=2)
      
      ggsave(grid, file=paste0(rd,"Figs/Profile/Profile_",n.pop,"-",n.ind,"-",n.gps,"_",n.sim,".jpeg"), device="jpeg", width=20, height=20, units="in")
      
      end.i <- Sys.time()
      
      time.i <- end.i - start.i
      
      times <- rbind(times, data.frame(time=time.i,
                                       start=start.i,
                                       end=end.i,
                                       pop=n.pop,
                                       ind=n.ind,
                                       gps=n.gps,
                                       sim=n.sim))
      
      write.csv(times, file=paste0(rd,"RunTimes.csv"), row.names=FALSE)
      
      print(paste0("COMPLETED ", n.sim, " of ", set.sim, " SIMULATIONS"))
      
    }

    }
    
  }
    
  
  print(paste0("COMPLETED ", f, " of ", nrow(settings), " SETTING COMBINATIONS"))
    
}


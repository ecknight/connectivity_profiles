#title: Figures from Knight et al. 2021. Comprehensive estimation of spatial and temporal migratory connectivity across the annual cycle to direct conservation efforts. Ecography ECOG-05111
#author: Elly C. Knight
#date: Jan. 17, 2021

options(scipen = 99999)

library(tidyverse)
library(sf)
library(maps)
library(Cairo)
library(gridExtra)
library(cowplot)
library(raster)
library(lubridate)
library(ggspatial)
library(gganimate)
library(ggmap)
library(ggforce)
library(gridExtra)
library(grid)

my.theme <- theme_classic() +
  theme(text=element_text(size=12, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(size=12, hjust = 0.5))

map.theme <- theme_nothing() +
  theme(text=element_text(size=12, family="Arial"),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.text = element_blank())

whemi <- map_data("world", region=c("Canada", 
                                    "USA", 
                                    "Mexico",
                                    "Guatemala", 
                                    "Belize", 
                                    "El Salvador",
                                    "Honduras", 
                                    "Nicaragua", 
                                    "Costa Rica",
                                    "Panama", 
                                    "Jamaica", 
                                    "Cuba", 
                                    "The Bahamas",
                                    "Haiti", 
                                    "Dominican Republic", 
                                    "Antigua and Barbuda",
                                    "Dominica", 
                                    "Saint Lucia", 
                                    "Saint Vincent and the Grenadines", 
                                    "Barbados",
                                    "Grenada",
                                    "Trinidad and Tobago",
                                    "Colombia",
                                    "Venezuela",
                                    "Guyana",
                                    "Suriname",
                                    "Ecuador",
                                    "Peru",
                                    "Brazil",
                                    "Bolivia",
                                    "Paraguay",
                                    "Chile",
                                    "Argentina",
                                    "Uruguay")) %>% 
  dplyr::filter(!group%in%c(258:264))

nam <- map_data("world", region=c("Canada", 
                                           "USA", 
                                           "Mexico",
                                           "Guatemala", 
                                           "Belize", 
                                           "El Salvador",
                                           "Honduras", 
                                           "Nicaragua", 
                                           "Costa Rica",
                                           "Panama", 
                                           "Jamaica", 
                                           "Cuba", 
                                           "The Bahamas",
                                           "Haiti", 
                                           "Dominican Republic", 
                                           "Antigua and Barbuda",
                                           "Dominica", 
                                           "Saint Lucia", 
                                           "Saint Vincent and the Grenadines", 
                                           "Barbados",
                                           "Grenada",
                                           "Trinidad and Tobago")) %>% 
  dplyr::filter(!group%in%c(258:264))

sam <- map_data("world", region=c(
                                    "Trinidad and Tobago",
                                    "Colombia",
                                    "Venezuela",
                                    "Guyana",
                                    "Suriname",
                                    "Ecuador",
                                    "Peru",
                                    "Brazil",
                                    "Bolivia",
                                    "Paraguay",
                                    "Chile",
                                    "Argentina",
                                    "Uruguay")) %>% 
  dplyr::filter(!group%in%c(258:264))

#1. Figure 1 - Theoretical figure----
country <- map_data("world", region=c("Canada", 
                                  "USA"))

state <- map_data("state")

lake <- map_data("lakes")



pop.breed <- data.frame(ID = c(1:3),
                    lat = c(45, 45, 45),
                    long = c(-85, -90, -95))

pop.stop.low <- data.frame(ID = 1,
                   lat = 30,
                   long = -90)

pop.stop.high <- data.frame(ID = c(1:3),
                        lat = rep(30, 3),
                        long = c(-85, -90, -95))

ind <- 25

ind.high.high <- data.frame(ID = c(rep(1, ind), rep(2, ind), rep(3,ind)),
                            lat = c(runif(ind, min=29.5, max=31),
                                        runif(ind, min=34, max=36),
                                        runif(ind, min=39, max=41)),
                            long = c(runif(ind, min=-96, max=-94),
                                     runif(ind, min=-91, max=-89),
                                     runif(ind, min=-86, max=-84)))

ind.high.low <- data.frame(ID = c(rep(1, ind), rep(2, ind), rep(3,ind)),
                           lat = c(runif(ind, min=34, max=36),
                                   runif(ind, min=34, max=36),
                                   runif(ind, min=34, max=36)),
                           long = c(runif(ind, min=-96, max=-94),
                                    runif(ind, min=-91, max=-89),
                                    runif(ind, min=-86, max=-84)))

ind.low.high <- data.frame(ID = c(rep(1, ind), rep(2, ind), rep(3,ind)),
                            lat = c(runif(ind, min=30, max=31),
                                    runif(ind, min=34, max=36),
                                    runif(ind, min=39, max=41)),
                            long = c(runif(ind, min=-94, max=-86),
                                     runif(ind, min=-93, max=-87),
                                     runif(ind, min=-89, max=-85)))

ind.low.low <- data.frame(ID = c(rep(1, ind), rep(2, ind), rep(3,ind)),
                          lat = c(runif(ind, min=34, max=36),
                                  runif(ind, min=34, max=36),
                                  runif(ind, min=34, max=36)),
                           long = c(runif(ind, min=-94, max=-86),
                                    runif(ind, min=-94, max=-86),
                                    runif(ind, min=-94, max=-86)))
scale <- ggplot() +
  scale_x_continuous(limits=c(0, 1),
                     breaks=c(0, 1)) +
  xlab("Migratory\nconnectivity") +
  ylim(c(0, 0)) +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=10))

highspat.hightemp <- ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray70", colour = "gray70", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  geom_segment(aes(x=-85, xend=-85, y=45, yend=41.2), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-85, xend=-85, y=38.8, yend=32), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_segment(aes(x=-90, xend=-90, y=45, yend=36.2), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-90, xend=-90, y=33.8, yend=32), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_segment(aes(x=-95, xend=-95, y=45, yend=32), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_point(data=pop.breed, aes(x=long, y=lat, fill=ID), colour="black", pch=21, size=14) +
  geom_point(data=pop.stop.high, aes(x=long, y=lat, fill=ID), colour="black", pch=21, size=14, alpha=0.3) +
  geom_point(data=pop.stop.high, aes(x=long, y=lat), fill=NA, colour="black", pch=21, size=14) +
  geom_point(data=ind.high.high, aes(x=long, y=lat, colour=ID), size=1) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(direction=-1) +
  theme(legend.position = "none")+
  coord_fixed(xlim = c(-100, -65),  ylim = c(25, 47), ratio = 1.3) +
  map.theme +
  annotation_custom(grob = ggplotGrob(scale), xmin=-80, xmax=-65, ymin=28, ymax=28) +
  geom_point(aes(x=-66.8, y=31.2), size=4, colour="black") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
#highspat.hightemp
  
highspat.lowtemp <-  ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray70", colour = "gray70", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  geom_segment(aes(x=-85, xend=-85, y=45, yend=36), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-85, xend=-85, y=34, yend=32), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_segment(aes(x=-90, xend=-90, y=45, yend=36), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-90, xend=-90, y=34, yend=32), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_segment(aes(x=-95, xend=-95, y=45, yend=36), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-95, xend=-95, y=34, yend=32), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_point(data=pop.breed, aes(x=long, y=lat, fill=ID), colour="black", pch=21, size=14) +
  geom_point(data=pop.stop.high, aes(x=long, y=lat, fill=ID), colour="black", pch=21, size=14, alpha=0.3) +
  geom_point(data=pop.stop.high, aes(x=long, y=lat), fill=NA, colour="black", pch=21, size=14) +
  geom_point(data=ind.high.low, aes(x=long, y=lat, colour=ID), size=1) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(direction=-1) +
  theme(legend.position = "none")+
  coord_fixed(xlim = c(-100, -65),  ylim = c(25, 47), ratio = 1.3) +
  map.theme +
  annotation_custom(grob = ggplotGrob(scale), xmin=-80, xmax=-65, ymin=28, ymax=28) +
  geom_point(aes(x=-66.8, y=31.2), size=4, colour="black") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
#highspat.lowtemp

lowspat.hightemp <- ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray70", colour = "gray70", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  geom_segment(aes(x=-85, xend=-86.1, y=45, yend=41.2), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-86.8, xend=-89, y=38.8, yend=31.5), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_segment(aes(x=-90, xend=-90, y=45, yend=36.2), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-90, xend=-90, y=33.8, yend=31.5), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_segment(aes(x=-95, xend=-91, y=45, yend=31.5), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_point(data=pop.breed, aes(x=long, y=lat, fill=ID), colour="black", pch=21, size=14) +
  geom_ellipse(aes(x0 = -90, y0 = 30, a = 5, b = 1.5, angle = 0), fill="black", alpha=0.3) +
  geom_point(data=ind.low.high, aes(x=long, y=lat, colour=ID), size=1) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(direction=-1) +
  theme(legend.position = "none")+
  coord_fixed(xlim = c(-100, -65),  ylim = c(25, 47), ratio = 1.3) +
  map.theme +
  annotation_custom(grob = ggplotGrob(scale), xmin=-80, xmax=-65, ymin=28, ymax=28) +
  geom_point(aes(x=-66.8, y=31.2), size=4, colour="black") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
#lowspat.hightemp

lowspat.lowtemp <- ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray70", colour = "gray70", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  geom_segment(aes(x=-85, xend=-87.6, y=45, yend=36), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-88.2, xend=-89, y=34, yend=31.5), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_segment(aes(x=-90, xend=-90, y=45, yend=36), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-90, xend=-90, y=34, yend=31.5), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_segment(aes(x=-95, xend=-92.4, y=45, yend=36), arrow=arrow(length=unit(0.15, "inches"))) +
  geom_segment(aes(x=-91.8, xend=-91, y=34, yend=31.5), arrow=arrow(length=unit(0.15, "inches")), linetype="dashed") +
  geom_point(data=pop.breed, aes(x=long, y=lat, fill=ID), colour="black", pch=21, size=14) +
  geom_ellipse(aes(x0 = -90, y0 = 30, a = 5, b = 1.5, angle = 0), fill="black", alpha=0.3) +
  geom_point(data=ind.low.low, aes(x=long, y=lat, colour=ID), size=1) +
  scale_fill_viridis_c() +
  scale_colour_viridis_c(direction=-1) +
  coord_fixed(xlim = c(-100, -65),  ylim = c(25, 47), ratio = 1.3) +
  map.theme +
  annotation_custom(grob = ggplotGrob(scale), xmin=-80, xmax=-65, ymin=28, ymax=28) +
  geom_point(aes(x=-77.7, y=31.2), size=4, colour="black") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
lowspat.lowtemp

text.highspat.1 <- grobTree(textGrob("High Spatial\nConnectivity",
                                      x=0.9,  y=0.85, just="right",
                                      gp=gpar(fontsize=14, fontface="bold")))
text.highspat.2 <- grobTree(textGrob("Individuals from\nmultiple breeding\npopulations visit\ndifferent stopover\nlocations",
                                      x=0.9,  y=0.57, just="right",
                                      gp=gpar(fontsize=10, fontface="italic")))
text.highspat.plot <- ggplot() +
  annotation_custom(text.highspat.1) +
  annotation_custom(text.highspat.2) +
  theme_nothing()

text.lowspat.1 <- grobTree(textGrob("Low Spatial\nConnectivity",
                                      x=0.9,  y=0.85, just="right",
                                      gp=gpar(fontsize=14, fontface="bold")))
text.lowspat.2 <- grobTree(textGrob("Individuals from\nmultiple breeding\npopulations share\nthe same stopover\nlocation",
                                      x=0.9,  y=0.57, just="right",
                                      gp=gpar(fontsize=10, fontface="italic")))
text.lowspat.plot <- ggplot() +
  annotation_custom(text.lowspat.1) +
  annotation_custom(text.lowspat.2) +
  theme_nothing()

text.hightemp.1 <- grobTree(textGrob("High Temporal Connectivity",
                                     x=0,  y=0.9, just="left",
                                     gp=gpar(fontsize=14, fontface="bold")))
text.hightemp.2 <- grobTree(textGrob("Individuals from multiple breeding\npopulations arrive at stopover\nlocations at different times",
                                     x=0,  y=0.5, just="left",
                                     gp=gpar(fontsize=10, fontface="italic")))
text.hightemp.plot <- ggplot() +
  annotation_custom(text.hightemp.1) +
  annotation_custom(text.hightemp.2) +
  theme_nothing()

text.lowtemp.1 <- grobTree(textGrob("Low Temporal Connectivity",
                                     x=0,  y=0.9, just="left",
                                     gp=gpar(fontsize=14, fontface="bold")))
text.lowtemp.2 <- grobTree(textGrob("Individuals from multiple breeding\npopulations arrive at stopover\nlocations at similar times",
                                     x=0,  y=0.5, just="left",
                                     gp=gpar(fontsize=10, fontface="italic")))
text.lowtemp.plot <- ggplot() +
  annotation_custom(text.lowtemp.1) +
  annotation_custom(text.lowtemp.2) +
  theme_nothing()

Fig1 <- grid.arrange(text.highspat.plot, highspat.hightemp, highspat.lowtemp,
             text.lowspat.plot, lowspat.hightemp, lowspat.lowtemp, 
             text.hightemp.plot, text.lowtemp.plot,
             widths=c(1.5,3,3),
             heights=c(2.5,2.5,1),
             layout_matrix=rbind(c(1,2,3),
                                 c(4,5,6),
                                 c(NA,7,8)))

ggsave(plot=Fig1, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs/Fig1HighLow.jpeg", width=7.5, height=7, units="in", device="jpeg")

#2. Figure 2 - Study area----
setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Data/Abundance&Trend/Cluster")
nat <- read_sf("NaturalPops.shp")

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data")
pop <- read.csv("tbl_population_abundance.csv") %>% 
  dplyr::filter(Region != "Florida") %>% 
  dplyr::mutate(Region = case_when(Region=="BC coast" ~ "Coastal BC",
                                   Region=="BC Okanagan" ~ "Southcentral BC",
                                   !is.na(Region) ~ as.character(Region))) %>% 
  arrange(-Lat) %>% 
  mutate(order=row_number()) 

tags1 <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  dplyr::select(PinpointID, Population) %>% 
  unique() %>% 
  left_join(pop)

tags2 <- table(tags1$Region) %>% 
  as.data.frame() %>% 
  rename(Region = Var1,
         NTagsAnalysis = Freq) %>% 
  mutate(Region = as.character(Region))

pop2 <- pop %>% 
  left_join(tags2) %>% 
  mutate(NTagsAnalysis = ifelse(is.na(NTagsAnalysis), 0, NTagsAnalysis))
  
pop.sf <- pop2 %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=crs(nat)) 

map <- map_data("world", region=c("Canada", 
                                  "USA", 
                                  "Mexico",
                                  "Guatemala", 
                                  "Belize", 
                                  "El Salvador",
                                  "Honduras", 
                                  "Nicaragua", 
                                  "Costa Rica",
                                  "Panama", 
                                  "Jamaica", 
                                  "Cuba", 
                                  "The Bahamas",
                                  "Haiti", 
                                  "Dominican Republic", 
                                  "Antigua and Barbuda",
                                  "Dominica", 
                                  "Saint Lucia", 
                                  "Saint Vincent and the Grenadines", 
                                  "Barbados",
                                  "Grenada",
                                  "Trinidad and Tobago")) %>% 
  dplyr::filter(!group%in%c(258:264)) %>% 
  dplyr::filter(!group%in%c(158:164))

base <- ggplot() +
  geom_polygon(data=map, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3)

range <- read_sf("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/GIS/CONI.shp")

studyarea <- ggplot() +
  geom_polygon(data=map, aes(x=long, y=lat, group=group), fill="gray70", colour = "gray85", size=0.3) +
   geom_sf(data=range, aes(fill=SCINAME), fill="gray10", alpha=0.3) +
  geom_point(data=pop, aes(x=Long, y=Lat), fill="white", colour="black", pch=21, size=6) +
  geom_sf_text(data=pop.sf, aes(label=order), nudge_y=0, nudge_x=0, size=3.5) +
  geom_point(data=pop, aes(x = Long, y=Lat, colour=factor(order)), size=-1) +
  xlab("") +
  ylab("") +
  xlim(c(-169, -52)) +
#  ylim(c(14, 85)) +
  scale_colour_manual(labels=paste0(pop.sf$order, " - ", pop.sf$Region),
                      values=rep("black", nrow(pop.sf)),
                      name="") +
  scale_fill_manual(labels=c("Breeding range"),
                    values=c("gray10"),
                    name="") +
  my.theme +
  theme(legend.position = "none",
        plot.margin = unit(c(0,1,-0.5,-0.5), "cm"),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10)) +
  theme(panel.grid.major = element_line(colour = "gray90"))
#studyarea

legend1map <- ggplot() +
  geom_sf(data=range, aes(fill=SCINAME), alpha=0.3) +
  scale_fill_manual(labels=c("Breeding range"),
                    values=c("gray10"),
                    name="") +
  my.theme +
  theme(legend.position=c(0.2, 2))
#legend1map

legend1 <- get_legend(legend1map)

legend2map <- ggplot() +
  geom_point(data=pop, aes(x = Long, y=Lat, colour=factor(order)), size=-1) +
  scale_colour_manual(labels=paste0(pop.sf$order, " - ", pop.sf$Region, " (", pop.sf$NTagsAnalysis,")"),
                      values=rep("black", nrow(pop.sf)),
                      name="Region\n(# tags transmitted)") +
  my.theme +
  theme(legend.position=c(0.45, 0.58),
        legend.text=element_text(size=10),
        legend.margin=margin(t = 0, unit='cm'))
legend2map

legend2 <- get_legend(legend2map)

ggplot2::ggsave(plot=grid.arrange(studyarea, legend1, legend2, 
                                  widths = c(6,1.75),
                                  heights = c(4, 0.5),
                                  layout_matrix = rbind(c(1,3),
                                                        c(1,2))),
                "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs/StudyArea.jpeg",
                width=7.5, height=5, units="in", device="jpeg")


#3. Figure 4 - Winter connectivity----
pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/tbl_population_abundance.csv") %>% 
  dplyr::filter(!Region %in% c("Quebec", "Florida")) %>% 
  dplyr::mutate(Region = case_when(Region=="BC coast" ~ "Coastal BC",
                                   Region=="BC Okanagan" ~ "Southcentral BC",
                                   !is.na(Region) ~ as.character(Region))) %>% 
  arrange(desc(Lat)) %>% 
  dplyr::select(Population, Region)

pop$Region <- factor(pop$Region, levels=c("Yukon", "Northwest Territories", "Alberta", "Saskatchewan", "Southcentral BC", "Coastal BC", "New Brunswick", "Ontario", "South Dakota", "Oregon", "Arizona", "Texas"))
classes <- length(unique(pop$Region))
clrs <- viridis::plasma(classes)

dat <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv")  %>% 
  mutate(Winter2 = case_when(PinpointID %in% c(81, 439, 443, 490, 825, 826, 828) & Season2=="Winter2" ~ 2,
                             PinpointID %in% c(81, 439, 443, 490, 825, 826, 828) & Season2=="Winter" ~ 3,
                             !PinpointID %in% c(81, 439, 443, 490, 825, 826, 828) & Season2=="Winter" ~ 1,
                             is.na(Season2) ~ 0)) %>% 
  dplyr::filter(PinpointID!=829) %>% 
  left_join(pop, by="Population") %>% 
  dplyr::select(Winter, Season, Season2, Region, PinpointID, BandLat, BandLong, Lat, Long) %>% 
  unique()

dat1 <- dat %>%
  dplyr::filter(Winter==1,
                Season=="Winter") %>% 
  dplyr::group_by(Region, PinpointID, BandLat, BandLong) %>% 
  summarize(WintLat = mean(Lat),
            WintLong = mean(Long)) %>% 
  ungroup() %>% 
  mutate(position=1)

dat2 <- dat %>%
  dplyr::filter(Winter==1,
                Season2=="Winter2") %>% 
  dplyr::group_by(Region, PinpointID, BandLat, BandLong) %>% 
  summarize(WintLat = mean(Lat),
            WintLong = mean(Long)) %>% 
  ungroup() %>% 
  mutate(position=2) %>% 
  rbind(dat1) %>% 
  arrange(PinpointID, position)

dat.id <- dat1 %>% 
  dplyr::select(Region, PinpointID) %>% 
  unique()
dat.source <- dat1 %>% 
  dplyr::select(PinpointID, BandLong, BandLat) %>% 
  dplyr::rename(Long=BandLong, Lat=BandLat) %>% 
  unique() %>% 
  dplyr::select(Long, Lat) %>% 
  as.matrix()
dat.dest <- dat1 %>% 
  dplyr::select(PinpointID, WintLong, WintLat) %>% 
  dplyr::rename(Long=WintLong, Lat=WintLat) %>% 
  unique() %>% 
  dplyr::select(Long, Lat) %>% 
  as.matrix()

#Great circle distance----
lines.gc <- vector(mode = "list", length = nrow(dat.id))
for (i in 1:nrow(dat.id)) {
  lines.gc[[i]] <- st_linestring(rbind(dat.source[i, ], dat.dest[i, ]))
}

lines.sf <- st_sfc(lines.gc, crs = 4326) %>% 
  st_segmentize(units::set_units(100, km)) 

routes.sf <- st_sf(dat.id, geometry = lines.sf)

plot.gcd <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  xlim(c(-170, -30)) +
  geom_sf(data = routes.sf, colour="grey50",
          alpha = 0.7,
          size=0.3) + 
  geom_point(aes(x = WintLong, y = WintLat,
                 colour=Region),
             data = dat1, 
             alpha = .6,
             size=2.5) +
  geom_point(aes(x = BandLong, y = BandLat,
                 colour=Region),
             data = dat1, 
             alpha = .6,
             size=2.5) +
  geom_rect(aes(xmin=(min(dat1$WintLong)-3),
                xmax=(max(dat1$WintLong)+3),
                ymin=(min(dat1$WintLat)-3),
                ymax=(max(dat1$WintLat)+3)),
            colour="black",
            fill="transparent",
            size=0.7) +
  geom_segment(aes(y = -10, x = -105, yend = -10, xend = min(dat1$WintLong)-3), size=0.7, colour="black") +
  scale_color_manual(values=clrs, name="") +
  my.theme +
  xlab("") +
  ylab("") +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_line(colour = "gray90")) +
  theme(legend.title = element_blank()) +
  guides(colour=guide_legend(nrow=3, ncol=4))
#plot.gcd

#Wintering ground inset----
#register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

center <- dat1 %>% 
  summarize(long=mean(WintLong),
            lat=mean(WintLat))

#map <- get_map(center, zoom=4, force=TRUE, maptype="satellite")

plot.wint <- ggmap(map, darken = c(1, "green")) +
  geom_point(aes(x = WintLong, y = WintLat,
                 fill=Region),
             shape = 21,
             colour="grey75",
             data = dat2, 
             alpha = 1,
             size=3) +
  geom_path(aes(x=WintLong, y=WintLat, group=PinpointID),
            colour="grey75",
            alpha = 0.7,
            size=0.5,
            data=dat2,
            arrow=arrow(length=unit(0.1, "inches"), type="closed")) +
  scale_fill_manual(values=clrs, name="") + 
  my.theme +
  coord_sf(datum = NA) +
  xlab("") +
  ylab("") +
  xlim(c(min(dat2$WintLong)-1, max(dat2$WintLong)+1)) +
  ylim(c(min(dat2$WintLat)-1, max(dat2$WintLat)+1)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        legend.position = "none",
        plot.background = element_rect(fill="transparent",
                                         colour =NA),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5))
plot.wint

#Put together----
plot.wintmc <- ggdraw() +
  draw_plot(plot.gcd) +
  draw_plot(plot.wint, 0.11, 0.185, 0.43, 0.43)

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
ggplot2::ggsave(plot=plot.wintmc, "WinteringConnectivity.jpeg", width=7.5, height=8, units="in", device="jpeg")


#4. Figure 3 - Connectivity Profiles----

my.theme <- theme_classic() +
  theme(text=element_text(size=12, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(size=12, hjust = 0.5))

profiles.theme <- my.theme +
  theme(text = element_text(size=10, family="Arial"),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10),
        legend.margin=margin(t = 0, unit='cm'),
        legend.spacing.x = unit(0.05, "cm"),
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(0.5, "cm"))

count <- 3

pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/tbl_population_abundance.csv") %>% 
  dplyr::filter(!Region %in% c("Quebec", "Florida")) %>% 
  dplyr::mutate(Region = case_when(Region=="BC coast" ~ "Coastal BC",
                                   Region=="BC Okanagan" ~ "Southcentral BC",
                                   !is.na(Region) ~ as.character(Region))) %>%
  arrange(desc(Lat)) %>% 
  dplyr::select(Population, Region) %>% 
  dplyr::rename(population = Population)

pop$Region <- factor(pop$Region, levels=c("Yukon", "Northwest Territories", "Alberta", "Saskatchewan", "Southcentral BC", "Coastal BC", "New Brunswick", "Ontario", "South Dakota", "Oregon", "Arizona", "Texas"))
classes <- length(unique(pop$Region))
clrs <- viridis::plasma(classes)

#A. Fall spatial connectivity map----

fall.lats <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_FallSpatial.csv") %>% 
  dplyr::filter(count.ind >= count) 

#Map
fs.mn <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Fall_Spatial_Mean.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(population = pop) %>% 
  left_join(pop, by="population") %>% 
  arrange(Region)

individuals <- fs.mn %>% 
  data.frame() %>% 
  dplyr::select(deployid, population) %>% 
  unique() %>% 
  group_by(population) %>% 
  summarize(inds=n()) %>% 
  ungroup() %>% 
  left_join(pop) %>% 
  mutate(label = paste0(Region, " (", inds, ")"))
individuals

clrs.fs <- clrs[sort(as.numeric(unique(individuals$Region)))]

plot.fallspat.map <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  xlim(c(-170, -30)) +
  ylim(c(-55, 75)) +
  geom_sf(data=fs.mn, size = 0.8, aes(color=Region), show.legend = "line") +
  scale_colour_manual(values=clrs.fs, name="") +
  xlab("") +
  ylab("") +
  profiles.theme +
  theme(legend.position="none") +
  theme(panel.grid.major = element_line(colour = "gray90")) + 
  theme(plot.margin = unit(c(-0.5,0,-1.5,-0.5), "cm"),
        plot.title=element_text(size=16, hjust = 0.5),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10)) +
  ggtitle("Fall migration")
plot.fallspat.map

#Inset profile
fall.spatial <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_Spatial_LeaveOneOut_MantelV3.csv") %>% 
  dplyr::rename(population = pop) %>% 
  dplyr::filter(lat %in% fall.lats$latr) %>% 
  dplyr::filter(!is.na(M))

fall.spatial.99 <- fall.spatial %>% 
  dplyr::filter(population==99)

fall.spatial.sum <- fall.spatial %>% 
  group_by(lat, population) %>% 
  summarize(Mmean=mean(M),
            Mlow83=quantile(M, probs=0.083),
            Mhigh83=quantile(M, probs=0.917),
            Mlow95=quantile(M, probs=0.025),
            Mhigh95=quantile(M, probs=0.975),
            Mlowsd=mean(M)-sd(M),
            Mhighsd=mean(M)+sd(M))  %>% 
  left_join(pop)

fall.spatial.sum.99 <- fall.spatial.sum %>% 
  dplyr::filter(population==99)

fall.spatial.gam <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_Spatial_GAM_Mantel.csv") %>% 
  dplyr::filter(lat %in% fall.lats$latr) %>% 
  dplyr::rename(population = pop)

fall.spatial.gam.99 <- fall.spatial.gam %>% 
  dplyr::filter(population==99)

fall.spatial.gam.pops <- fall.spatial.gam %>% 
  dplyr::filter(population!=99)  %>% 
  left_join(pop)

plot.fallspat <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=fall.spatial.sum.99) +
  geom_line(aes(x=lat, y=fit), data=fall.spatial.gam.99, colour="black") +
  scale_x_continuous(breaks=c(0,20,40)) +
  scale_y_continuous(limits = c(-1,1), breaks=c(-1,0,1.0), labels=c(-1, 0, 1)) +
  coord_flip() +
  ylab("Spatial\nconnectivity") +
  xlab("") +
  profiles.theme +
  theme(plot.margin = unit(c(0,0,0,-0.8), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.y.left = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x.bottom = element_text(size=10))

#Put together
grob.fallspat = ggplotGrob(plot.fallspat)
mapspatfall = plot.fallspat.map + 
  annotation_custom(grob = grob.fallspat, xmin=-175, xmax=-135, ymin=-40, ymax=61)
mapspatfall

#B. Fall spatial connectivity profile----
plot.fallspatmain <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=fall.spatial.sum.99) +
  geom_line(aes(x=lat, y=fit, colour=Region), size = 0.5, data=fall.spatial.gam.pops) +
  geom_line(aes(x=lat, y=fit), colour="black", size=1,  data=fall.spatial.gam.99) +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  scale_x_reverse() +
  scale_y_continuous(limits = c(-1,1)) +
  scale_fill_continuous(low = "grey80", high = "grey20") +
  scale_colour_manual(values=clrs.fs, name="Population left out\n(# individuals)",
                      labels=c("Northwest Territories (2)",
                               "Alberta (8)",
                               "Saskatchewan (2)",
                               "Southcentral BC (1)",
                               "Coastal BC (1)",
                               "Ontario (2)",
                               "South Dakota (2)",
                               "Oregon (1)",
                               "Arizona (1)",
                               "Texas (2)")) +
  ylab("Spatial connectivity") +
  xlab("Latitude") +
  profiles.theme +
  theme(plot.margin = unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10))
plot.fallspatmain


#C. Fall temporal connectivity profile----
fall.doys <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_FallTemporal.csv") %>% 
  dplyr::filter(count.ind >= count) 

fall.temporal <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_Temporal_LeaveOneOut_MantelV3.csv") %>%
  dplyr::filter(doy %in% fall.doys$doy) %>% 
  dplyr::filter(!is.na(M))

fall.temporal.99 <- fall.temporal %>% 
  dplyr::filter(day==99)

fall.temporal.sum <- fall.temporal %>% 
  group_by(doy, day) %>% 
  summarize(Mmean=mean(M),
            Mlowq=quantile(M, probs=0.083),
            Mhighq=quantile(M, probs=0.917),
            Mlow83=quantile(M, probs=0.083),
            Mhigh83=quantile(M, probs=0.917),
            Mlow95=quantile(M, probs=0.025),
            Mhigh95=quantile(M, probs=0.975),
            Mlowsd=mean(M)-sd(M),
            Mhighsd=mean(M)+sd(M)) %>% 
  ungroup()

fall.temporal.sum.99 <- fall.temporal.sum %>% 
  dplyr::filter(day==99)

fall.temporal.gam <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_temporal_GAM_Mantel.csv") %>% 
  dplyr::filter(doy %in% fall.doys$doy)

fall.temporal.gam.99 <- fall.temporal.gam %>% 
  dplyr::filter(day==99)

fall.temporal.gam.pops <- fall.temporal.gam %>% 
  dplyr::filter(day!=99)

fall.temporal.inds <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MovementModelPredictions_FallTemporal.csv")

individuals <- fall.temporal.inds %>% 
  dplyr::select(deployid, firstdoy, lastdoy) %>% 
  unique() %>% 
  gather(key=firstlast, value=doy, firstdoy:lastdoy) %>%
  group_by(doy) %>% 
  summarize(inds=n()) %>% 
  ungroup()
individuals

unique(fall.temporal$day)

plot.falltempmain <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlowq, ymax=Mhighq), alpha=0.4, data=fall.temporal.sum.99) +
  geom_line(aes(x=doy, y=fit, colour=factor(day)), size = 0.5, data=fall.temporal.gam.pops) +
  geom_line(aes(x=doy, y=fit), colour="black", size=1,  data=fall.temporal.gam.99) +
  geom_vline(aes(xintercept=297), linetype="longdash") +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  scale_x_continuous(breaks=c(243, 273, 304), labels=c("Sept", "Oct", "Nov")) +
  scale_y_continuous(limits = c(-1,1)) +
  scale_fill_continuous(low = "grey80", high = "grey20") +
  scale_colour_viridis_d(name="Start/end date left out\n(# individuals)",, option="D",
                         labels=c("Sept 10 (6)", "Sept 20 (3)", "Oct 21 (2)", "Oct 26 (6)", "Nov 11 (1)", "Nov 22")) +
  ylab("Temporal connectivity") +
  xlab("Day of year") +
  profiles.theme +
  theme(plot.margin = unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.margin=margin(t = 0, unit='cm'),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10))
plot.falltempmain


#D. Spring spatial connectivity map----
spring.lats <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_SpringSpatial.csv") %>% 
  dplyr::filter(count.ind >= count)

#Map
ss.mn <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Mean.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(population = pop) %>% 
  left_join(pop, by="population") %>% 
  arrange(Region)

individuals <- ss.mn %>% 
  data.frame() %>% 
  dplyr::select(deployid, population) %>% 
  unique() %>% 
  group_by(population) %>% 
  summarize(inds=n()) %>% 
  ungroup() %>% 
  left_join(pop) %>% 
  mutate(label = paste0(Region, " (", inds, ")"))
individuals

clrs.fs <- clrs[sort(as.numeric(unique(individuals$Region)))]

plot.springspat.map <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  xlim(c(-170, -30)) +
  ylim(c(-55, 75)) +
  geom_sf(data=ss.mn, size = 0.8, aes(color=Region), show.legend = "line") +
  scale_colour_manual(values=clrs, name="", drop=FALSE) +
  xlab("") +
  ylab("") +
  profiles.theme +
  theme(legend.position="none") +
  theme(panel.grid.major = element_line(colour = "gray90")) + 
  theme(plot.margin = unit(c(-0.5,0,-1.5,-0.5), "cm"),
        plot.title=element_text(size=16, hjust = 0.5),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10)) +
  ggtitle("Spring migration")
plot.springspat.map

#Inset profile
spring.spatial <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Spatial_LeaveOneOut_MantelV3.csv") %>% 
  dplyr::rename(population = pop) %>% 
  dplyr::filter(lat %in% spring.lats$latr) %>% 
  dplyr::filter(!is.na(M))

spring.spatial.99 <- spring.spatial %>% 
  dplyr::filter(population==99)

spring.spatial.sum <- spring.spatial %>% 
  group_by(lat, population) %>% 
  summarize(Mmean=mean(M),
            Mlowq=quantile(M, probs=0.083),
            Mhighq=quantile(M, probs=0.917),
            Mlowsd=mean(M)-sd(M),
            Mhighsd=mean(M)+sd(M))  %>% 
  left_join(pop)

spring.spatial.sum.99 <- spring.spatial.sum %>% 
  dplyr::filter(population==99) %>% 
  dplyr::filter(lat %in% spring.lats$latr)

spring.spatial.gam <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Spatial_GAM_Mantel.csv") %>% 
  dplyr::rename(population = pop)

spring.spatial.gam.99 <- spring.spatial.gam %>% 
  dplyr::filter(population==99)

spring.spatial.gam.pops <- spring.spatial.gam %>% 
  dplyr::filter(population!=99)  %>% 
  left_join(pop)

plot.springspat <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlowq, ymax=Mhighq), alpha=0.4, data=spring.spatial.sum.99) +
  geom_line(aes(x=lat, y=fit), data=spring.spatial.gam.99, colour="black") +
  scale_x_continuous(breaks=c(0,20,40)) +
  scale_y_continuous(limits = c(-1,1), breaks=c(-1,0,1.0), labels=c(-1, 0, 1)) +
  coord_flip() +
  ylab("Spatial\nconnectivity") +
  xlab("") +
  profiles.theme +
  theme(plot.margin = unit(c(0,0,0,-0.8), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.y.left = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x.bottom = element_text(size=10))

#Put together
grob.springspat = ggplotGrob(plot.springspat)
mapspatspring = plot.springspat.map + 
  annotation_custom(grob = grob.springspat, xmin=-175, xmax=-135, ymin=-37, ymax=60.5) +
  geom_segment(aes(y = 6, x = -152, yend = 6, xend = -78), linetype="longdash") +
  geom_segment(aes(y = 43.5, x = -145, yend = 43.5, xend = -124), linetype="longdash")
mapspatspring

annotation_custom(grob = grob.fallspat, xmin=-178, xmax=-135, ymin=-40, ymax=61)

#E. Spring spatial connectivity profile----
plot.springspatmain <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlowq, ymax=Mhighq), alpha=0.4, data=spring.spatial.sum.99) +
  geom_line(aes(x=lat, y=fit, colour=Region), size = 0.5, data=spring.spatial.gam.pops) +
  geom_line(aes(x=lat, y=fit), colour="black", size=1,  data=spring.spatial.gam.99) +
  geom_vline(aes(xintercept=c(6, 44)), linetype="longdash") +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  scale_x_continuous() +
  scale_y_continuous(limits = c(-1,1)) +
  scale_fill_continuous(low = "grey80", high = "grey20") +
  scale_colour_manual(values=clrs.fs, name="Population left out\n(# individuals)",
                      labels=c("Yukon (1)",
                               "Northwest Territories (1)",
                               "Alberta (8)",
                               "Saskatchewan (3)",
                               "Southcentral BC (2)",
                               "Coastal BC (1)",
                               "New Brunswick (1)",
                               "Ontario (3)",
                               "South Dakota (1)",
                               "Oregon (1)",
                               "Arizona (1)",
                               "Texas (1)")) +
  ylab("") +
  xlab("Latitude") +
  profiles.theme +
  theme(plot.margin = unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.margin=margin(t = 0, unit='cm'),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10))
plot.springspatmain

#F. Spring temporal connectivity profile----
spring.doys <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_SpringTemporal.csv") %>% 
  dplyr::filter(count.ind >= count) 

spring.temporal <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Temporal_LeaveOneOut_MantelV3.csv") %>% 
  dplyr::filter(doy %in% spring.doys$doy,
                !day %in% c(116, 119))

spring.temporal.99 <- spring.temporal %>% 
  dplyr::filter(day==99)

spring.temporal.sum <- spring.temporal %>% 
  group_by(doy, day) %>% 
  summarize(Mmean=mean(M),
            Mlow83=quantile(M, probs=0.083),
            Mhigh83=quantile(M, probs=0.917),
            Mlowsd=mean(M)-sd(M),
            Mhighsd=mean(M)+sd(M))

spring.temporal.sum.99 <- spring.temporal.sum %>% 
  dplyr::filter(day==99)

spring.temporal.gam <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_temporal_GAM_Mantel.csv") %>% 
  dplyr::filter(doy %in% spring.doys$doy,
                !day %in% c(116, 119))

spring.temporal.gam.99 <- spring.temporal.gam %>% 
  dplyr::filter(day==99)

spring.temporal.gam.pops <- spring.temporal.gam %>% 
  dplyr::filter(day!=99)

spring.temporal.inds <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MovementModelPredictions_SpringTemporal.csv")

individuals <- spring.temporal.inds %>% 
  dplyr::select(deployid, firstdoy, lastdoy) %>% 
  unique() %>% 
  gather(key=firstlast, value=doy, firstdoy:lastdoy) %>%
  group_by(doy) %>% 
  summarize(inds=n()) %>% 
  ungroup()
individuals

unique(spring.temporal$day)

plot.springtempmain <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=spring.temporal.sum.99) +
  geom_line(aes(x=doy, y=fit, colour=factor(day)), size = 0.5, data=spring.temporal.gam.pops) +
  geom_line(aes(x=doy, y=fit), colour="black", size=1,  data=spring.temporal.gam.99) +
  geom_vline(aes(xintercept=c(127)), linetype="longdash") +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  scale_x_continuous(breaks=c(60, 91, 121, 152), labels=c("March", "April", "May", "June")) +
  scale_y_continuous(limits = c(-1,1)) +
  scale_fill_continuous(low = "grey80", high = "grey20") +
  scale_colour_viridis_d(name="Start/end date left out\n(# individuals)", option="D",
                         labels=c("Mar 10 (1)", "Mar 17 (1)", "Mar 28 (2)", "Apr 1 (2)", "Apr 11 (6)", "Apr 28 (1)", "May 17 (6)", "May 27 (10)")) +
  ylab("") +
  xlab("Day of year") +
  profiles.theme +
  theme(plot.margin = unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.margin=margin(t = 0, unit='cm'),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10))
plot.springtempmain

#H. Put together----

profiles <- grid.arrange(plot.fallspatmain, plot.springspatmain,
                         mapspatfall, mapspatspring,
                         plot.falltempmain, plot.springtempmain,
                         widths = c(8,8),
                         heights = c(7, 4, 4),
                         layout_matrix = rbind(c(3,4),
                                               c(1,2),
                                               c(5,6)))

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
ggplot2::ggsave(plot=profiles,
                "ConnectivityProfiles.jpeg", width=11, height=12, units="in", device="jpeg")


setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
ggplot2::ggsave(plot=grid.arrange(plot.fallspatmain, plot.springspatmain,
                                  mapspatfall, mapspatspring,
                                  plot.falltempmain, plot.springtempmain,
                                  widths = c(8,8),
                                  heights = c(7, 4, 4),
                                  layout_matrix = rbind(c(3,4),
                                                        c(1,2),
                                                        c(5,6))),
                "ConnectivityProfiles.jpeg", width=11, height=12, units="in", device="jpeg")


#5. Figure 5 - Position of birds at temporal peaks----

#A. Fall migration day 297----

#Inset profile
fall.temp <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=fall.temporal.sum.99) +
  geom_line(aes(x=doy, y=fit), data=fall.temporal.gam.99, colour="black") +
  scale_x_reverse(breaks=c(243, 273, 304), labels=c("Sept", "Oct", "Nov")) +
  scale_y_continuous(breaks=c(-1, 0, 1), limits=c(-1, 1)) +
  coord_flip() +
  ylab("Temporal\nconnectivity") +
  xlab("") +
  my.theme +
  theme(plot.margin = unit(c(0,0,0,-0.8), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size=12))
fall.temp

#Map
fall.peaks <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/FallTemporalConnectivityPeaks_Mantel.csv") %>% 
  left_join(pop)

fall.peaks.errorup <- fall.peaks %>% 
  st_as_sf(coords=c("up.x", "up.y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(up.x=X, up.y=Y) 

fall.peaks.errorlw <- fall.peaks %>% 
  st_as_sf(coords=c("lw.x", "lw.y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(lw.x=X, lw.y=Y)

fall.peaks.error <- fall.peaks %>% 
  dplyr::select(X, Y) %>% 
  cbind(fall.peaks.errorup, fall.peaks.errorlw)
  
clrs.peak <- pop %>% 
  dplyr::select(Region) %>% 
  unique() %>% 
  cbind(clrs) %>% 
  dplyr::filter(Region %in% fall.peaks$Region)
clrs.peak <- as.character(clrs.peak$clrs)

plot.fall.peaks <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_point(aes(x=X, y=Y,
                 fill=Region,
                 colour=Region),
             shape = 21,
             alpha = 0.6,
             size=4,
            data = fall.peaks) +
  geom_segment(data=fall.peaks.error, aes(x=lw.x, xend=up.x, y=Y, yend=Y)) +
  geom_segment(data=fall.peaks.error, aes(x=X, xend=X, y=lw.y, yend=up.y)) +
  geom_rect(aes(xmin=(min(fall.peaks$X)-5),
                xmax=(max(fall.peaks$X)+5),
                ymin=(min(fall.peaks$Y)-5),
                ymax=(max(fall.peaks$Y)+5)),
            colour="black",
            fill="transparent",
            size=0.5,
            linetype="longdash") +
  scale_color_manual(values=clrs.peak, name="") +
  scale_fill_manual(values=clrs.peak, name="") +
  my.theme +
  xlab("") +
  ylab("") +
  theme(legend.pos = "none") +
  theme(panel.grid.major = element_line(colour = "gray90")) +
  theme(plot.background = element_rect(colour="transparent", fill="transparent")) +
  theme(plot.margin = unit(c(-0.5,0,-1.5,-0.5), "cm"),
        plot.title=element_text(size=16, hjust = 0.5)) +
  ggtitle("Fall migration - October 24") +
  coord_sf(xlim=c(-170, -30), ylim=c(-55, 75), expand = FALSE, crs=4326)
plot.fall.peaks

#Put together
grob.falltemp = ggplotGrob(fall.temp)
maptempfall = plot.fall.peaks + annotation_custom(grob = grob.falltemp, xmin=-165, xmax=-125, ymin=-46, ymax=43.5) +
  geom_segment(aes(y = -9.5, x = -139, yend = -9.5, xend = -90), linetype="longdash")
maptempfall


#B. Spring migration day 127----

#Inset profile
spring.temp <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=spring.temporal.sum.99) +
  geom_line(aes(x=doy, y=fit), data=spring.temporal.gam.99, colour="black") +
  scale_x_continuous(breaks=c(60, 91, 121, 152), labels=c("March", "April", "May", "June")) +
  scale_y_continuous(breaks=c(-1, 0, 1), limits=c(-1, 1)) +
  coord_flip() +
  ylab("Temporal\nconnectivity") +
  xlab("") +
  my.theme +
  theme(plot.margin = unit(c(0,0,0,-0.8), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size=12))
spring.temp

#Map
spring.peaks <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/SpringTemporalConnectivityPeaks_Mantel.csv") %>% 
  left_join(pop)

spring.peaks.errorup <- spring.peaks %>% 
  st_as_sf(coords=c("up.x", "up.y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(up.x=X, up.y=Y) 

spring.peaks.errorlw <- spring.peaks %>% 
  st_as_sf(coords=c("lw.x", "lw.y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(lw.x=X, lw.y=Y)

spring.peaks.error <- spring.peaks %>% 
  dplyr::select(X, Y) %>% 
  cbind(spring.peaks.errorup, spring.peaks.errorlw)

clrs.peak <- pop %>% 
  dplyr::select(Region) %>% 
  unique() %>% 
  cbind(clrs) %>% 
  dplyr::filter(Region %in% spring.peaks$Region)
clrs.peak <- as.character(clrs.peak$clrs)

plot.spring.peaks <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_point(aes(x=X, y=Y,
                 fill=Region,
                 colour=Region),
             shape = 21,
             alpha = 0.6,
             size=4,
             data = spring.peaks) +
  geom_segment(data=spring.peaks.error, aes(x=lw.x, xend=up.x, y=Y, yend=Y)) +
  geom_segment(data=spring.peaks.error, aes(x=X, xend=X, y=lw.y, yend=up.y)) +
  geom_rect(aes(xmin=(min(spring.peaks$X)-5),
                xmax=(max(spring.peaks$X)+5),
                ymin=(min(spring.peaks$Y)-5),
                ymax=(max(spring.peaks$Y)+5)),
            colour="black",
            fill="transparent",
            size=0.5,
            linetype="longdash") +
  scale_color_manual(values=clrs.peak, name="") +
  scale_fill_manual(values=clrs.peak, name="") +
  my.theme +
  xlab("") +
  ylab("") +
  theme(legend.position="none") +
  theme(panel.grid.major = element_line(colour = "gray90")) + 
  theme(plot.margin = unit(c(-0.5,0,-1.5,-0.5), "cm"),
        plot.title=element_text(size=16, hjust = 0.5)) +
  coord_sf(xlim=c(-170, -30), ylim=c(-55, 75), expand = FALSE, crs=4326) +
  ggtitle("Spring migration - May 7")
plot.spring.peaks


grob.springtemp = ggplotGrob(spring.temp)
maptempspring = plot.spring.peaks + annotation_custom(grob = grob.springtemp, xmin=-165, xmax=-125, ymin=-46, ymax=43.5) +
  geom_segment(aes(y = 18.5, x = -138, yend = 18.5, xend = -102), linetype="longdash")
maptempspring

#D. Legend----
temp.peaks <- rbind(fall.peaks, spring.peaks)

clrs.peak <- pop %>% 
  dplyr::select(Region) %>% 
  unique() %>% 
  cbind(clrs) %>% 
  dplyr::filter(Region %in% temp.peaks$Region)
clrs.peak <- as.character(clrs.peak$clrs)

plot.temp.legend <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  xlim(c(-170, -30)) +
  ylim(c(-55, 75)) +
  geom_point(aes(x=X, y=Y,
                 fill=Region,
                 colour=Region),
             shape = 21,
             alpha = 0.6,
             size=4,
             data = temp.peaks) +
  geom_rect(aes(xmin=(min(temp.peaks$X)-5),
                xmax=(max(temp.peaks$X)+5),
                ymin=(min(temp.peaks$Y)-5),
                ymax=(max(temp.peaks$Y)+5)),
            colour="black",
            fill="transparent",
            size=0.5,
            linetype="longdash") +
  scale_color_manual(values=clrs.peak, name="") +
  scale_fill_manual(values=clrs.peak, name="") +
  my.theme +
  xlab("") +
  ylab("") +
  theme(panel.grid.major = element_line(colour = "gray90")) +
  theme(plot.background = element_rect(colour="transparent", fill="transparent")) +
  theme(plot.margin = unit(c(0,0,-1.5,-0.5), "cm"),
        legend.position = "bottom") +
  guides(col = guide_legend(nrow = 3))
#plot.temp.legend

legend.temp <- get_legend(plot.temp.legend)

#E. Put together----
temporallocs <- grid.arrange(maptempfall, 
                             maptempspring, 
                             legend.temp,
                         widths = c(8),
                         heights = c(8,8,1),
                         layout_matrix = rbind(c(1),
                                               c(2),
                                               c(3)))

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
ggplot2::ggsave(plot=temporallocs, "Fig4TemporalPeakLocations.jpeg", width=7.5, height=16, units="in", device="jpeg")



#6. Figure 6: Simulation----

sim.theme <- theme_classic() +
  theme(text=element_text(size=10, family="Arial"),
        axis.text.x.bottom =element_text(size=10),
        axis.text.y.left=element_text(size=10),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        plot.title=element_text(size=10, hjust = 0.5))


#Step 1: Simulate random points----

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs/Simulation")
ptslow <- read.csv("SimulationStartingPoints_3-3-10_1.csv") %>% 
  sf::st_as_sf(coords = c("X","Y")) %>% 
  sf::st_set_crs(3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  dplyr::rename(lat=Y, long=X) %>% 
  cbind(read.csv("SimulationStartingPoints_3-3-10_1.csv")) %>% 
  mutate(pop=as.factor(pop))

boxes <- data.frame(xmin=c(-123, -120, -115, -105),
                    xmax=c(-96, -96, -96, -98),
                    ymin=c(50, 40, 30, 20),
                    ymax=c(65, 45, 35, 25),
                    colour=c("black", "grey20", "grey50", "grey80"),
                    label=c("Breeding", "Stopover 1", "Stopover 2", "Nonbreeding"))

segments <- data.frame(x=c(-96, -96, -96, -98),
                       xend=c(-80, -80, -80, -80),
                       y=c(57.5, 45, 32.5, 20),
                       yend=c(57.5, 45, 32.5, 20))

annotations <- data.frame(xmin=segments$xend + 20,
                          y=segments$y,
                          label=c("Breeding\n(high connectivity)",
                                  "Stopover 1\n(low connectivity)",
                                  "Stopover 2\n(high connectivity)",
                                  "Nonbreeding\n(low connectivity)"))

step1low <- ggplot() +
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_rect(data=boxes,
            aes(xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax),
            fill="transparent",
            colour="black",
            size=0.7) +
  geom_segment(data=segments, aes(x=x, xend=xend, y=y, yend=yend)) +
  geom_label(data=annotations, aes(x=xmin, y=y, label=label), nudge_x = 0, size=3) +
  geom_point(aes(x = long, y = lat,
                 colour=pop),
             data = ptslow, 
             alpha = .6,
             size=2.5,
             show.legend = FALSE) +
  scale_colour_viridis_d() +
  my.theme +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        axis.line = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  theme(plot.margin = unit(c(0,0,0,-0.8), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size=12))
step1low

#Step 2. Predict from random points----
pred1mnlow <- st_read("Predictions1_Mean_3-3-10_1.shp")
pred1lowlow <- st_read("Predictions1_Lower_3-3-10_1.shp")
pred1uplow <- st_read("Predictions1_Upper_3-3-10_1.shp")

step2low <- ggplot()+
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_sf(data=pred1mnlow, size = 0.8, aes(color=pop)) +
  geom_sf(data=pred1lowlow, size = 0.3, colour="black") +
  geom_sf(data=pred1uplow, size = 0.3, colour="black") +
  geom_point(aes(x = long, y = lat,
                 colour=pop),
             data = ptslow, 
             alpha = .6,
             size=2.5) +
  scale_colour_viridis_d() +
  my.theme +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        axis.line = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  theme(plot.margin = unit(c(0,0,0,-0.8), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size=12))
step2low

#Step 3. Simulate from posterior----
simlow <- st_read("Simulation_Mean_3-3-10_1.shp")

step3low <- ggplot()+
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_sf(data=simlow, size = 0.8, aes(color=pop)) +
  geom_point(aes(x = long, y = lat,
                 colour=pop),
             data = ptslow, 
             alpha = .6,
             size=2.5) +
  scale_colour_viridis_d() +
  my.theme +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        axis.line = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  theme(plot.margin = unit(c(0,0,0,-0.8), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size=12))
step3low

#Step 4. Randomly select GPS points----
gpsspatlow <- read.csv("SampledPoints_Spatial_3-3-10_1.csv") %>% 
  sf::st_as_sf(coords = c("X","Y")) %>% 
  sf::st_set_crs(3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  dplyr::rename(lat=Y, long=X) %>% 
  cbind(read.csv("SampledPoints_Spatial_3-3-10_1.csv")) %>% 
  mutate(pop=as.factor(pop))

step4spatlow <- ggplot() +
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_point(aes(x = long, y = lat,
                 colour=pop),
             data = gpsspatlow, 
             alpha = .6,
             size=2.5) +
  scale_colour_viridis_d() +
  my.theme +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        axis.line = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  theme(plot.margin = unit(c(0,0,0,-0.8), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size=12))
step4spatlow

#5. Step 5. Predict mean migration pathway from GPS points---

pred2spatlow <- st_read("Predictions2_Spatial_Mean_3-3-10_1.shp")

step5spatlow <- ggplot()+
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_sf(data=pred2spatlow, size = 0.8, aes(color=pop)) +
  geom_point(aes(x = long, y = lat,
                 colour=pop),
             data = gpsspatlow, 
             alpha = .6,
             size=2.5) +
  scale_colour_viridis_d() +
  my.theme +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        axis.line = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  theme(plot.margin = unit(c(0,0,0,-0.8), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size=12))
step5spatlow

#Step 6. Create connectivity profiles----

#Spatial
mcspat <- read.csv("MC_Spatial_3-3-10_1.csv") %>% 
  dplyr::rename(M=mantel)

mcspat99 <- mcspat %>% 
  dplyr::filter(l1o==99)

mcspatsum <- mcspat %>% 
  group_by(lat, l1o) %>% 
  summarize(Mmean=mean(M),
            Mlow83=quantile(M, probs=0.083),
            Mhigh83=quantile(M, probs=0.917))

mcspatsum99 <- mcspatsum %>% 
  dplyr::filter(l1o==99)

mcspatgam <- read.csv("GAM_Spatial3-3-10_1.csv")

mcspatgam99 <- mcspatgam %>% 
  dplyr::filter(l1o==99)

mcspatgaml1o <- mcspatgam %>% 
  dplyr::filter(l1o!=99)

step6spatlw <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=mcspatsum99) +
  geom_line(aes(x=lat, y=fit, colour=factor(l1o)), size = 0.5, data=mcspatgaml1o) +
  geom_line(aes(x=lat, y=fit), colour="black", size=1,  data=mcspatgam99) +
  scale_x_reverse() +
  scale_fill_continuous(low = "grey80", high = "grey20") +
  scale_colour_viridis_d(option="plasma", name="Population\nleft out") +
  ylab("Spatial\nconnectivity") +
  xlab("Latitude") +
  sim.theme +
  theme(plot.margin = unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))
step6spatlw

#Temporal
mctemp <- read.csv("MC_Temporal_3-3-10_1.csv") %>% 
  dplyr::rename(M=mantel)

mctemp99 <- mctemp %>% 
  dplyr::filter(l1o==99)

mctempsum <- mctemp %>% 
  group_by(doy, l1o) %>% 
  summarize(Mmean=mean(M),
            Mlow83=quantile(M, probs=0.083),
            Mhigh83=quantile(M, probs=0.917))

mctempsum99 <- mctempsum %>% 
  dplyr::filter(l1o==99)

mctempgam <- read.csv("GAM_Temporal3-3-10_1.csv")

mctempgam99 <- mctempgam %>% 
  dplyr::filter(l1o==99)

mctempgaml1o <- mctempgam %>% 
  dplyr::filter(l1o!=99)

step6templw <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=mctempsum99) +
  geom_line(aes(x=doy, y=fit, colour=factor(l1o)), size = 0.5, data=mctempgaml1o) +
  geom_line(aes(x=doy, y=fit), colour="black", size=1,  data=mctempgam99) +
  scale_fill_continuous(low = "grey80", high = "grey20") +
  scale_colour_viridis_d(option="viridis", name="Start/end\ndayleft out") +
  ylab("Temporal\nconnectivity") +
  xlab("Day of year") +
  sim.theme +
  theme(plot.margin = unit(c(0,0,1,0), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))
step6templw

#Step 7. Identify and validate peaks----
spatpeaklines <- data.frame(x=c(27, 32),
                            linetype=c("solid", "dotted"),
                            label=c("True peak",
                                    "False peak"))

step7spatlw <- ggplot() +
  geom_ribbon(aes(x=lat, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=mcspatsum99) +
  geom_line(aes(x=lat, y=fit), colour="black", size=1,  data=mcspatgam99) +
  geom_vline(dat=spatpeaklines, aes(xintercept=x, linetype=linetype)) +
  scale_x_reverse() +
  scale_fill_continuous(low = "grey80", high = "grey20") +
  scale_colour_viridis_d(option="plasma", name="Population\nleft out") +
  scale_linetype(name="", labels=spatpeaklines$label) +
  ylab("Spatial\nconnectivity") +
  xlab("Latitude") +
  sim.theme +
  theme(plot.margin = unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))
step7spatlw

temppeaklines <- data.frame(x=c(39, 57),
                            linetype=c("solid", "dotted"),
                            label=c("True peak",
                                    "False peak"))

step7templw <- ggplot() +
  geom_ribbon(aes(x=doy, ymin=Mlow83, ymax=Mhigh83), alpha=0.4, data=mctempsum99) +
  geom_line(aes(x=doy, y=fit), colour="black", size=1,  data=mctempgam99) +
  geom_vline(dat=temppeaklines, aes(xintercept=x, linetype=linetype)) +
  scale_fill_continuous(low = "grey80", high = "grey20") +
  scale_linetype(name="", labels=temppeaklines$label) +
  ylab("Temporal\nconnectivity") +
  xlab("Day of year") +
  sim.theme +
  theme(plot.margin = unit(c(0,0,1,0), "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))
step7templw

#Make word grobs----
step1text <- grobTree(textGrob("Step 1.\nSelect random points\nfor four stages\nwith high and\nlow connectivity",
                               x=0.5,  y=0.5, just="centre",
                                     gp=gpar(fontsize=12)))

step2text <- grobTree(textGrob("Step 2.\nCreate movement\nmodel from random\npoints",
                               x=0.5,  y=0.5, just="centre",
                               gp=gpar(fontsize=12)))

step3text <- grobTree(textGrob("Step 3.\nSimulate migration\npathway from posterior\nof movement model",
                               x=0.5,  y=0.5, just="centre",
                               gp=gpar(fontsize=12)))

step4text <- grobTree(textGrob("Step 4.\nRandomly sample\nlocations from\nsimulated migration\npathway",
                               x=0.5,  y=0.5, just="centre",
                               gp=gpar(fontsize=12)))

step5text <- grobTree(textGrob("Step 5.\nPredict migration\npathway from\nlocations",
                               x=0.5,  y=0.5, just="centre",
                               gp=gpar(fontsize=12)))

step6text <- grobTree(textGrob("Step 6.\nEstimate connectivity\nprofile from\nmean pathway\nand spatial error",
                               x=0.5,  y=0.5, just="centre",
                               gp=gpar(fontsize=12)))

step7text <- grobTree(textGrob("Step 7.\nIdentify and\nvalidate peaks in\nconnectivity profile",
                               x=0.5,  y=0.5, just="centre",
                               gp=gpar(fontsize=12)))


#Put together----
setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
ggplot2::ggsave(plot=grid.arrange(step1text, step1low, 
                                  step2text, step2low,
                                  step3text, step3low,
                                  step4text, step4spatlow,
                                  step5text, step5spatlow,
                                  step6text, step6spatlw, step6templw,
                                  step7text, step7spatlw, step7templw,
                                  widths = c(3,6,3,6),
                                  heights = c(6,3,3,3,3,6),
                                  layout_matrix = rbind(c(1,2,9,10),
                                                        c(3,4,11,12),
                                                        c(3,4,11,13),
                                                        c(5,6,14,15),
                                                        c(5,6,14,16),
                                                        c(7,8,NA,NA))),
                "Fig6Simulation.jpeg", width=11, height=11, units="in", device="jpeg")

 #7. Appendix 2: CRAWL output----
pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/tbl_population_abundance.csv") %>% 
  dplyr::filter(!Region %in% c("Quebec", "Florida")) %>% 
  dplyr::mutate(Region = case_when(Region=="BC coast" ~ "Coastal BC",
                                   Region=="BC Okanagan" ~ "Southcentral BC",
                                   !is.na(Region) ~ as.character(Region))) %>% 
  arrange(desc(Lat)) %>% 
  dplyr::select(Population, Region)

pop$Region <- factor(pop$Region, levels=c("Yukon", "Northwest Territories", "Alberta", "Saskatchewan", "Southcentral BC", "Coastal BC", "New Brunswick", "Ontario", "South Dakota", "Oregon", "Arizona", "Texas"))
classes <- length(unique(pop$Region))
clrs <- viridis::plasma(classes)

#Fall spatial----
fs.mn <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Fall_Spatial_Mean.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

fs.up <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Fall_Spatial_Upper.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

fs.lw <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Fall_Spatial_Lower.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

fs.pt <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Fall_Spatial_Points.csv") %>% 
  sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(4326) %>% 
  dplyr::rename(Population = population) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region) %>% 
  mutate(deployid=as.factor(deployid))

clrs.crawl <- pop %>% 
  dplyr::select(Region) %>% 
  unique() %>% 
  cbind(clrs) %>% 
  dplyr::filter(Region %in% fs.pt$Region)
clrs.crawl <- as.character(clrs.crawl$clrs)

app4.fs <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75") +
  xlim(c(-170, -30)) +
  geom_sf(data=fs.up, size = 0.3, colour="black") +
  geom_sf(data=fs.lw, size = 0.3, colour="black") +
  geom_sf(data=fs.mn, size = 0.8, aes(color=Region), show.legend = "line") +
  geom_sf(data=fs.pt, size = 1, col="black") +
  facet_wrap(deployid~., ncol=5) +
  scale_color_manual(values=clrs.crawl, name="") +
  xlab("") +
  ylab("") +
  my.theme +
  ggtitle("Fall migration spatial models") +
  theme(plot.title=element_text(size=32, hjust=0.5),
        plot.margin = unit(c(1,-1,0,-1.5), "cm"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text=element_text(size=24),
        strip.text.x=element_text(size=24),
        axis.text=element_text(size=10)) +
  theme(panel.grid.major = element_line(colour = "gray50")) +
  theme(legend.key.size = unit(4,"line")) +
  guides(col = guide_legend(nrow = 2))

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
Cairo::CairoPNG(width=20, height=24, units="in", dpi=300,
                "CRAWLFallSpatial.jpeg")
app4.fs
dev.off()

#Fall temporal----
ft.mn <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Fall_Temporal_Mean.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

ft.up <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Fall_Temporal_Upper.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

ft.lw <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Fall_Temporal_Lower.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

ft.pt <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Fall_Temporal_Points.csv") %>% 
  sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(4326) %>% 
  dplyr::rename(Population = population) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region) %>% 
  mutate(deployid=as.factor(deployid))

clrs.crawl <- pop %>% 
  dplyr::select(Region) %>% 
  unique() %>% 
  cbind(clrs) %>% 
  dplyr::filter(Region %in% ft.pt$Region)
clrs.crawl <- as.character(clrs.crawl$clrs)


app4.ft <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75") +
  xlim(c(-170, -30)) +
  geom_sf(data=ft.up, size = 0.3, colour="black") +
  geom_sf(data=ft.lw, size = 0.3, colour="black") +
  geom_sf(data=ft.mn, size = 0.8, aes(color=Region), show.legend = "line") +
  geom_sf(data=ft.pt, size = 1, col="black") +
  facet_wrap(deployid~., ncol=5) +
  scale_color_manual(values=clrs.crawl, name="") +
  xlab("") +
  ylab("") +
  my.theme +
  ggtitle("Fall migration temporal models") +
  theme(plot.title=element_text(size=32, hjust=0.5),
        plot.margin = unit(c(1,-1,0,-1.5), "cm"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text=element_text(size=24),
        strip.text.x=element_text(size=24),
        axis.text=element_text(size=10)) +
  theme(panel.grid.major = element_line(colour = "gray50")) +
  theme(legend.key.size = unit(4,"line")) +
  guides(col = guide_legend(nrow = 2))

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
Cairo::CairoPNG(width=20, height=20, units="in", dpi=300,
                "CRAWLFallTemporal.jpeg")
app4.ft
dev.off()

#Spring spatial----
ss.mn <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Mean.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

ss.up <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Upper.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

ss.lw <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Lower.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

ss.pt <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Points.csv") %>% 
  sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(4326) %>% 
  dplyr::rename(Population = population) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region) %>% 
  mutate(deployid=as.factor(deployid))

clrs.crawl <- pop %>% 
  dplyr::select(Region) %>% 
  unique() %>% 
  cbind(clrs) %>% 
  dplyr::filter(Region %in% ss.pt$Region)
clrs.crawl <- as.character(clrs.crawl$clrs)

app4.fs <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75") +
  xlim(c(-170, -30)) +
  geom_sf(data=ss.up, size = 0.3, colour="black") +
  geom_sf(data=ss.lw, size = 0.3, colour="black") +
  geom_sf(data=ss.mn, size = 0.8, aes(color=Region), show.legend = "line") +
  geom_sf(data=ss.pt, size = 1, col="black") +
  facet_wrap(deployid~., ncol=5) +
  scale_color_manual(values=clrs.crawl, name="") +
  xlab("") +
  ylab("") +
  my.theme +
  ggtitle("Spring migration spatial models") +
  theme(plot.title=element_text(size=32, hjust=0.5),
        plot.margin = unit(c(1,-1,0,-1.5), "cm"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text=element_text(size=24),
        strip.text.x=element_text(size=24),
        axis.text=element_text(size=10)) +
  theme(panel.grid.major = element_line(colour = "gray50")) +
  theme(legend.key.size = unit(4,"line")) +
  guides(col = guide_legend(nrow = 2))

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
Cairo::CairoPNG(width=20, height=24, units="in", dpi=300,
                "CRAWLSpringSpatial.jpeg")
app4.fs
dev.off()

#Spring temporal----
st.mn <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Temporal_Mean.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

st.up <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Temporal_Upper.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

st.lw <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Temporal_Lower.shp") %>% 
  mutate(pop=as.numeric(as.character(pop))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region)

st.pt <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Temporal_Points.csv") %>% 
  sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(4326) %>% 
  dplyr::rename(Population = population) %>% 
  left_join(pop, by="Population") %>% 
  arrange(Region) %>% 
  mutate(deployid=as.factor(deployid))

clrs.crawl <- pop %>% 
  dplyr::select(Region) %>% 
  unique() %>% 
  cbind(clrs) %>% 
  dplyr::filter(Region %in% st.pt$Region)
clrs.crawl <- as.character(clrs.crawl$clrs)


app4.ft <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75") +
  xlim(c(-170, -30)) +
  geom_sf(data=st.up, size = 0.3, colour="black") +
  geom_sf(data=st.lw, size = 0.3, colour="black") +
  geom_sf(data=st.mn, size = 0.8, aes(color=Region), show.legend = "line") +
  geom_sf(data=st.pt, size = 1, col="black") +
  facet_wrap(deployid~., ncol=5) +
  scale_color_manual(values=clrs.crawl, name="") +
  xlab("") +
  ylab("") +
  my.theme +
  ggtitle("Spring migration temporal models") +
  theme(plot.title=element_text(size=32, hjust=0.5),
        plot.margin = unit(c(1,-1,0,-1.5), "cm"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text=element_text(size=24),
        strip.text.x=element_text(size=24),
        axis.text=element_text(size=10)) +
  theme(panel.grid.major = element_line(colour = "gray50")) +
  theme(legend.key.size = unit(4,"line")) +
  guides(col = guide_legend(nrow = 2))

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
Cairo::CairoPNG(width=20, height=24, units="in", dpi=300,
                "CRAWLSpringTemporal.jpeg")
app4.ft
dev.off()


#8. Appendix 4: Stripchart-----
pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/tbl_population_abundance.csv") %>% 
  arrange(desc(Lat)) %>% 
  mutate(order=row_number(),
         poporder=ifelse(nchar(order)==1, paste0("0", order), order)) %>% 
  dplyr::select(Population, Region, poporder, order) %>% 
  dplyr::mutate(Region = case_when(Region=="BC coast" ~ "Coastal BC",
                                   Region=="BC Okanagan" ~ "Southcentral BC",
                                   !is.na(Region) ~ as.character(Region))) %>% 
  mutate(label = paste0(order, " - ", Region))
  

dat <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  filter(Type != "Band") %>% 
  left_join(pop) %>% 
  dplyr::mutate(Season=ifelse(Season2 %in% c("Breed1", "Breed2"), "Breeding", as.character(Season2)),
                Season=ifelse(Season2 %in% c("Winter1", "Winter2", "WinterMig"), "Winter", as.character(Season)),
                ID=paste0(poporder, Region, PinpointID)) %>% 
  arrange(desc(ID)) %>% 
  group_by(ID) %>% 
  mutate(rowID=row_number(),
         max=max(rowID)) %>% 
  ungroup() %>% 
  mutate(ID = as.factor(ID),
         ID = factor(ID, levels = rev(levels(ID))))

dat$Season <- factor(dat$Season, levels=c("Breeding", "FallMig", "Winter", "SpringMig"),
                     labels=c("Breeding", "Fall Migration", "Winter", "Spring Migration"))


labels <- dat %>% 
  arrange(ID) %>% 
  dplyr::select(label, ID) %>% 
  unique() %>% 
  group_by(label) %>% 
  mutate(order=row_number(),
         max=max(order)) %>% 
  ungroup() %>% 
  mutate(label=ifelse(order==max, as.character(label), ""),
         roworder=row_number())

plot1 <- ggplot(dat) +
  geom_point(aes(x=doy.order, y=factor(ID), colour=Season), size=2) +
  ylab("") +
  xlab("") +
  my.theme +
  scale_colour_viridis_d(name="", option="D") +
  scale_y_discrete(labels=labels$label) +
  scale_x_continuous(breaks=c(22, 113, 203, 295),
                     labels=c("Sep", "Dec", "Mar", "June")) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y  = element_text(angle=0)) +
  theme(plot.margin = unit(c(0,0,-0.6,-0.6), "cm"),
        legend.position = "none")
plot1

legend <- get_legend(plot1 + theme(legend.position="bottom",
                                   legend.text=element_text(size=14)))

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs")
Cairo::CairoPNG(width=8, height=9, units="in", dpi=300,
                "Stripchart.png")
grid.arrange(plot1, legend, 
             widths = c(1),
             heights = c(1, 0.1),
             layout_matrix = rbind(c(1),
                                   c(2)))
dev.off()

#9. Summary stats----
setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data")
pop <- read.csv("tbl_population_abundance.csv") %>% 
  dplyr::filter(Region != "Florida") %>% 
  dplyr::mutate(Region = case_when(Region=="BC coast" ~ "Coastal BC",
                                   Region=="BC Okanagan" ~ "Southcentral BC",
                                   !is.na(Region) ~ as.character(Region))) %>% 
  arrange(-Lat) %>% 
  mutate(order=row_number()) 

dat <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv")

tags <- dat %>% 
  dplyr::select(PinpointID, Population) %>% 
  unique() %>% 
  left_join(pop)

table(tags$Region)

#Ratio of males to females
sex <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  dplyr::select(PinpointID, Population, Sex) %>% 
  unique() %>% 
  group_by(Sex) %>% 
  summarize(count = n())
sex

population <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  left_join(pop) %>% 
  dplyr::select(PinpointID, Region, Sex) %>% 
  unique() %>% 
  group_by(Region) %>% 
  summarize(count = n())
population

#Mean migration connectivity
fall.spatial <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_Spatial_LeaveOneOut_MantelV4.csv") %>% 
  dplyr::rename(population = pop) %>% 
#  dplyr::filter(lat %in% fall.lats$latr) %>% 
  dplyr::filter(!is.na(M))

fall.spatial.sum.99 <- fall.spatial %>% 
  group_by(lat, population) %>% 
  summarize(Mmean=mean(M),
            Mlow83=quantile(M, probs=0.083),
            Mhigh83=quantile(M, probs=0.917),
            Mlow95=quantile(M, probs=0.025),
            Mhigh95=quantile(M, probs=0.975),
            Mlowsd=mean(M)-sd(M),
            Mhighsd=mean(M)+sd(M)) %>% 
  ungroup() %>% 
  dplyr::filter(population==99)

mean(fall.spatial.sum.99$Mmean)
mean(fall.spatial.sum.99$Mlow95)
mean(fall.spatial.sum.99$Mhigh95)
tail(fall.spatial.sum.99)
head(fall.spatial.sum.99)

spring.spatial <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Spatial_LeaveOneOut_MantelV4.csv") %>% 
  dplyr::rename(population = pop) %>% 
  dplyr::filter(lat %in% spring.lats$latr) %>% 
  dplyr::filter(!is.na(M))

spring.spatial.sum.99 <- spring.spatial %>% 
  group_by(lat, population) %>% 
  summarize(Mmean=mean(M),
            Mlow83=quantile(M, probs=0.083),
            Mhigh83=quantile(M, probs=0.917),
            Mlow95=quantile(M, probs=0.025),
            Mhigh95=quantile(M, probs=0.975),
            Mlowsd=mean(M)-sd(M),
            Mhighsd=mean(M)+sd(M)) %>% 
  ungroup() %>% 
  dplyr::filter(population==99)

spring.spatial.sum.99.50 <- spring.spatial.sum.99 %>% 
  dplyr::filter(Mmean > 0.5)

mean(spring.spatial.sum.99$Mmean)
mean(spring.spatial.sum.99$Mlow95)
mean(spring.spatial.sum.99$Mhigh95)
tail(spring.spatial.sum.99)
head(spring.spatial.sum.99)

fall.spatial.sum.99.50 <- fall.spatial.sum.99 %>%
  dplyr::filter(Mmean > 0.5)

fall.temporal <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Fall_Temporal_LeaveOneOut_MantelV4.csv") %>%
  dplyr::filter(!is.na(M))

fall.temporal.sum.99 <- fall.temporal %>% 
  group_by(doy, day) %>% 
  summarize(Mmean=mean(M),
            Mlow83=quantile(M, probs=0.083),
            Mhigh83=quantile(M, probs=0.917),
            Mlow95=quantile(M, probs=0.025),
            Mhigh95=quantile(M, probs=0.975),
            Mlowsd=mean(M)-sd(M),
            Mhighsd=mean(M)+sd(M)) %>% 
  ungroup() %>% 
  dplyr::filter(day==99)

fall.temporal.sum.99.50 <- fall.temporal.sum.99 %>% 
  dplyr::filter(Mmean > 0.5)

mean(fall.temporal.sum.99$Mmean)
mean(fall.temporal.sum.99$Mlow95)
mean(fall.temporal.sum.99$Mhigh95)
tail(fall.temporal.sum.99)
head(fall.temporal.sum.99)

spring.temporal.doys <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/PopulationCount_SpringTemporal.csv") %>% 
  dplyr::filter(count.ind >= count) 

spring.temp <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/MigratoryConnectivity_Spring_Temporal_LeaveOneOut_MantelV4.csv") %>%
  dplyr::filter(!is.na(M),
                doy %in% spring.temporal.doys$doy)

spring.temporal.sum.99 <- spring.temporal %>% 
  group_by(doy, day) %>% 
  summarize(Mmean=mean(M),
            Mlow83=quantile(M, probs=0.083),
            Mhigh83=quantile(M, probs=0.917),
            Mlow95=quantile(M, probs=0.025),
            Mhigh95=quantile(M, probs=0.975),
            Mlowsd=mean(M)-sd(M),
            Mhighsd=mean(M)+sd(M)) %>% 
  ungroup() %>% 
  dplyr::filter(day==99)

spring.temporal.sum.99.50 <- spring.temporal.sum.99 %>% dplyr::filter(Mmean > 0.5)

mean(spring.temporal.sum.99$Mmean)
mean(spring.temporal.sum.99$Mlow95)
mean(spring.temporal.sum.99$Mhigh95)
tail(spring.temporal.sum.99)
head(spring.temporal.sum.99)

#10. Appendix 3 - Differential Migration----
dat <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  mutate(DateTime = ymd_hms(DateTime),
         BandDate = ymd(BandDate),
         Year = year(DateTime),
         YearDep = year(BandDate)) %>% 
  mutate(Season=Season2)

#Winter----
dat.wint <- dat %>% 
  dplyr::filter(Season=="Winter") %>% 
  mutate(Sex=factor(Sex, levels=c("M", "F")))

plot.wint.sex <- ggplot() +
  geom_polygon(data=sam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_point(aes(x = Long, y = Lat,
                 colour=Sex),
             data = dat.wint, 
             alpha = .6,
             size=3) +
  scale_colour_viridis_d(option="E", labels=c("Male", "Female")) +
  my.theme +
  xlab("") +
  ylab("") +
  theme(panel.grid.major = element_line(colour = "gray90")) +
  theme(legend.title = element_blank()) +
  coord_sf(xlim=c(-90, -30), ylim=c(-60, 15), expand = FALSE, crs=4326)
plot.wint.sex

ggsave(plot=plot.wint.sex, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs/Appendix3_Winter.jpeg", width=6, height=6, units="in", device="jpeg")

#Spring migration spatial----
dat.sex <- dat %>% 
  dplyr::select(PinpointID, Sex) %>% 
  rename(deployid = PinpointID) %>% 
  unique() %>% 
  mutate(Sex=factor(Sex, levels=c("M", "F")))

ss.mn <- st_read("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Results/CRAWL_MigConnectivity_Spring_Spatial_Mean.shp") %>% 
  mutate(pop=as.numeric(as.character(pop)),
         deployid=as.numeric(as.character(deployid))) %>% 
  dplyr::rename(Population = pop) %>% 
  left_join(pop) %>% 
  arrange(Region) %>% 
  left_join(dat.sex)

plot.spring.sex.spat <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  xlim(c(-170, -30)) +
  ylim(c(-55, 75)) +
  geom_sf(data=ss.mn, size = 0.8, aes(color=Sex), show.legend = "line") +
  scale_colour_viridis_d(option="E", labels=c("Male", "Female")) +
  xlab("") +
  ylab("") +
  my.theme +
  theme(panel.grid.major = element_line(colour = "gray90")) +
  theme(legend.title = element_blank()) +
  theme(plot.margin = unit(c(0,0,-1.5,-0.5), "cm"),
        plot.title=element_text(size=24, hjust = 0.5))
plot.spring.sex.spat

ggsave(plot=plot.spring.sex.spat, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs/Appendix3_SpringSpat.jpeg", width=7, height=6, units="in", device="jpeg")

#Spring migration temporal----
pred <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/springmigratepred.csv")

plot.spring.sex.temp <- ggplot(pred) +
  geom_line(aes(x=doy, y=fit, colour=Sex)) +
  geom_ribbon(aes(x=doy, ymax=upr, ymin=lwr, colour=Sex), alpha=0.5) +
  geom_point(aes(x=doy, y=Lat, colour=factor(Sex))) +
  scale_colour_viridis_d(option="E", labels=c("Male", "Female")) +
  scale_x_continuous(breaks=c(91, 121, 152), labels=c("April 1", "May 1", "June 1")) +
  xlab("") +
  ylab("Latitude (degrees north)") +
  my.theme
plot.spring.sex.temp

ggsave(plot=plot.spring.sex.temp, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs/Appendix3_SpringTemp.jpeg", width=7, height=6, units="in", device="jpeg")

#11. Alternative figure - simulation results----

pred.ind.fn <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Connectivity/Simulation_MC5/SimulationResults_FN_Individual.csv") %>% 
  mutate(ind="fn",
         val="ind") %>% 
  rename(n = n.ind) %>% 
  dplyr::select(fit, se, n, con, ind, val) %>% 
  mutate(se=ifelse(con=="temporal", NA, se),
         fit=ifelse(con=="temporal", NA, fit))

pred.gps.fn <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Connectivity/Simulation_MC5/SimulationResults_FN_GPS.csv") %>% 
  mutate(ind="fn",
         val="gps") %>% 
  rename(n = n.gps) %>% 
  dplyr::select(fit, se, n, con, ind, val) %>% 
  dplyr::filter(n<=50)

pred.ind.fp <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Connectivity/Simulation_MC5/SimulationResults_FP_Individual.csv") %>% 
  mutate(ind="fp",
         val="ind") %>% 
  rename(n = n.ind) %>% 
  dplyr::select(fit, se, n, con, ind, val) %>% 
  mutate(se=ifelse(con=="temporal", NA, se),
         fit=ifelse(con=="temporal", NA, fit))

pred.gps.fp <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Connectivity/Simulation_MC5/SimulationResults_FP_GPS.csv") %>% 
  mutate(ind="fp",
         val="gps") %>% 
  rename(n = n.gps) %>% 
  dplyr::select(fit, se, n, con, ind, val) %>% 
  mutate(se=ifelse(con=="spatial", NA, se),
         fit=ifelse(con=="spatial", NA, fit)) %>% 
  dplyr::filter(n<=50)

pred <- rbind(pred.ind.fn, pred.gps.fn, pred.ind.fp, pred.gps.fp)

my.theme <- theme.classic() +
  theme(text = element_text(size=10, family="Arial"),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10),
        legend.margin=margin(t = 0, unit='cm'),
        legend.spacing.x = unit(0.05, "cm"),
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(0.5, "cm"))

plot.ind.fn <- ggplot(pred.ind.fn) +
  geom_ribbon(aes(x=n, ymax=fit+2*se, ymin=fit-2*se, group=con), alpha=0.5) +
  geom_line(aes(x=n, y=fit, colour=con)) +
  scale_colour_viridis_d(option="D", name="Migratory\nconnectivity\nprofile") +
  labs(x="", y="Probability of detecting a true peak") +
  my.theme +
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0,0,0,0), "cm")) +
  ylim(c(0.65, 1))

plot.ind.fp <- ggplot(pred.ind.fp) +
  geom_ribbon(aes(x=n, ymax=fit+2*se, ymin=fit-2*se, group=con), alpha=0.5) +
  geom_line(aes(x=n, y=fit, colour=con)) +
  scale_colour_viridis_d(option="D", name="Migratory\nconnectivity\nprofile") +
  labs(x="Individuals per population", y="Probability of detecting a false peak") +
  my.theme +
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0,0,0,0), "cm")) +
  ylim(c(0, 0.5))

plot.gps.fn <- ggplot(pred.gps.fn) +
  geom_ribbon(aes(x=n, ymax=fit+2*se, ymin=fit-2*se, group=con), alpha=0.5) +
  geom_line(aes(x=n, y=fit, colour=con)) +
  scale_colour_viridis_d(option="D", name="Migratory\nconnectivity\nprofile") +
  labs(x="", y="") +
  my.theme +
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0,0,0,0), "cm")) +
  ylim(c(0.65, 1))

plot.gps.fp <- ggplot(pred.gps.fp) +
  geom_ribbon(aes(x=n, ymax=fit+2*se, ymin=fit-2*se, group=con), alpha=0.5) +
  geom_line(aes(x=n, y=fit, colour=con)) +
  scale_colour_viridis_d(option="D", name="Migratory\nconnectivity\nprofile") +
  labs(x="Locations per individual", y="") +
  my.theme +
  theme(legend.position = "none") +
  theme(plot.margin=unit(c(0,0,0,0), "cm")) +
  ylim(c(0, 0.5))

plot.legend <- ggplot(pred.gps.fn) +
  geom_ribbon(aes(x=n, ymax=fit+2*se, ymin=fit-2*se, group=con), alpha=0.5) +
  geom_line(aes(x=n, y=fit, colour=con)) +
  scale_colour_viridis_d(option="D", name="", labels=c("Spatial connectivity", "Temporal connectivity")) +
  labs(x="Locations per individual", y="") +
  my.theme +
  theme(legend.position = "bottom",
        legend.text = element_text(size=12))
legend <- get_legend(plot.legend)

plot.sim1 <- grid.arrange(plot.ind.fn, plot.gps.fn, 
             plot.ind.fp, plot.gps.fp, legend,
             widths = c(2,2,2,2),
             heights = c(4,4,1),
             layout_matrix=rbind(c(1,1,2,2),
                                 c(3,3,4,4),
                                 c(NA,5,5,NA)))

ggsave(plot=plot.sim1, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Figs/SimulationResults.jpeg", width=6.5, height=7, units="in", device="jpeg", dpi=300)


plot.sim2 <- ggplot(pred) +
  geom_ribbon(aes(x=n, ymax=fit+2*se, ymin=fit-2*se, group=con), alpha=0.5) +
  geom_line(aes(x=n, y=fit, colour=con)) +
  scale_colour_viridis_d(option="E", name="Migratory\nconnectivity\nprofile") +
  facet_grid(ind~val, scales="free")
plot.sim2

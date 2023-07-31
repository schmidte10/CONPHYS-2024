#--- Loading packages ---#
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(tidyverse)
library(cowplot)
library(grid)
library(ggpubr)
library(ggrepel)
#--- Loading shape data ---#
GBR_data <- st_read("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Statistics_wrkshp/Maps/TS_AIMS_NESP_Torres_Strait_Features_V1b_with_GBR_Features")
GBR_data <- st_read("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Statistics_wrkshp/Maps/TS_AIMS_NESP_Torres_Strait_Features_V1b_with_GBR_Features")
class(GBR_data)

#--- set working directory ---# 
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/map")

#--- Filtering islands in Australia only ---#
# Sampling 500 random reefs from the dataset so it won't be too heavy
GBR_data_sub2 <- GBR_data %>% 
  filter(DATASET == "GBR Features" & 
           LEVEL_1 == "Reef" & 
           Country == "Australia" &
           GBR_ID %in% sample(unique(GBR_data$GBR_ID)) |   
           GBR_ID == "17001A" | 
           GBR_ID == "16044B" | 
           GBR_ID == "16026" |
           GBR_ID == "20308" | 
           QLD_NAME == "Arlington Reef") %>% 
  filter(SHAPE_AREA >= 3) 

myreefs.core <- GBR_data %>% 
  filter(DATASET == "GBR Features" & 
           LEVEL_1 == "Reef" & 
           Country == "Australia" & 
          QLD_NAME != "U/N Reef" & 
           QLD_NAME != "U/N Rock") %>% 
  filter(#QLD_NAME == "Sudbury Reef" | 
           #QLD_NAME == "Vlasoff Reef" | 
           #QLD_NAME == "Tongue Reef" | 
           QLD_NAME == "Arlington Reef" | 
           QLD_NAME == "Morinda Shoal" | 
           #QLD_NAME == "Maori Reef" | 
           #QLD_NAME == "Oyster Reef" | 
           QLD_NAME == "Otter Reef" | 
           QLD_NAME == "Satellite Reef") %>%
  add_row(QLD_NAME = "Pretty Patches Reef", X_COORD = 146.042, Y_COORD =-16.622)  %>% 
  add_row(QLD_NAME = "Russell Island", X_COORD = 146.094, Y_COORD =-17.229) %>% 
  #add_row(QLD_NAME = "Keswick Island", X_COORD = 149.406, Y_COORD =-20.908) %>% 
  add_row(QLD_NAME = "Cockermouth Island", X_COORD = 149.398, Y_COORD =-20.772) %>% 
  add_row(QLD_NAME = "Low Wooded Island", X_COORD = 145.38, Y_COORD =-15.094)
  
  # add locations 
city_names <- c("Cairns", "Cooktown", "Mackay", "Townsville")
cities <-  GBR_data %>% 
  filter(DATASET == "GBR Features" & 
           LEVEL_1 == "Reef" & 
           Country == "Australia" & 
           QLD_NAME != "U/N Reef" & 
           QLD_NAME != "U/N Rock") %>% 
  add_row(QLD_NAME = "Cooktown", X_COORD = 145.25050, Y_COORD =-15.46814)  %>% 
  add_row(QLD_NAME = "Cairns", X_COORD = 145.754120, Y_COORD =-16.925491) %>% 
  add_row(QLD_NAME = "Townsville", X_COORD = 146.816956, Y_COORD =-19.258965)  %>% 
  add_row(QLD_NAME = "Mackay", X_COORD = 149.186813, Y_COORD =-21.144337) %>% 
  filter(QLD_NAME %in% city_names)
  
  #add_row(QLD_NAME = "Chauvel Reef (Southern)", X_COORD = 150.363, Y_COORD =-20.863) #%>% 
  #add_row(QLD_NAME = "Low Wooded Island", X_COORD = 145.367, Y_COORD = -15.084)

#myreefs.core.sample.size <- myreefs.core %>% 
  #add_column(sample.size= c("1","1","1","4","2","1","5"))

# Putting the gbr dataset in the same coordinate system as the AUS dataset
world <- ne_countries(scale = "medium", returnclass = "sf")
AUS <- world %>% 
  filter(name_long == "Australia") 
coords <- sf::st_transform(GBR_data_sub2,
                           crs = st_crs(AUS)$proj4string)

map_aus <- ggplot(data = AUS) +
  geom_sf(size = 0.1, fill = "#333333") +
  coord_sf(expand = F, xlim = c(112,154), ylim = c(-45,-10)) +
  geom_rect(aes(xmin=142, xmax=153, ymin=-25, ymax=-10),fill="transparent", linetype = "dotted", color="firebrick") +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border= element_blank(),
        axis.text.x= element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

#--- Plot ---#
core.reefs <- ggplot() +
  geom_sf(data = coords,fill = "gray90", color = "darkblue", size = 0.5) +
  geom_sf(data = AUS,fill = "gray90", color = "#333333", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(143, 155), ylim = c(-22, -14), expand = F) +
  theme(panel.background = element_rect(fill = "lightblue")) +
  xlab("")+ylab("")+ #ggtitle("Core region reefs") +
  geom_point(data = myreefs.core, aes(x = X_COORD, y = Y_COORD), size = 5, 
             shape = 21, fill = "#2f3544") + 
  geom_label_repel(data = myreefs.core, aes(X_COORD, Y_COORD, label=QLD_NAME), 
                   fill = "white", 
                   nudge_x = 4.0, 
                   nudge_y = 1.0, 
                   max.iter = 5000)+
  geom_point(data = cities, aes(x = X_COORD, y = Y_COORD), size =4, shape = 18, color = "red") + 
  geom_label_repel(data = cities, aes(X_COORD, Y_COORD, label=QLD_NAME), 
                   fill = "gray90", 
                   nudge_x = -0.5); core.reefs #%>%
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        panel.border = element_rect(colour = "#2f3544", fill=NA, size=2),
        axis.ticks.y=element_blank()); core.reefs

#---final figure ---# 
pdf("population_map_HolsworthReefs.pdf")
core.reefs
dev.off() 

ggsave("population_map_HolsworthReefs.jpeg", plot = core.reefs, width = 10, height = 7, units = "in", dpi = 300)

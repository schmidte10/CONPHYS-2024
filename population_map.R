#--- Loading packages ---#
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(tidyverse)
library(cowplot)
#--- Loading shape data ---#
GBR_data <- st_read("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Statistics_wrkshp/Maps/TS_AIMS_NESP_Torres_Strait_Features_V1b_with_GBR_Features")
class(GBR_data)

#--- set working directory ---# 
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation")

#--- Filtering islands in Australia only ---#
#--- GBR data ---#
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
  filter(SHAPE_AREA >= 3) %>%
  add_row(QLD_NAME = "Cockermouth Island", X_COORD = 149.398, Y_COORD =-20.772)  %>% 
  add_row(QLD_NAME = "Keswick Island", X_COORD = 149.406, Y_COORD =-20.908)


myreefs.core <- GBR_data %>% 
  filter(DATASET == "GBR Features" & 
           LEVEL_1 == "Reef" & 
           Country == "Australia" & 
          QLD_NAME != "U/N Reef" & 
           QLD_NAME != "U/N Rock") %>% 
  filter(GBR_ID == "17001A" | 
           GBR_ID == "16044B" | 
           GBR_ID == "16026")  

myreefs.leading <- GBR_data %>% 
  filter(DATASET == "GBR Features" & 
           LEVEL_1 == "Reef" & 
           Country == "Australia" & 
           QLD_NAME != "U/N Reef" & 
           QLD_NAME != "U/N Rock") %>% 
  filter(GBR_ID == "20308")%>% 
  add_row(GBR_NAME = "Cockermouth Island", X_COORD = 149.398, Y_COORD =-20.772)  %>% 
  add_row(GBR_NAME = "Keswick Island", X_COORD = 149.406, Y_COORD =-20.908)



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
p1 <- ggplot() +
  geom_sf(data = coords,fill = "gray90", color = "darkblue", size = 0.5) +
  geom_sf(data = AUS,fill = "gray90", color = "#333333", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(142, 153), ylim = c(-25, -10), expand = F) + 
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "lightblue"), 
        panel.grid.major = element_line(colour = "lightblue")) + 
  geom_point(data = myreefs.core, aes(x = X_COORD, y = Y_COORD), size = 3, 
             shape = 21, fill = "#2f3544")+ 
  geom_point(data = myreefs.leading, aes(x = X_COORD, y = Y_COORD), size = 3, 
             shape = 21, fill = "#4c8494") + 
  geom_rect(aes(xmin=145, xmax=147, ymin=-17.5, ymax=-16),fill="transparent", linetype = "dotted", color="#2f3544") + 
  geom_rect(aes(xmin=148.5, xmax=151, ymin=-21.7, ymax=-20.2),fill="transparent", linetype = "dotted", color="#4c8494") + 
  annotate("text", x = 149.8, y = -14, label = "Coral \nSea", fontface = "italic", size = 6);p1 

insert1 <- ggplot() +
  geom_sf(data = coords,fill = "darkblue", color = "darkblue", size = 0.5) +
  geom_sf(data = AUS,fill = "gray90", color = "#333333", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(145, 147), ylim = c(-17.6, -16), expand = F) + 
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "lightblue")) + 
  annotate("segment", x = 145.7926, xend = 145.7926+0.3, y = -16.32093, yend = -16.32093+0.28, colour = "black", size = 1)+
  annotate("segment", x = 146.2049, xend = 146.2049+0.3, y = -16.99838, yend = -16.99838+0.28, colour = "black", size = 1)+
  annotate("segment", x = 145.9929, xend = 145.9929+0.3, y = -16.65520, yend = -16.65520+0.28, colour = "black", size = 1)+
  geom_point(data = myreefs.core, aes(x = X_COORD, y = Y_COORD), size = 6, 
             shape = 21, fill = "#4c8494") + 
  geom_text(data = myreefs.core, aes(X_COORD+.3, Y_COORD+.3, label=QLD_NAME)); insert1

insert2 <- ggplot() +
  geom_sf(data = coords,fill = "gray90", color = "darkblue", size = 0.5) +
  geom_sf(data = AUS,fill = "gray90", color = "#333333", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(148.5, 151), ylim = c(-21.7, -20.2), expand = F) + 
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "lightblue"), 
        panel.grid.major = element_line(colour = "lightblue")) +  
  geom_point(data = myreefs.leading, aes(x = X_COORD, y = Y_COORD), size = 3, 
             shape = 21, fill = "#2f3544")+ 
  geom_text(data = myreefs.leading, aes(X_COORD, Y_COORD-0.03, label=GBR_NAME)); insert2

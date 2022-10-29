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
class(GBR_data)

#--- set working directory ---# 
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/map/")

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
  filter(QLD_NAME == "Sudbury Reef" | 
           QLD_NAME == "Vlasoff Reef" | 
           QLD_NAME == "Tongue Reef" | 
           QLD_NAME == "Arlington Reef" | 
           QLD_NAME == "Morinda Shoal" | 
           QLD_NAME == "Maori Reef" | 
           QLD_NAME == "Oyster Reef" | 
           QLD_NAME == "Otter Reef" | 
           QLD_NAME == "Satellite Reef") %>%
  add_row(QLD_NAME = "Pretty Patches Reef", X_COORD = 146.042, Y_COORD =-16.622)  %>% 
  add_row(QLD_NAME = "Russell Island", X_COORD = 146.094, Y_COORD =-17.229)

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
  coord_sf(xlim = c(142, 153), ylim = c(-25, -10), expand = F) +
  coord_sf(xlim = c(145, 150), ylim = c(-20, -15), expand = F) +
  theme(panel.background = element_rect(fill = "lightblue")) +
  xlab("")+ylab("")+ ggtitle("Core region reefs") +
  geom_point(data = myreefs.core, aes(x = X_COORD, y = Y_COORD), size = 5, 
             shape = 21, fill = "#2f3544") + 
  geom_label_repel(data = myreefs.core, aes(X_COORD, Y_COORD, label=QLD_NAME), 
                   fill = "white", 
                   nudge_x = 1.2, 
                   nudge_y = 0.5, 
                   max.iter = 5000)+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        panel.border = element_rect(colour = "#2f3544", fill=NA, size=2),
        axis.ticks.y=element_blank())
core.reefs

#---final figure ---# 
pdf("population_map_allReefs.pdf")
core.reefs
dev.off() 

ggsave("population_map_allReefs.jpeg", plot = core.reefs, width = 10, height = 7, units = "in", dpi = 300)

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
GBR_data <- st_read("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Statistics_wrkshp/Maps/TS_AIMS_NESP_Torres_Strait_Features_V1b_with_GBR_Features")
class(GBR_data)

#--- set working directory ---# 
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation")

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
  geom_sf(data = coords,fill = "gray90", color = "orange", size = 0.5) +
  geom_sf(data = AUS,fill = "gray90", color = "#333333", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(141, 153), ylim = c(-25, -10), expand = F) + 
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "lightblue"), 
        panel.grid.major = element_line(colour = "lightblue")) + 
  geom_point(data = myreefs.core, aes(x = X_COORD, y = Y_COORD), size = 3, 
             shape = 21, fill = "#DA3A36")+ 
  geom_point(data = myreefs.leading, aes(x = X_COORD, y = Y_COORD), size = 3, 
             shape = 21, fill = "#0D47A1") + 
  geom_rect(aes(xmin=145, xmax=147, ymin=-17.5, ymax=-16),fill="transparent", linetype = "dotted", color="#DA3A36") + 
  geom_rect(aes(xmin=148.5, xmax=151, ymin=-21.7, ymax=-20.2),fill="transparent", linetype = "dotted", color="#0D47A1") + 
  annotate("text", x = 149.8, y = -14, label = "Coral \nSea", fontface = "italic", size = 6);p1 

insert1 <- ggplot() +
  geom_sf(data = coords,fill = "gray90", color = "orange", size = 0.5) +
  geom_sf(data = AUS,fill = "gray90", color = "#333333", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  coord_sf(xlim = c(145, 147), ylim = c(-17.6, -16), expand = F) +
  theme(panel.background = element_rect(fill = "lightblue"), 
        panel.grid.major = element_line(colour = "lightblue")) +
  xlab("")+ylab("")+
  annotate("segment", x = 145.7926, xend = 145.7926+0.3, y = -16.32093, yend = -16.32093+0.18, colour = "black", size = 1)+
  annotate("segment", x = 146.2049, xend = 146.2049+0.3, y = -16.99838, yend = -16.99838+0.18, colour = "black", size = 1)+
  annotate("segment", x = 145.9929, xend = 145.9929+0.3, y = -16.65520, yend = -16.65520+0.18, colour = "black", size = 1)+
  geom_point(data = myreefs.core, aes(x = X_COORD, y = Y_COORD), size = 5, 
             shape = 21, fill = "#DA3A36") + 
  geom_label(data = myreefs.core, aes(X_COORD+.3, Y_COORD+.2, label=QLD_NAME), fill = "white")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        panel.border = element_rect(colour = "#DA3A36", fill=NA, size=2),
        axis.ticks.y=element_blank()); insert1

insert2 <- ggplot() +
  geom_sf(data = coords,fill = "gray90", color = "orange", size = 0.5) +
  geom_sf(data = AUS,fill = "gray90", color = "#333333", size = 0.1) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  coord_sf(xlim = c(148, 151), ylim = c(-21.7, -20.2), expand = F) +
  theme(panel.background = element_rect(fill = "lightblue"), 
        panel.grid.major = element_line(colour = "lightblue")) +  
  geom_point(data = myreefs.leading, aes(x = X_COORD, y = Y_COORD), size = 5, 
             shape = 21, fill = "#0D47A1")+ 
  xlab("") + ylab("") +
  geom_label_repel(data = myreefs.leading, aes(X_COORD, Y_COORD, label=GBR_NAME), fill = "white")+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        panel.border = element_rect(colour = "#0D47A1", fill=NA, size=2),
        axis.ticks.y=element_blank());insert2

#---final figure ---# 

cb <- ggarrange(insert1, insert2,
          ncol = 1, nrow = 2, labels = c("B","C"), 
          align = "hv",
          heights = c (1,1), 
          widths = c (1,1)) 

map <- ggarrange(p1, cb, 
                 ncol=2, 
                 nrow=1, 
                 labels = c("A","")); map

pdf("population_map.pdf")
map2
dev.off() 

ggsave("population_map2.jpeg", width = 11, height = 8, units = "in", dpi = 360)
ggsave("population_map2.pdf", width = 11, height = 8, units = "in", dpi = 360)

#--- 3MT final figure ---# 
mt_figure <- ggplot() +
  geom_sf(data = coords,fill = "NA", color = "white", size = 0.5) +
  geom_sf(data = AUS,fill = "NA", color = "white", size = 0.1) +
  #annotation_scale(location = "bl", width_hint = 0.5) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
                         #pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
                         #style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(141, 153), ylim = c(-25, -10), expand = F) + 
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "transparent", 
                                        color = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), 
        panel.grid.major = element_line(colour = "NA")) + 
  geom_point(data = myreefs.core, aes(x = X_COORD, y = Y_COORD), size = 3, 
             shape = 21, fill = "#DA3A36")+ 
  geom_point(data = myreefs.leading, aes(x = X_COORD, y = Y_COORD), size = 3, 
             shape = 21, fill = "#00FFFF") + 
  geom_rect(aes(xmin=145, xmax=147, ymin=-17.5, ymax=-16),fill="transparent", linetype = "dotted", color="#DA3A36") + 
  geom_rect(aes(xmin=148.5, xmax=151, ymin=-21.7, ymax=-20.2),fill="transparent", linetype = "dotted", color="#00FFFF") + 
  #annotate("text", x = 149.3, y = -17.5, label = "Coral \nSea", size = 7, fontface = "italic", color = "black", alpha = 0.3) + 
  annotate("text", x = 144.8, y = -21, label = "AUSTRALIA", size = 9, fontface = "italic", color = "black", alpha = 0.3) 
  #annotate("text", x = 146, y = -19.2589, label = "Townsville", fontface = "italic", size = 4) + 
  #annotate("text", x = 145.4, y = -16.92366, label = "Cairns", fontface = "italic", size = 4) + 
 # annotate("text", x = 148.3, y = -21.15345, label = "Mackay", fontface = "italic", size = 4);mt_figure 


ggsave('mt_figure.png',mt_figure,bg='transparent', 
       width = 8, height = 10, units = "in", dpi=300, device = 'png')


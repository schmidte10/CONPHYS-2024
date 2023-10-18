#--- load libraries ---# 
library(tidyverse)
library(ggpubr)

#--- recreate figure above but with legend ---# 
reefs <- rbind(CairnsTemp2, MackayTemp2) 
save(reefs, file="./sampled_reefs_temp_data.RData")
temperature.plot2 <- ggplot(reefs, aes(x=cal_val, fill=REGION)) + 
  geom_density(alpha = 0.8, 
               adjust=1.5) + 
  scale_fill_manual(values = c("#0D47A1","#DA3A36")) + 
  geom_vline(xintercept = c(27,28.5,30,31.5), 
             linetype = "dashed", 
             color ="black") +
  scale_x_continuous(breaks = seq(19.5,33,1.5))+
  theme_classic(); temperature.plot2 

months=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
temperature.plot3 <- ggplot(reefs, aes(x=MONTH, y=cal_val, color=REGION)) + 
  geom_smooth(size=2) +
  scale_x_continuous(labels=months, breaks=seq(1,12,1))+ 
  scale_color_manual(values = c("#0D47A1","#DA3A36"),  name = "Latitude") + 
  ggplot2::scale_y_continuous(limits=c(19.5,31.5),breaks = seq(19.5,33,1.5)) +
  geom_hline(yintercept = c(27,28.5,30,31.5), 
             linetype = "dashed", 
             color ="black")+
  xlab("Month")+ylab("Temperature (°C)") +
  theme_classic(); temperature.plot3

yplot <- ggdensity(reefs, "cal_val", fill="REGION") + 
  scale_fill_manual(values = c("#0D47A1","#DA3A36"),  name = "Latitude") + 
  ggplot2::scale_x_continuous(limits=c(19.5,31.5),breaks = seq(19.5,33,1.5)) +
  geom_vline(xintercept = c(27,28.5,30,31.5), 
             linetype = "dashed", 
             color ="black") +
  rotate() + clean_theme()
env.temp.plot <- ggpubr::ggarrange(temperature.plot3, yplot,
                                   ncol = 2, nrow=1, align = "hv", 
                                   widths = c(4, 1), heights = c(1, 2),
                                   common.legend = TRUE)

ggsave("env_temp_plot.pdf", width=20, height=12, units = "cm", dpi=360)
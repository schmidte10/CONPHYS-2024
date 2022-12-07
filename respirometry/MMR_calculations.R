#--- loading packages ---#
library(tidyverse)
library(respR)
library(janitor)

#--- import data ---#
Cycle_1 <- read_delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Resp/Dell/Experiment_ 26 June 2022 03 41PM/All slopes/Cycle_1.txt", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)
mass = 0.01844 
FISH_ID="INSERT_FISH_ID_HERE"
#--- format data ---# 
#make sure you are looking at the correct chamber
apoly <- Cycle_1 %>%
  dplyr::rename(time = `Seconds from start for linreg`, 
                oxygen = `ch2 po2`) %>%
  select(c("time","oxygen"))  

#inspect data
apoly.inspect <- inspect(apoly) 
#--- look at background respiration rates ---# 
# import background resp data file 
pre <- read_delim("./Resp/Dell/Experiment_ 26 June 2022 02 38PM/All slopes/Cycle_1.txt", 
                  delim = ";", escape_double = FALSE, trim_ws = TRUE)

pre <-  pre %>% 
  dplyr::rename(time = Time, 
         oxygen = `ch2 po2`) %>%
  select(c("time","oxygen"))  

bg_pre <- subset(pre, from = pre[1,1], to = tail(pre$time, n=1), by ="time") %>% 
  calc_rate.bg()
bg_pre

#--- calculate MMR ---# 
# wait time was 5 second added a 10 sec buffer to start measurements 
# 15 seconds after the fish was put  in to ensure proper mixing of water 
# water in chamber completed a circuit once every 15 seconds
#buffer =5
#measure = 240 

# get mmr 
apoly_mmr <- auto_rate(apoly, width = 0.25, method = "rolling")
apoly_mmr <- auto_rate(apoly, width = 60, by = "time", method = "rolling") 

# summary 
summary(apoly_mmr)
which.min(apoly_mmr$rate)

mymaxrate <- which.min(apoly_mmr$rate);

which.min(apoly_mmr$rate);plot(apoly_mmr, pos = mymaxrate)
min(apoly_mmr$rate)

# plot 
plot(apoly_mmr, pos = mymaxrate)
plot(apoly_mmr, pos = 1)

#--- adjusting MMR value for background resp ---# 
apoly_mmr_adj <- adjust_rate(apoly_mmr, 
                             by= bg_pre, 
                             by2= bg_pre,
                             method="linear")

summary(apoly_mmr_adj)

#--- convert MMR value to units that you will report values in ---# 

chamber.vol= 1.5
apoly_mmr_conv <- convert_rate(apoly_mmr_adj, 
                               oxy.unit = "%Air", 
                               time.unit = "secs", 
                               output.unit = "mg/h/kg", 
                               mass = mass,
                               volume = chamber.vol-mass, # volumr of chamber
                               S = 35, # salinity 
                               t = 27) # temperature 

summary(apoly_mmr_conv)
which.min(apoly_mmr_conv$rate.output) # which regression has the highest rate 
# note because oxygen declines the highest rate will be the lowest value
print(FISH_ID);min(apoly_mmr_conv$rate.output);min(apoly_mmr_conv$rate.output*mass);apoly_mmr_conv$summary$rsq[1] #lowest value/highest rate 



mmr_data <- summary(apoly_mmr_conv$summary$rsq, pos = 1, export = TRUE);mmr_data





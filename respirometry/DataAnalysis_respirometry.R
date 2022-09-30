#--- load libraries ---#  
library(tidyverse) 
library(plyr)
library(dplyr)

#--- set working directory ---#
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation")

#--- load data ---# 
resp <- read.delim("./respirometry/SummaryData_2022_resp.txt")

# data seems to have loaded with two extract columns at the end 
# remove extract columns by name 

resp <- resp %>% select(-c("X","X.1"))

#
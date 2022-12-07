#--- set working directory ---# 
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/hematocrit-hemoglobin") 

#--- load libraries ---# 
library(tidyverse)
library(readr)
library(ggplot2)


#--- load data ---# 
hema <- read_delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/import_files/HematocritHemoglobin.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE) 

#--- Data Manipulation ---# 
hema <-  hema %>% 
  mutate(PERC_RBC = as.numeric(PERC_RBC), 
         MASS = as.numeric(MASS),
         MASS_CENTERED = scale(MASS, center = TRUE, scale = FALSE))
#--- exploratory analysis ---# 
ggplot(hema, aes(y=PERC_RBC, x=REGION)) + 
  geom_boxplot() + 
  theme_classic() 

ggplot(hema, aes(y=HEMOGLOBIN, x= REGION)) + 
  geom_boxplot() + 
  theme_classic()

hist(hema$HEMOGLOBIN)
hist(hema$PERC_RBC)

hema %>% ggplot(aes(x=PERC_RBC)) + 
  geom_density() +
  facet_wrap(~REGION)

hema %>% ggplot(aes(x=HEMOGLOBIN)) + 
  geom_density() +
  facet_wrap(~REGION)

#--- model exploration ---# 

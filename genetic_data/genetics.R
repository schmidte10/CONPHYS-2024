#--- set working directory ---# 
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/genetic_data") 

#--- libraries ---# 
library(vcfR) 

#--- importing data ---# 
vcf_file <- read.vcfR("SCH11429_finalgatk.vcf") 

# this can take a couple minutes 
# go grab a coffee or something :) 
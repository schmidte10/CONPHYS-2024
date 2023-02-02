#--- load libraries ---#  
library(tidyverse) 
library(plyr)
library(dplyr)
library(lubridate) 
library(ggplot2)
library(lme4)
library(nlme) 
library(glmmTMB)
library(sjPlot)
library(gridExtra)
library(performance)
library(car)
library(DHARMa)
library(MuMIn)
library(kableExtra) 
library(broom)
library(emmeans)
library(ggeffects)
library(moments)
library(ggsignif)

#--- set working directory ---# 
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/resp-immunocompetence")

#--- load data ---# 
resp <- read.delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/SummaryData_2022_resp_updated.txt")

#--- import hematocrit data ---# 
hema <- read_delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/import_files/HematocritHemoglobin.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE) %>% 
  select(c(1:9))



#--- Data Manipulation - hematocrit ---# 
hema <-  hema %>% 
  mutate(PERC_RBC = as.numeric(PERC_RBC), 
         MASS = as.numeric(MASS),
         MASS_CENTERED = scale(MASS, center = TRUE, scale = FALSE)) %>% 
  drop_na(PERC_RBC) %>% 
  dplyr::rename(FISH_ID = FISHID)

#--- preparing data - respiration ---# 
resp2 = resp %>% 
  dplyr::rename(EXP_FISH_ID = FISH_ID) %>%
  separate(EXP_FISH_ID, c("FISH_ID"), remove = FALSE) %>%
  mutate(FISH_ID = factor(FISH_ID), 
         POPULATION = factor(POPULATION), 
         REGION = factor(REGION), 
         TEMPERATURE = factor(TEMPERATURE),
         RESTING_DATE = factor(RESTING_DATE), 
         RESTING_CHAMBER = factor(RESTING_CHAMBER), 
         RESTING_SYSTEM = factor(RESTING_SYSTEM), 
         RESTING_SUMP = factor(RESTING_SUMP), 
         RESTING_AM_PM = factor(RESTING_AM_PM), 
         RESTING_START_TIME = hms(RESTING_START_TIME),
         MAX_DATE = factor(MAX_DATE), 
         MAX_CHAMBER = factor(MAX_CHAMBER), 
         MAX_SYSTEM = factor(MAX_SYSTEM), 
         MAX_SUMP = factor(MAX_SUMP), 
         MAX_AM_PM = factor(MAX_AM_PM), 
         MAX_START_TIME = hms(MAX_START_TIME), 
         Swim.performance = factor(Swim.performance), 
         MgO2.hr_Net = as.numeric(MgO2.hr_Net)) %>% 
  dplyr::rename(MASS = DRY_WEIGHT) %>% 
  mutate(MASS_CENTERED = scale(MASS, scale = FALSE, center = TRUE)) %>% 
  drop_na(MASS)%>%
  drop_na(MgO2.hr_Net)

#--- remove individuals where min.max data is unreliable ---# 
resp3 <- resp2 %>% 
  subset(  
    EXP_FISH_ID !="LCHA127_27" & # deceased during experiment
      EXP_FISH_ID !="LCHA132_27" & # deceased during experiment
      EXP_FISH_ID !="LKES168_27" # poor data quality
    
  ) 

#--- remove individuals where ONLY max data is unreliable ---# 
resp4 <- resp3 %>% 
  subset(
    EXP_FISH_ID !="CSUD008_27" &  # poor swim
      EXP_FISH_ID !="CSUD008_30" &  # poor swim 
      EXP_FISH_ID !="CSUD008_31.5" & # poor swim
      EXP_FISH_ID !="CSUD018_31.5" & # poor swim 
      EXP_FISH_ID !="CSUD026_30" & # max. value low 
      EXP_FISH_ID !="CSUD074_28.5" & # fas value low 
      EXP_FISH_ID !="CSUD079_30" &
      EXP_FISH_ID !="CVLA052_27" & #nas value low 
      EXP_FISH_ID !="CVLA054_28.5" & # low max value? 
      EXP_FISH_ID !="LCHA113_27" & # poor data quality 
      EXP_FISH_ID !="LCHA113_30" & # poor swim 
      EXP_FISH_ID !="LCHA127_27" & # deceased during experiment
      EXP_FISH_ID !="LCHA114_28.5"  # poor swim 
  )    

#--- hematocrit was only collected at 31.5C so remove all other data points 
#--- from resp4 dataframe 

resp5 <- resp4 %>% 
  subset(TEMPERATURE == "31.5") 

#--- join data frames ---# 
hema.nas <- resp5 %>% 
  full_join(select(hema, c("FISH_ID","PERC_RBC")), by="FISH_ID") %>% 
  drop_na(EXP_FISH_ID) %>% 
  mutate(sPERC_RBC = as.numeric(scale(PERC_RBC)), 
         sMASS = as.numeric(scale(MASS)))

#--- model selection ---# 
model1.1 <- lm(MAX_MgO2.hr_RESPR ~ 1 + sPERC_RBC + REGION, 
               data = hema.nas)

model1.2 <- lm(MAX_MgO2.hr_RESPR ~ 1 + sPERC_RBC*REGION, 
               data = hema.nas)

model1.3 <- lm(MAX_MgO2.hr_RESPR ~ 1 + sPERC_RBC + REGION + sMASS, 
               data = hema.nas)

model1.4 <- lm(MAX_MgO2.hr_RESPR ~ 1 + sPERC_RBC * REGION + sMASS, 
               data = hema.nas)

AIC(model1.1, model1.2, model1.3, model1.4, k=2)
# model1.3 has the lowest AICc score 

#--- check model performance ---# 
model1.4 %>% check_model() # model distribution looks a bit off try a glm
pha.resid <-  model1.4 %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
model1.4 %>% DHARMa::testResiduals() 

#--- model selection --# 
model.gamma1 <- glm(MAX_MgO2.hr_RESPR ~ 1 + sPERC_RBC, 
                    data = hema.nas, 
                    family = Gamma(link = 'log'))


model.gamma2 <- glm(MAX_MgO2.hr_RESPR ~ 1 + sPERC_RBC + sMASS, 
               data = hema.nas, 
               family = Gamma(link="log")) 

model.gamma3 <- glm(MAX_MgO2.hr_RESPR ~ 1 + REGION + sPERC_RBC + sMASS, 
                    data = hema.nas, 
                    family = Gamma(link = 'log')) 

model.gamma4 <- glm(MAX_MgO2.hr_RESPR ~ 1 + sPERC_RBC * REGION + sMASS, 
                    data = hema.nas, 
                    family = Gamma(link = 'log')) 

AIC(model.gamma1, model.gamma2, model.gamma3, model.gamma4, k=2)
saveRDS(model.gamma3, "hema-nas_model_gamma_3.RDS")
#--- model performance ---# 
model.gamma4 %>% check_model()
pha.resid <-  model.gamma4 %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
model.gamma4 %>% DHARMa::testResiduals() 
# model looks good

#--- results ---#
model.gamma4 %>% plot_model(type="pred", terms=c('sPERC_RBC'), show.data=TRUE)
model.gamma4 %>% plot_model(type='eff',  terms=c('sPERC_RBC',"REGION"), show.data=TRUE)
figure1 <- model.gamma4 %>% plot_model(type='eff',  terms=c('sPERC_RBC',"REGION"), show.data=TRUE)
model.gamma4 %>% plot_model(type='est')

model.gamma4 %>% summary()
model.gamma4 %>% confint()

model.gamma4 %>% emtrends(var = "sPERC_RBC", type = "response")  %>% regrid() %>% summary(infer=TRUE)

pdf("hema_max.pdf", width = 7, height = 5)
print(figure1)
dev.off()

jpeg("hema_max.jpeg", units="in", width=7, height=5, res=300) 
print(figure1)
dev.off()

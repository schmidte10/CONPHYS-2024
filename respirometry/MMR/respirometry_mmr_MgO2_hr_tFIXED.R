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
#--- set working directory ---#
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/respirometry")

#--- load data ---# 
resp <- read.delim("./SummaryData_2022_resp.txt")
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/MMR/")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/respirometry/MMR")

#--- preparation of data ---# 
# data seems to have loaded with two extract columns at the end 
# remove extract columns by name 

resp <- resp %>% select(-c("X","X.1"))

#--- preparing data ---# 
resp2 = resp %>% 
  dplyr::rename(EXP_FISH_ID = FISH_ID) %>%
  separate(EXP_FISH_ID, c("FISH_ID"), remove = FALSE) %>%
  mutate(FISH_ID = factor(FISH_ID), 
         POPULATION = factor(POPULATION), 
         REGION = factor(REGION), 
         TEMPERATURE = factor(TEMPERATURE), #run with temperature as a factor
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
         NAS = as.numeric(NAS), 
         FAS = as.numeric(FAS), 
         MgO2.hr_Net = as.numeric(MgO2.hr_Net), 
         RESTING_RUNTIME_SECONDS = as.numeric(hms(RESTING_RUNTIME))) %>% 
  dplyr::rename(MASS = WEIGHT) %>% 
  mutate(MASS_CENTERED = scale(MASS, scale = FALSE, center = TRUE))

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
      EXP_FISH_ID !="LCHA127_27" # deceased during experiment
  )   


#--- exploratory data analysis ---# 
hist(resp4$MAX_MgO2.hr); shapiro.test(resp4$MAX_MgO2.hr) 

resp4 %>% 
  group_by(REGION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(MAX_MgO2.hr), 
                   Max. = max(MAX_MgO2.hr), 
                   Mean = mean(MAX_MgO2.hr)) 

#--- model formula ---# 
#--- base model ---# 
mmr.1 <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4,
                 REML = FALSE) 


#--- experimental mmr equipment hypothesis ---#
mmr.2 <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MAX_SUMP + MAX_CHAMBER + 
                   MAX_AM_PM, 
                 family=gaussian(),
                 data = resp4,
                 REML = FALSE)  

mmr.3 <-  glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + steepest_slope, 
              family=gaussian(),
              data = resp4,
              REML = FALSE)  

AIC(mmr.1, mmr.2, mmr.3, k=2)
#followed by random effects
mmr.1 <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4,
                 REML = TRUE) 

mmr.1a <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE) 

mmr.1b <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|POPULATION/FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

mmr.1c <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID) + (REGION|POPULATION), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

AIC(mmr.1, mmr.1a, mmr.1b, mmr.1c,  k=2)

#--- Final model ---# 
mmr.1b <- glmmTMB(MAX_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|POPULATION/FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

#--- saving model ---#
saveRDS(mmr.1b, file = "mmr_1b.RDS") 

#--- load model ---# 
mmr.1b <- readRDS("mmr_1b.RDS")

#--- investigate model ---#
#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")
check_model(mmr.1b)

mmr.1b %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
mmr.1b %>% plot_model(type='est')

mmr.1b %>% summary()
mmr.1b %>% confint()
mmr.1b  %>% performance::r2_nakagawa()

mmr.1b %>% emmeans(~ TEMPERATURE*REGION)
mmr.1b %>% emmeans(~ TEMPERATURE*REGION) %>% pairs(by = "TEMPERATURE") %>% summary(infer=TRUE)
mmr.1b %>% emmeans(~ TEMPERATURE*REGION) %>% pairs(by = "REGION") %>% summary(infer=TRUE)

#--- plot ---#
mmr.newdata <- mmr.1b %>% ggemmeans(~TEMPERATURE|REGION) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

mmr.g1 <- ggplot(mmr.newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); mmr.g1

predict(mmr.1b, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(mmr.1b, type='response') 
#data points based on month/situation/random effects - to get the data points
obs <-  resp4 %>% 
  mutate(Pred=predict(mmr.1b, re.form=NA),
         Resid = residuals(mmr.1b, type='response'),
         Fit = Pred + Resid)
obs %>% head() 
mmr.g2 <- ggplot(mmr.newdata, aes(y=predicted, x=TEMPERATURE, color=group))+
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape=19,
                  size=1,
                  position=position_dodge(0.2)) + 
  #scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+ 
  scale_y_continuous(limits = c(11,18), breaks = seq(11, 18, by = 2)) +
  theme_classic() + ylab("MAXIMUM METABOLIC RATE (MMR: MgO2/hr)")+
  scale_color_manual(values=c("#DA3A36", "#0D47A1")) + 
  theme(legend.position = 'none')+ 
  geom_signif( 
    y_position = c(17.5,17.7), 
    xmin = c(2.75,3.75), 
    xmax = c(3.15,4.15), 
    annotations = c("*","*"), 
    tip_length = 0, 
    color = "black"); mmr.g2
#+ geom_signif(
  #y_position = c(13.91+0.5, 15.10+0.5,15.80+0.5,16.06+0.5), xmin = c(0.8, 1.8,2.8,3.8), xmax = c(1.2,2.2,3.2,4.2),
  #annotation = c("ns", "ns", "**\np =0.020", "**\np =0.010"), tip_length = 0.025, color = "black"); g2

pdf("mmr_1b.pdf", width= 7, height = 5)
print(g2)
dev.off()

jpeg("mmr_1b.jpeg", units="in", width=7, height=5, res=300) 
print(g2)
dev.off()
#########################################################################################























resp4 %>% 
  group_by(POPULATION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(NAS), 
                   Max. = max(NAS), 
                   Mean = mean(NAS)) %>% 
  print(n = 24)

#--- model formula ---# 
#POPULATION - MAX metablic rate
pop.mmr <- glmmTMB(MAX_MgO2.hr ~ 1+ POPULATION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                   family=gaussian(),
                   data = resp4,
                   REML = FALSE)

#--- save model ---# 
saveRDS(pop.mmr, file = "glmmTMB_mmr_p3_population_MgO2.hr.RDS") 

#--- load model ---#
#pop.mmr.p3_MgO2.hr <- readRDS("glmmTMB_mmr_p3_population_MgO2.hr.RDS")

#--- investigate model ---#
pop.mmr.p2 %>% plot_model(type='eff',  terms=c('TEMPERATURE','POPULATION'), show.data=TRUE)
pop.mmr.p2 %>% ggemmeans(~TEMPERATURE|POPULATION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
pop.mmr.p2 %>% plot_model(type='est')

pop.mmr.p2 %>% summary()
pop.mmr.p2 %>% confint()
pop.mmr.p2  %>% r.squaredGLMM()

pop.mmr.p2 %>% emmeans("POPULATION", var="TEMPERATURE") %>% pairs() %>% summary(infer=TRUE)

##########################################################################################


#--- custome contrast ---# 

emm1 = emmeans(pop.mmr, specs = ~ POPULATION*TEMPERATURE); emm1

#--- temperature treatment - 27C - ---#
core27 = c(0,0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
mackay27 = c(0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
chauvel27 = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 

contrast27 = emmeans(pop.mmr, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core27 - mackay27, 
                                                                                           "core - chauvel" = core27 - chauvel27, 
                                                                                           "mackay - chauvel" = mackay27 -chauvel27)) %>% summary(infer=TRUE) 

#--- temperature treatment - 28.5C - ---#
core28.5 =    c(0,0,0,0,0,0,0,0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0,0,0) 
mackay28.5 =  c(0,0,0,0,0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
chauvel28.5 = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 

contrast28.5 = emmeans(pop.mmr, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core28.5 - mackay28.5, 
                                                                                             "core - chauvel" = core28.5 - chauvel28.5, 
                                                                                             "mackay - chauvel" = mackay28.5 -chauvel28.5)) %>% summary(infer=TRUE) 

#--- temperature treatment - 30.0C - ---#
core30 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/3,1/3,1/3,0,0,0,0,0,0) 
mackay30 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0) 
chauvel30 = c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0) 

contrast30 = emmeans(pop.mmr, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core30 - mackay30, 
                                                                                           "core - chauvel" = core30 - chauvel30, 
                                                                                           "mackay - chauvel" = mackay30 -chauvel30)) %>% summary(infer=TRUE) 

#--- temperature treatment - 31.5C - ---#
core31.5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/3,1/3,1/3) 
mackay31.5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/2,1/2,0,0,0) 
chauvel31.5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0) 

contrast31.5 = emmeans(pop.mmr, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core31.5 - mackay31.5, 
                                                                                             "core - chauvel" = core31.5 - chauvel31.5, 
                                                                                             "mackay - chauvel" = mackay31.5 - chauvel31.5)) %>% summary(infer=TRUE)


#--- contrasts ---# 
print("TEMPERATURE - 27"); contrast27; print("TEMPERATURE - 28.5"); contrast28.5;print("TEMPERATURE - 30"); contrast30; print("TEMPERATURE - 31.5"); contrast31.5



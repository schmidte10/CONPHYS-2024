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
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry")

#--- load data ---# 
resp <- read.delim("./SummaryData_2022_resp.txt")
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/MMR/")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/MMR")

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
         Swim.performance = factor(Swim.performance)) %>% 
  dplyr::rename(MASS = WEIGHT) %>% 
  mutate(MASS_CENTERED = scale(MASS, scale = FALSE, center = TRUE))

#--- remove individuals where min.max data is unreliable ---# 
resp3 <- resp2 %>% 
  subset(EXP_FISH_ID !="LCHA132_27" & 
           EXP_FISH_ID != "LKES168_27" & 
           EXP_FISH_ID != "LCHA113_30" & 
           EXP_FISH_ID != "CSUD088_27" & 
           EXP_FISH_ID != "CTON062_27") 

#--- remove individuals where ONLY max data is unreliable ---# 
resp4 <- resp3 %>% 
  subset(EXP_FISH_ID !="CSUD014_27" & 
           EXP_FISH_ID != "CTON065_27" & 
           EXP_FISH_ID != "CTON069_30" & 
           EXP_FISH_ID != "CVLA054_27" & 
           EXP_FISH_ID != "LCHA129_27"& 
           EXP_FISH_ID != "LCHA135_27" & 
           EXP_FISH_ID != "LCKM162_27" & 
           EXP_FISH_ID != "LCKM165_27" & 
           EXP_FISH_ID != "LCKM180_27"& 
           EXP_FISH_ID != "CSUD026_30" & 
           EXP_FISH_ID != "CTON067_28.5" & 
           EXP_FISH_ID != "CVLA054_28.5" & 
           EXP_FISH_ID != "LCHA114_28.5"& 
           EXP_FISH_ID != "LCKM163_28.5" & 
           EXP_FISH_ID != "LCKM154_31.5"& 
           EXP_FISH_ID != "CSUD079_30" & 
           EXP_FISH_ID != "CSUD079_31.5"& 
           EXP_FISH_ID != "LCHA125_30")

#--- exploratory data analysis ---# 
hist(resp4$MAX_MgO2.hr); shapiro.test(resp4$MAX_MgO2.hr) 

resp4 %>% 
  group_by(REGION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(MAX_MgO2.hr), 
                   Max. = max(MAX_MgO2.hr), 
                   Mean = mean(MAX_MgO2.hr)) 

#--- model formula ---# 
#max metablic rate
mmr_MgO2.hr_tfixed <- glmmTMB(MAX_MgO2.hr ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + MAX_CHAMBER + MAX_SUMP + (1|REGION:POPULATION) + (1|FISH_ID), 
               family=gaussian(),
               data = resp4,
               REML = TRUE)


#--- saving model ---#
saveRDS(mmr_MgO2.hr_tfixed, file = "glmmTMB_mmr_MgO2_hr_tfixed.RDS") 

#--- load model ---# 
#mmr.p3_MgO2.hr <- readRDS("glmmTMB_mmr_p3_MgO2.hr.RDS")

#--- investigate model ---#
#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")
check_model(mmr_MgO2.hr_tfixed)

mmr_MgO2.hr_tfixed %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
mmr_MgO2.hr_tfixed %>% plot_model(type='est')

mmr_MgO2.hr_tfixed %>% summary()
mmr_MgO2.hr_tfixed %>% confint()
mmr_MgO2.hr_tfixed  %>% r.squaredGLMM()
mmr_MgO2.hr_tfixed  %>% performance::r2_nakagawa()

mmr_MgO2.hr_tfixed %>% emmeans(~ TEMPERATURE*REGION) %>% pairs(by = "TEMPERATURE") %>% summary(infer=TRUE)

#--- plot ---#
newdata <- mmr_MgO2.hr_tfixed %>% ggemmeans(~TEMPERATURE|REGION) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); g1

predict(mmr_MgO2.hr_tfixed, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(mmr_MgO2.hr_tfixed, type='response') 
#data points based on month/situation/random effects - to get the data points
obs <-  resp4 %>% 
  mutate(Pred=predict(mmr_MgO2.hr_tfixed, re.form=NA),
         Resid = residuals(mmr_MgO2.hr_tfixed, type='response'),
         Fit = Pred + Resid)
obs %>% head() 
g2 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group))+
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape=19,
                  size=1,
                  position=position_dodge(0.2)) + 
  #scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+ 
  scale_y_continuous(limits = c(11,18), breaks = seq(11, 18, by = 2)) +
  theme_classic() + ylab("MAXIMUM METABOLIC RATE (MMR: MgO2/hr)")+
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), name = "Regions")+ geom_signif(
  y_position = c(13.91+0.5, 15.10+0.5,15.80+0.5,16.06+0.5), xmin = c(0.8, 1.8,2.8,3.8), xmax = c(1.2,2.2,3.2,4.2),
  annotation = c("ns", "ns", "**\np =0.020", "**\np =0.010"), tip_length = 0.025, color = "black"); g2

pdf("MMR_MgO2_hr_tFIXED2.pdf", width= 7, height = 5)
print(g2)
dev.off()

jpeg("MMR_MgO2_hr_tFIXED2.jpeg", units="in", width=7, height=5, res=300) 
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
pop.mmr_MgO2.hr <- glmmTMB(MAX_MgO2.hr ~ 1+ POPULATION * TEMPERATURE + MASS_CENTERED + MAX_CHAMBER + MAX_SUMP + (1|FISH_ID), 
                   family=gaussian(),
                   data = resp3,
                   REML = TRUE)



check_model(pop.mmr.p3_MgO2.hr) 

#--- save model ---# 
saveRDS(pop.mmr.p3_MgO2.hr, file = "glmmTMB_mmr_p3_population_MgO2.hr.RDS") 

#--- load model ---#
#pop.mmr.p3_MgO2.hr <- readRDS("glmmTMB_mmr_p3_population_MgO2.hr.RDS")

#--- investigate model ---#
pop.mmr.p3_MgO2.hr %>% plot_model(type='eff',  terms=c('TEMPERATURE','POPULATION'), show.data=TRUE)
pop.mmr.p3_MgO2.hr %>% ggemmeans(~TEMPERATURE|POPULATION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
pop.mmr.p3_MgO2.hr %>% plot_model(type='est')

pop.mmr.p3_MgO2.hr %>% summary()
pop.mmr.p3_MgO2.hr %>% confint()
pop.mmr.p3_MgO2.hr  %>% r.squaredGLMM()

pop.mmr.p3_MgO2.hr %>% emtrends("POPULATION", var="TEMPERATURE") %>% pairs() %>% summary(infer=TRUE)

##########################################################################################





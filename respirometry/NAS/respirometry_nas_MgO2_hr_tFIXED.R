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
# personal computer
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry")
#--- load data ---# 
resp <- read.delim("./SummaryData_2022_resp_updated.txt")
# personal computer
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/NAS/")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/NAS")


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
         RESTING_END_TIME = hms(RESTING_ENDTIME),
         MAX_DATE = factor(MAX_DATE), 
         MAX_CHAMBER = factor(MAX_CHAMBER), 
         MAX_SYSTEM = factor(MAX_SYSTEM), 
         MAX_SUMP = factor(MAX_SUMP), 
         MAX_AM_PM = factor(MAX_AM_PM), 
         MAX_START_TIME = hms(MAX_START_TIME), 
         Swim.performance = factor(Swim.performance), 
         MgO2.hr_Net = as.numeric(MgO2.hr_Net), 
         RESTING_RUNTIME_SECONDS = as.numeric(hms(RESTING_RUNTIME))) %>% 
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


#--- exploratory data analysis ---# 
hist(resp4$MgO2.hr_Net); shapiro.test(resp4$MgO2.hr_Net); skewness(resp4$MgO2.hr_Net) #postive/right skewed
resp4 %>% ggplot(aes(x=TEMPERATURE, y=MgO2.hr_Net, fill = REGION)) + geom_boxplot() 
resp4 %>% ggplot(aes(x=MgO2.hr_Net, fill=TEMPERATURE)) + geom_density(alpha = 0.5) + 
  facet_wrap(~TEMPERATURE)

hist(sqrt(resp4$MgO2.hr_Net)); shapiro.test(sqrt(resp4$MgO2.hr_Net)); skewness(sqrt(resp4$MgO2.hr_Net)) # sqrt transforming data seems to fix skewness
resp4 %>% ggplot(aes(x=TEMPERATURE, y=sqrt(MgO2.hr_Net), fill = REGION)) + geom_boxplot() 
resp4 %>% ggplot(aes(x=sqrt(MgO2.hr_Net), fill=TEMPERATURE)) + geom_density(alpha = 0.5) + 
  facet_wrap(~TEMPERATURE)
#--- make sqrt transformed column ---# 
resp4 <- resp4 %>% 
  mutate(sqrt.MgO2.hr_NET = sqrt(MgO2.hr_Net))

resp4 %>% 
  group_by(REGION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(MgO2.hr_Net), 
                   Max. = max(MgO2.hr_Net), 
                   Mean = mean(MgO2.hr_Net))

pop.sample.size <- resp4 %>% 
  group_by(POPULATION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(MgO2.hr_Net), 
                   Max. = max(MgO2.hr_Net), 
                   Mean = mean(MgO2.hr_Net))

#--- model formula ---# 
#--- base model ---# 
nas.1 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
                              family=gaussian(),
                              data = resp4,
                              REML = FALSE) 
 
#--- experimental resting equipment hypothesis ---#
nas.2 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_SUMP + RESTING_RUNTIME_SECONDS + 
                                RESTING_AM_PM, 
                              family=gaussian(),
                              data = resp4,
                              REML = FALSE) 

#--- experimental max equipment hypothesis ---#
nas.3 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MAX_SUMP + MAX_CHAMBER + 
                   MAX_AM_PM, 
                 family=gaussian(),
                 data = resp4,
                 REML = FALSE) 

AICc(nas.1, nas.2, nas.3, k=2)
#--- followed by how to treat temperature
nas.1 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4,
                 REML = TRUE) 

nas.1q <-  glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED, 
                   family=gaussian(),
                   data = resp4,
                   REML = FALSE)  

nas.1p <-  glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE, 3) + MASS_CENTERED, 
                   family=gaussian(),
                   data = resp4,
                   REML = FALSE) 

AICc(nas.1, nas.1q, nas.1p, k=2) 

#--- followed by inclusion of random variables
nas.1a <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                 family=gaussian(),
                 data = resp4,
                 REML = TRUE) 

nas.1b <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID) + (1|POPULATION), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

nas.1c <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID) + (1|REGION/POPULATION), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

nas.1d <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID) + (REGION|POPULATION), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

AICc(nas.1, nas.1a, nas.1b, nas.1c, nas.1d, k=2)

#--- Final model ---# 
nas.1a <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

#--- saving model ---#
saveRDS(nas.1a, file = "nas_1a.RDS") 

#--- load model ---# 
# nas.1a <- readRDS("nas_1a.RDS") 

#--- investigate model ---#
#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")
check_model(nas.1a)
pha.resid <-  nas.1a %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
nas.1a %>% DHARMa::testResiduals() 

#sim <- simulateResiduals(nas.1a)
#which(residuals(sim) == 1 | residuals(sim) == 0)

nas.1a %>% plot_model(type='eff',  terms=c('TEMPERATURE','REGION'), show.data=TRUE)
nas.1a %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
nas.1a %>% plot_model(type='est')

nas.1a %>% summary()
nas.1a %>% confint()
nas.1a  %>% r.squaredGLMM()
nas.1a  %>% performance::r2_nakagawa()

nas.1a %>% emmeans(~ TEMPERATURE*REGION, type = "response")  %>% summary(infer=TRUE)
nas.1a %>% emmeans(~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(infer=TRUE)
nas.1a %>% emmeans(~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "REGION") %>% summary(infer=TRUE)
#--- plot ---#
newdata <- nas.1a %>% ggemmeans(~TEMPERATURE|REGION) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); g1

predict(nas.1a, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(nas.1a, type='response')
#data points based on month/situation/random effects - to get the data points
obs <-  resp4 %>% 
  mutate(Pred=predict(nas.1a, re.form=NA),
         Resid = residuals(nas.1a, type='response'),
         Fit = Pred + Resid)
obs %>% head() 
g2 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) + 
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape=19,
                  size=1,
                  position=position_dodge(0.2)) + 
  scale_y_continuous(limits = c(6,13), breaks = seq(6, 13, by = 2)) + 
  #scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+
  theme_classic() + ylab("NET AEROBIC SCOPE (NAS: MgO2/hr)") +
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Cairns (north)","Mackay (south)"),
                     name = "Regions")#+ 
  #geom_signif(
    #y_position = c(4.11+1.5, 5.18+1.5,5.15+1.5,4.66+1.5), xmin = c(0.8, 1.8,2.8,3.8), xmax = c(1.2,2.2,3.2,4.2),
   # annotation = c("ns", "ns", "**\np =0.046", "ns"), tip_length = 0.025, color = "black"); g2

pdf("nas_1a.pdf", width = 7, height = 5)
print(g2)
dev.off()

jpeg("nas_1a.jpeg", units="in", width=7, height=5, res=300) 
print(g2)
dev.off()
##########################################################################################


##########################################################################################
##########################################################################################
####################                                                    ##################
####################                 population analysis                ##################
###################                                                     ##################
##########################################################################################


# thus indicating that there is little difference in variation explained by populations

#--- custome contrast ---# 
nas.pop <- glmmTMB(MgO2.hr_Net ~ 1+ TEMPERATURE*POPULATION + MASS_CENTERED + (1|FISH_ID), 
                   family=gaussian(),
                   data = resp4,
                   REML = FALSE) 
saveRDS(nas.pop, "population.RDS")


emm1 = emmeans(nas.pop, specs = ~ POPULATION*TEMPERATURE); emm1

#--- temperature treatment - 27C - ---#
core27 = c(0,0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
mackay27 = c(0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
chauvel27 = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 

contrast27 = emmeans(nas.pop, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core27 - mackay27, 
                                                                                     "core - chauvel" = core27 - chauvel27, 
                                                                                     "mackay - chauvel" = mackay27 -chauvel27)) %>% summary(infer=TRUE) 
contrast.27 <- as.data.frame(contrast27)
#--- temperature treatment - 28.5C - ---#
core28.5 =    c(0,0,0,0,0,0,0,0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0,0,0) 
mackay28.5 =  c(0,0,0,0,0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
chauvel28.5 = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 

contrast28.5 = emmeans(nas.pop, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core28.5 - mackay28.5, 
                                                                                           "core - chauvel" = core28.5 - chauvel28.5, 
                                                                                           "mackay - chauvel" = mackay28.5 -chauvel28.5)) %>% summary(infer=TRUE) 
contrast.28.5 <- as.data.frame(contrast28.5)
#--- temperature treatment - 30.0C - ---#
core30 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/3,1/3,1/3,0,0,0,0,0,0) 
mackay30 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0) 
chauvel30 = c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0) 

contrast30 = emmeans(nas.pop, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core30 - mackay30, 
                                                                                           "core - chauvel" = core30 - chauvel30, 
                                                                                           "mackay - chauvel" = mackay30 -chauvel30)) %>% summary(infer=TRUE) 
contrast.30 <- as.data.frame(contrast30)
#--- temperature treatment - 31.5C - ---#
core31.5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/3,1/3,1/3) 
mackay31.5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/2,1/2,0,0,0) 
chauvel31.5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0) 

contrast31.5 = emmeans(nas.pop, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core31.5 - mackay31.5, 
                                                                                           "core - chauvel" = core31.5 - chauvel31.5, 
                                                                                           "mackay - chauvel" = mackay31.5 - chauvel31.5)) %>% summary(infer=TRUE)

contrast.31.5 <- as.data.frame(contrast31.5)
#--- contrasts ---# 
print("TEMPERATURE - 27"); contrast27; print("TEMPERATURE - 28.5"); contrast28.5;print("TEMPERATURE - 30"); contrast30; print("TEMPERATURE - 31.5"); contrast31.5

#--- figure ---# 
resp7 <- resp4 |> 
  mutate(LOCATION = case_when(POPULATION == "Sudbury Reef" ~ "Core", 
                              POPULATION == "Vlassof Cay" ~ "Core",
                              POPULATION == "Tongue Reef" ~ "Core",
                              POPULATION == "Cockermouth Island" ~ "Mackay (inshore)",
                              POPULATION == "Keswick Island" ~ "Mackay (inshore)",
                              POPULATION == "Chauvel Reef" ~ "Chauvel")) 

#--- plot ---# 
cc.df <- rbind(contrast.27, contrast.28.5, contrast.30, contrast.31.5)
nas.pop2 <- glmmTMB(MgO2.hr_Net ~ 1+ TEMPERATURE*LOCATION + MASS_CENTERED + (1|FISH_ID), 
                   family=gaussian(),
                   data = resp7,
                   REML = FALSE)
nas.pop2 %>% plot_model(type='eff',  terms=c('TEMPERATURE','LOCATION'), show.data=TRUE)

newdata <- nas.pop2 %>% ggemmeans(~TEMPERATURE|LOCATION) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x) 


g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); g1

predict(nas.pop2, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(nas.pop2, type='response')
#data points based on month/situation/random effects - to get the data points
obs <-  resp4 %>% 
  mutate(Pred=predict(nas.pop2, re.form=NA),
         Resid = residuals(nas.pop2, type='response'),
         Fit = Pred + Resid)
obs %>% head() 
g2 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) + 
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape=19,
                  size=1,
                  position=position_dodge(0.2)) + 
  scale_y_continuous(limits = c(4,14), breaks = seq(4, 14, by = 2)) + 
  #scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+
  theme_classic() + ylab("NET AEROBIC SCOPE (NAS: MgO2/hr)") +
  scale_color_manual(values=c("#DA3A36", "orange", "#0D47A1"), labels = c("Cairns", "Chauvel Reef (southern)","Mackay (inshore)"),
                     name = "Regions"); g2
#+ 
#geom_signif(
#y_position = c(4.11+1.5, 5.18+1.5,5.15+1.5,4.66+1.5), xmin = c(0.8, 1.8,2.8,3.8), xmax = c(1.2,2.2,3.2,4.2),
# annotation = c("ns", "ns", "**\np =0.046", "ns"), tip_length = 0.025, color = "black"); g2

pdf("population_nas.pdf", width = 7, height = 5)
print(g2)
dev.off()

jpeg("population_nas.jpeg", units="in", width=7, height=5, res=300) 
print(g2)
dev.off()

# may be worth investigating individual performance curves. 
######################################################################################################################################


##########################################################################################
##########################################################################################
####################                                                    ##################
####################       re-running model without chauvel reef        ##################
###################                                                     ##################
##########################################################################################
resp5 <- resp4 %>% 
  subset(POPULATION != "Chauvel Reef")


nas.1a.woc <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                      family=gaussian(),
                      data = resp5,
                      REML = TRUE) 

#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")
check_model(nas.1a.woc)
pha.resid <-  nas.1a.woc %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
nas.1a.woc %>% DHARMa::testResiduals() 

#sim <- simulateResiduals(nas.1a.woc)
#which(residuals(sim) == 1 | residuals(sim) == 0)

nas.1a.woc %>% plot_model(type='eff',  terms=c('TEMPERATURE','REGION'), show.data=TRUE)
nas.1a.woc %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
nas.1a.woc %>% plot_model(type='est')

nas.1a.woc %>% summary()
nas.1a.woc %>% confint()
nas.1a.woc  %>% r.squaredGLMM()
nas.1a.woc  %>% performance::r2_nakagawa()

nas.1a.woc %>% emmeans(~ TEMPERATURE*REGION, type = "response")  %>% summary(infer=TRUE)
nas.1a.woc %>% emmeans(~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(infer=TRUE)

# CONCLUSION: Significant differences remain unchanged, however results are slightly 'less significant'

#################################################################################################

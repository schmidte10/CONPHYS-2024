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
library(vtable)
#--- set working directory ---#
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry")
#--- load data ---# 
resp <- read.delim("./SummaryData_2022_resp_updated.txt")
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/RMR/")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/RMR")

# data seems to have loaded with two extract columns at the end 

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

#--- remove individuals where data is irregular ---# 
resp3 <- resp2 %>% 
  subset(  
    EXP_FISH_ID !="LCHA127_27" & # deceased during experiment
      EXP_FISH_ID !="LCHA132_27" & # deceased during experiment
      EXP_FISH_ID !="LKES168_27" # poor data quality
    
  ) 
###--- EXPLORATORY ANALYSIS ----####
table(resp3$REGION, resp3$RESTING_CHAMBER, resp3$TEMPERATURE)
#--- exploratory data analysis: covariates ---# 
ggplot(resp3, aes(MASS, RESTING_MgO2.hr_RESPR)) + 
  geom_point() + 
  geom_smooth(method = "lm")
# negative relationships between mass and resting metabolic rate  

ggplot(resp3, aes(MASS, MAX)) + 
  geom_point() + 
  geom_smooth(method = "lm")
# positive relationship between mass and maximum metabolic rate

ggplot(resp3, aes(MASS, NAS)) + 
  geom_point() + 
  geom_smooth(method = "lm")
# No relationship between mass and net aerobic cope

ggplot(resp3, aes(MASS, FAS)) + 
  geom_point() + 
  geom_smooth(method = "lm")
# slight positive relationship between mass and factorial aerobic scope

# mass should be taken into account and used as a covariate in models
# when looking at resting and maximum metabolic rate, and factorial aerobic scope

#--- exploratory data analysis: populations ---# 
ggplot(resp3, aes(REGION, RESTING_MgO2.hr_RESPR)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE)

ggplot(resp3, aes(REGION, MAX)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE)

ggplot(resp3, aes(MASS, RESTING_MgO2.hr_RESPR, color = REGION)) + 
  geom_point() +
  theme_classic() + 
  geom_smooth(method = "lm")

ggplot(resp3, aes(REGION, RESTING_MgO2.hr_RESPR)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE)

ggplot(resp3, aes(REGION, MAX_MgO2.hr)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE) 

ggplot(resp3, aes(REGION, MgO2.hr_Net)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE) 

ggplot(resp3, aes(REGION, FAS)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE) 
#################################################################################################################
#--- exploratory data analysis ---# 
hist(resp3$RESTING_MgO2.hr); shapiro.test(resp3$RESTING_MgO2.hr)

resp3 %>% 
  group_by(REGION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(RESTING_MgO2.hr_RESPR), 
                   Max. = max(RESTING_MgO2.hr_RESPR), 
                   Mean = mean(RESTING_MgO2.hr_RESPR)) 

#--- model formula ---# 
#--- base model ---#
rmr.1 <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp3,
                 REML = FALSE) 

#--- experimental rmr equipment hypothesis ---#
rmr.2 <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + RESTING_SUMP + 
                   RESTING_AM_PM + RESTING_RUNTIME_SECONDS, 
                 family=gaussian(),
                 data = resp3,
                 REML = FALSE) 
summary(rmr.2)
#--- base model and resting runtime ---#
rmr.3 <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS, 
                 family=gaussian(),
                 data = resp3,
                 REML = FALSE)

AICc(rmr.1, rmr.2, rmr.3, k=2)
#followed by random effects
rmr.3 <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS, 
                 family=gaussian(),
                 data = resp3,
                 REML = TRUE) 

rmr.3a <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp3,
                  REML = TRUE) 

rmr.3b <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID) + (1|POPULATION), 
                  family=gaussian(),
                  data = resp3,
                  REML = TRUE)

rmr.3c <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID) + (1|REGION/POPULATION), 
                  family=gaussian(),
                  data = resp3,
                  REML = TRUE)

rmr.3d <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID) + (REGION|POPULATION), 
                  family=gaussian(),
                  data = resp3,
                  REML = TRUE)

AICc(rmr.3, rmr.3a, rmr.3b, rmr.3c, rmr.3d, k=2)

#--- Final model ---# 
rmr.3a <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp3,
                  REML = TRUE)

#--- saving model ---#
saveRDS(rmr.3a, file = "rmr_3a.RDS") 


#--- load model ---#
rmr.3a <- readRDS("rmr_3a.RDS")

#--- investigate model ---#

rmr.3a %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
rmr.3a %>% plot_model(type='est')

rmr.3a %>% summary()
rmr.3a %>% confint()
rmr.3a  %>% performance::r2_nakagawa()

rmr.3a %>% emmeans(~ TEMPERATURE*REGION)
rmr.3a %>% emmeans(~ TEMPERATURE*REGION) %>% pairs(by = "TEMPERATURE") %>% summary(infer=TRUE)

#--- plot ---#
newdata <- rmr.3a %>% 
  ggemmeans(~TEMPERATURE|REGION) %>% 
  as.data.frame() %>% 
  dplyr::rename(TEMPERATURE = x)

g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); g1

predict(rmr.3a, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(rmr.3a, type='response') 
#data points based on month/situation/random effects - to get the data points
obs <-  resp3 %>% 
  mutate(Pred=predict(rmr.3a, re.form=NA),
         Resid = residuals(rmr.3a, type='response'),
         Fit = Pred + Resid)
obs %>% head() 
g2 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group, fill = group)) + 
  geom_line(method="lm", stat="smooth",
            formula=y ~ poly(x, 3, raw=FALSE), 
            aes(color = group), 
            size = 1, 
            alpha = 0.5)+
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape=21,
                  size=1,
                  position=position_dodge(0.2)) + 
  #scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+
  theme_classic() + ylab("RESTING_MgO2.hr (RMR: MgO2/hr)")+
  scale_fill_manual(values=c("#DA3A36", "#0D47A1"), 
                    labels = c("Cairns (north)","Mackay (south)"), 
                    name = "Regions")+ 
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), 
                     labels = c("Cairns (north)","Mackay (south)"), 
                     name = "Regions"); g2

pdf("rmr_3a.pdf")
print(g2)
dev.off()

jpeg("rmr_3a.jpeg", units="in", width=7, height=5, res=300)
print(g2)
dev.off()

##########################################################################################
################# POPULATION POPULATION POPULATION ####################################### 
##########################################################################################

resp3 %>% 
  group_by(POPULATION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(RESTING_MgO2.hr), 
                   Max. = max(RESTING_MgO2.hr), 
                   Mean = mean(RESTING_MgO2.hr)) %>% 
  print(n = 24)

#--- model formula ---#############################################################################################
# Script below has not been run/edited - in current state - DO NOT USE
#POPULATION - resting metablic rate
pop.rest_MgO2.hr <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ POPULATION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID), 
                    family=gaussian(),
                    data = resp3,
                    REML = TRUE)


pop.rest.poly2_MgO2.hr <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ POPULATION * poly(TEMPERATURE, 2) + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID), 
                          family=gaussian(),
                          data = resp3,
                          REML = TRUE)

pop.rest.poly3_MgO2.hr <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ POPULATION * poly(TEMPERATURE, 3) + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID), 
                          family=gaussian(),
                          data = resp3,
                          REML = TRUE)

AICc(pop.rest_MgO2.hr, pop.rest.poly2_MgO2.hr, pop.rest.poly3_MgO2.hr, k = 2, REML = TRUE)
check_model(pop.rest.poly3_MgO2.hr) 

#--- save model ---# 
saveRDS(pop.rest.poly3_MgO2.hr, file = "glmmTMB_rest_p3_population_MgO2_hr.RDS") 

#--- load model ---#
#pop.rest.poly3_MgO2.hr <- readRDS("glmmTMB_rest_p3_population_MgO2_hr.RDS")


#--- investigate model ---#
pop.rest.poly3_MgO2.hr %>% plot_model(type='eff',  terms=c('TEMPERATURE','POPULATION'), show.data=TRUE)
pop.rest.poly3_MgO2.hr %>% ggemmeans(~TEMPERATURE|POPULATION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
pop.rest.poly3_MgO2.hr %>% plot_model(type='est')

pop.rest.poly3_MgO2.hr %>% summary()
pop.rest.poly3_MgO2.hr %>% confint()
pop.rest.poly3_MgO2.hr  %>% r.squaredGLMM()

pop.rest.poly3_MgO2.hr %>% emtrends("POPULATION", var="TEMPERATURE") %>% pairs() %>% summary(infer=TRUE)

###########################################################################################
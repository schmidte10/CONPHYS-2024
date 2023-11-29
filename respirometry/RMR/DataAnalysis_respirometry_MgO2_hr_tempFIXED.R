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
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation")
#--- load data ---# 
resp <- read.delim("./import_files/SummaryData_2022_resp_updated.txt")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/respirometry/RMR/")

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
        RESTING_END_TIME = hms(RESTING_ENDTIME),
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
ggplot(resp3, aes(MASS,MgO2.hr_Net)) + 
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

#--- base model and resting runtime ---#
rmr.3 <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS, 
                 family=gaussian(),
                 data = resp3,
                 REML = FALSE)

AIC(rmr.1, rmr.2, rmr.3, k=2)

#--- polynomials ---# 
rmr.3.p2 <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED + RESTING_RUNTIME_SECONDS, 
                 family=gaussian(),
                 data = resp3,
                 REML = FALSE) 

rmr.3.p3 <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * poly(TEMPERATURE, 3) + MASS_CENTERED + RESTING_RUNTIME_SECONDS, 
                 family=gaussian(),
                 data = resp3,
                 REML = FALSE)

AIC(rmr.3, rmr.3.p2, rmr.3.p3, k=2)
#followed by random effects
rmr.3.p2a <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp3,
                  REML = TRUE) 


rmr.3.p2b <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|POPULATION/FISH_ID), 
                  family=gaussian(),
                  data = resp3,
                  REML = TRUE)

rmr.3.p2c <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID) + (1|POPULATION), 
                  family=gaussian(),
                  data = resp3,
                  REML = TRUE)

AIC(rmr.3.p2, rmr.3.p2a, rmr.3.p2b, rmr.3.p2c,  k=2)

#--- Final model ---# 
rmr.3a <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp3,
                  REML = TRUE)

#--- saving model ---#
saveRDS(rmr.3.p2a, file = "rmr_3_p2a.RDS") 


#--- load model ---#
rmr.3.p2 <- readRDS("rmr_3_p2a.RDS")

#--- investigate model ---#
check_model(rmr.3.p2)
pha.resid <-  rmr.3.p2 %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
rmr.3.p2a %>% DHARMa::testResiduals() 

rmr.3.p2a %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
rmr.3.p2a %>% plot_model(type='est')

rmr.3.p2a %>% summary()
rmr.3.p2a %>% Anova()
rmr.3.p2a %>% confint()
rmr.3.p2a %>% performance::r2_nakagawa()

rmr.3.p2a %>% emmeans(~ TEMPERATURE*REGION, type = "response")  %>% summary(infer=TRUE)
rmr.3.p2a %>% emmeans(pairwise ~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
rmr.3.p2a %>% emtrends(var = "TEMPERATURE", type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
#--- plot ---#
rmr.emm <- emmeans(rmr.3.p2a, ~ TEMPERATURE*REGION, 
                  at = list(TEMPERATURE = seq(from=27, to = 31.5, by=.1)))
rmr.emm.df=as.data.frame(rmr.emm)

rmr.newdata <- rmr.3.p2a %>% 
  ggemmeans(terms = c("TEMPERATURE[all]","REGION")) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

rmr.g1 <- ggplot(rmr.newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high), shape=21,
                  position=position_dodge(0.2))+
  theme_classic(); rmr.g1

predict(rmr.3a, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(rmr.3a, type='response') 
#data points based on month/situation/random effects - to get the data points
rmr.obs <-  resp3 %>% 
  mutate(Pred=predict(rmr.3.p2a, re.form=NA),
         Resid = residuals(rmr.3.p2a, type='response'),
         Fit = Pred + Resid)

rmr.g2 <- ggplot(rmr.newdata, aes(y=predicted, x=TEMPERATURE, color=group))+
  geom_jitter(data=rmr.obs, aes(y=Fit, color=REGION), width=0.05, alpha = 0.3) +
  stat_smooth(method = "lm", 
              formula =y ~ poly(x, 2, raw=TRUE)) + 
  geom_ribbon(aes(x=TEMPERATURE, ymin= conf.low, ymax= conf.high, fill = group), 
              alpha = 0.2, color=NA)+ 
  scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+ 
  scale_y_continuous(limits = c(2,12), breaks = seq(2, 12, by = 2)) +
  theme_classic() + ylab("RESTING METABOLIC RATE (MMR: MgO2/hr)")+ xlab("")+
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Low-latitude","High-latitude"),
                     name = "Regions") +
  scale_fill_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Low-latitude","High-latitude"),
                    name = "Regions") + 
  theme(legend.position = "top") + 
  annotate("text", x=31, y= 11.5, label="P =0.62", fontface="italic", size=4); rmr.g2

rmr.g2 <- ggplot(rmr.emm.df, aes(y=emmean, x=TEMPERATURE, color=REGION))+
  geom_jitter(data=rmr.obs, aes(y=Fit, color=REGION), width=0.05, alpha = 0.3) +
  stat_smooth(method = "lm", 
              formula =y ~ poly(x, 2, raw=TRUE)) + 
  geom_ribbon(aes(x=TEMPERATURE, ymin= lower.CL, ymax= upper.CL, fill = REGION), 
              alpha = 0.2, color=NA)+ 
  scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+ 
  scale_y_continuous(limits = c(2,12), breaks = seq(2, 12, by = 2)) +
  theme_classic() + ylab("RESTING METABOLIC RATE (MMR: MgO2/hr)")+ xlab("")+
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Low-latitude","High-latitude"),
                     name = "Regions") +
  scale_fill_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Low-latitude","High-latitude"),
                    name = "Regions") + 
  theme(legend.position = "top") + 
  annotate("text", x=31, y= 11.5, label="P =0.62", fontface="italic", size=4); rmr.g2

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
pop.rest <- glmmTMB(RESTING_MgO2.hr_RESPR ~ 1+ POPULATION * TEMPERATURE + MASS_CENTERED + RESTING_RUNTIME_SECONDS + (1|FISH_ID), 
                    family=gaussian(),
                    data = resp3,
                    REML = TRUE)


check_model(pop.rest) 

#--- save model ---# 
saveRDS(pop.rest.poly3_MgO2.hr, file = "glmmTMB_rest_p3_population_MgO2_hr.RDS") 

#--- load model ---#
#pop.rest.poly3_MgO2.hr <- readRDS("glmmTMB_rest_p3_population_MgO2_hr.RDS")


#--- investigate model ---#
pop.rest %>% plot_model(type='eff',  terms=c('TEMPERATURE','POPULATION'), show.data=TRUE)
pop.rest %>% ggemmeans(~TEMPERATURE|POPULATION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
pop.rest %>% plot_model(type='est')

pop.rest %>% summary()
pop.rest %>% confint()
pop.rest  %>% r.squaredGLMM()


###########################################################################################
#--- custome contrast ---# 

emm1 = emmeans(pop.rest, specs = ~ POPULATION*TEMPERATURE); emm1

#--- temperature treatment - 27C - ---#
core27 = c(0,0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
mackay27 = c(0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
chauvel27 = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 

contrast27 = emmeans(pop.rest, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core27 - mackay27, 
                                                                                           "core - chauvel" = core27 - chauvel27, 
                                                                                           "mackay - chauvel" = mackay27 -chauvel27)) %>% summary(infer=TRUE) 

#--- temperature treatment - 28.5C - ---#
core28.5 =    c(0,0,0,0,0,0,0,0,0,1/3,1/3,1/3,0,0,0,0,0,0,0,0,0,0,0,0) 
mackay28.5 =  c(0,0,0,0,0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 
chauvel28.5 = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) 

contrast28.5 = emmeans(pop.rest, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core28.5 - mackay28.5, 
                                                                                             "core - chauvel" = core28.5 - chauvel28.5, 
                                                                                             "mackay - chauvel" = mackay28.5 -chauvel28.5)) %>% summary(infer=TRUE) 

#--- temperature treatment - 30.0C - ---#
core30 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/3,1/3,1/3,0,0,0,0,0,0) 
mackay30 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,1/2,1/2,0,0,0,0,0,0,0,0,0) 
chauvel30 = c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0) 

contrast30 = emmeans(pop.rest, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core30 - mackay30, 
                                                                                           "core - chauvel" = core30 - chauvel30, 
                                                                                           "mackay - chauvel" = mackay30 -chauvel30)) %>% summary(infer=TRUE) 

#--- temperature treatment - 31.5C - ---#
core31.5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/3,1/3,1/3) 
mackay31.5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/2,1/2,0,0,0) 
chauvel31.5 = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0) 

contrast31.5 = emmeans(pop.rest, specs = ~ POPULATION*TEMPERATURE) %>% contrast(method = list("core - mackay" = core31.5 - mackay31.5, 
                                                                                             "core - chauvel" = core31.5 - chauvel31.5, 
                                                                                             "mackay - chauvel" = mackay31.5 - chauvel31.5)) %>% summary(infer=TRUE)


#--- contrasts ---# 
print("TEMPERATURE - 27"); contrast27; print("TEMPERATURE - 28.5"); contrast28.5;print("TEMPERATURE - 30"); contrast30; print("TEMPERATURE - 31.5"); contrast31.5



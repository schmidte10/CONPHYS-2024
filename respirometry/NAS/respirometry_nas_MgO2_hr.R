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
#--- set working directory ---#
# personal computer
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry")
#--- load data ---# 
resp <- read.delim("./SummaryData_2022_resp.txt")
# personal computer
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/NAS/")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/NAS")
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
         TEMPERATURE = as.numeric(TEMPERATURE),
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
           EXP_FISH_ID != "LCHA125_30" & 
           EXP_FISH_ID != "CVLA045_27")

#--- exploratory data analysis ---# 
hist(resp4$MgO2.hr_Net); shapiro.test(resp4$MgO2.hr_Net); skewness(resp4$MgO2.hr_Net) #postive/right skewed
hist(sqrt(resp4$MgO2.hr_Net)); shapiro.test(sqrt(resp4$MgO2.hr_Net)); skewness(sqrt(resp4$MgO2.hr_Net)) # sqrt transforming data seems to fix skewness

skewness(resp4$MgO2.hr_Net)
#--- make sqrt transformed column ---# 


resp4 %>% 
  group_by(REGION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(MgO2.hr_Net), 
                   Max. = max(MgO2.hr_Net), 
                   Mean = mean(MgO2.hr_Net))

#--- model formula ---# 
#net aerobic scope
MgO2.hr_NET <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED +
                 RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP +
                 (1|REGION:POPULATION) + (1|FISH_ID), 
               family=gaussian(),
               control=glmmTMBControl(optimizer=optim,
                                      optArgs = list(method='BFGS')),
               data = resp4,
               REML = FALSE)

nas_p2_MgO2.hr_NET <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED +
                         RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP +
                         (1|REGION:POPULATION) + (1|FISH_ID), 
                       family=gaussian(),
                       data = resp4,
                       control=glmmTMBControl(optimizer=optim,
                                              optArgs = list(method='BFGS')),
                       REML = FALSE) 

nas_p3_MgO2.hr_NET <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE, 3) + MASS_CENTERED +
                         RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP +
                         (1|REGION:POPULATION) + (1|FISH_ID), 
                       family=gaussian(),
                       data = resp4,
                       control=glmmTMBControl(optimizer=optim,
                                              optArgs = list(method='BFGS')),
                       REML = FALSE)

#--- inverse gamma with power link ---# 
nas_p2_MgO2.hr_NET_gamma <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED +
                                RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP +
                                (1|REGION:POPULATION) + (1|FISH_ID), 
                              family=Gamma(link = "sqrt"),
                              data = resp4,
                              control=glmmTMBControl(optimizer=optim,
                                                     optArgs = list(method='BFGS')),
                              REML = FALSE) 

nas_p3_MgO2.hr_NET_gamma <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE, 3) + MASS_CENTERED +
                                RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP +
                                (1|REGION:POPULATION) + (1|FISH_ID), 
                              family=Gamma(link = "sqrt"),
                              data = resp4,
                              control=glmmTMBControl(optimizer=optim,
                                                     optArgs = list(method='BFGS')),
                              REML = FALSE)

#--- model compairson ---#
AICc(nas_p2_MgO2.hr_NET, nas_p3_MgO2.hr_NET, 
     nas_p2_MgO2.hr_NET_gamma, nas_p3_MgO2.hr_NET_gamma,
     k = 2, REML = TRUE) 

#--- final model ---# 
nas_p2_MgO2.hr_NET <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED +
                                RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP +
                                (1|REGION:POPULATION) + (1|FISH_ID), 
                              family=gaussian(),
                              data = resp4,
                              control=glmmTMBControl(optimizer=optim,
                                                     optArgs = list(method='BFGS')),
                              REML = TRUE)

#--- saving model ---#
saveRDS(nas_p2_MgO2.hr_NET, file = "glmmTMB_MgO2_hr_Net_p2.RDS") 

#--- load model ---# 
#sqrt.MgO2.hr_NET.p3 <- readRDS("glmmTMB_sqrt.MgO2.hr_NET_p3.RDS") 

#--- investigate model ---#
#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")
check_model(nas_p2_MgO2.hr_NET)

nas_p2_MgO2.hr_NET %>% plot_model(type='eff',  terms=c('TEMPERATURE','REGION'), show.data=TRUE)
nas_p2_MgO2.hr_NET %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
nas_p2_MgO2.hr_NET %>% plot_model(type='est')

nas_p2_MgO2.hr_NET %>% summary()
nas_p2_MgO2.hr_NET %>% confint()
nas_p2_MgO2.hr_NET  %>% r.squaredGLMM()
nas_p2_MgO2.hr_NET  %>% performance::r2_nakagawa()

nas_p2_MgO2.hr_NET %>% emtrends("REGION", var="TEMPERATURE") %>% pairs() %>% summary(infer=TRUE) 
nas_p2_MgO2.hr_NET %>% emtrends(~REGION|TEMPERATURE, "TEMPERATURE", max.degree = 2) %>% pairs() %>% summary(infer=TRUE)
#--- plot ---#
newdata <- nas_p2_MgO2.hr_NET %>% ggemmeans(~TEMPERATURE|REGION) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); g1

predict(nas_p2_MgO2.hr_NET, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(nas_p2_MgO2.hr_NET, type='response') 
#data points based on month/situation/random effects - to get the data points
obs <-  resp4 %>% 
  mutate(Pred=predict(nas_p2_MgO2.hr_NET, re.form=NA),
         Resid = residuals(nas_p2_MgO2.hr_NET, type='response'),
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
  theme_classic() + ylab("NET AEROBIC SCOPE (MMR)")+
  scale_fill_manual(values=c("#DA3A36", "#0D47A1"))+ 
  scale_color_manual(values=c("#DA3A36", "#0D47A1")) + 
  theme_bw(); g2

pdf("MgO2.hr_p2_NET.pdf")
print(g2)
dev.off()

###########################################################################################
resp4 %>% 
  group_by(POPULATION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(sqrt.MgO2.hr_NET), 
                   Max. = max(sqrt.MgO2.hr_NET), 
                   Mean = mean(sqrt.MgO2.hr_NET)) %>% 
  print(n = 24)

#--- model formula ---# 
#net aerobic scope looking at differences between populations
pop.sqrt.MgO2.hr_NET <- glmmTMB(sqrt.MgO2.hr_NET ~ 1+ POPULATION * TEMPERATURE + MASS_CENTERED +
                     RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP + (1|FISH_ID), 
                   family=gaussian(),
                   #control=glmmTMBControl(optimizer=optim,
                                          #optArgs = list(method='BFGS')),
                   data = resp4,
                   REML = TRUE)


pop.sqrt.MgO2.hr_NET.p2 <- glmmTMB(sqrt.MgO2.hr_NET ~ 1+ POPULATION * poly(TEMPERATURE, 2) + MASS_CENTERED +
                        RESTING_CHAMBER + RESTING_SUMP + 
                        MAX_CHAMBER + MAX_SUMP + (1|FISH_ID), 
                      #control=glmmTMBControl(optimizer=optim,
                                             #optArgs = list(method='BFGS')),
                      family=gaussian(),
                      data = resp4,
                      REML = TRUE)

pop.sqrt.MgO2.hr_NET.p3 <- glmmTMB(sqrt.MgO2.hr_NET ~ 1+ POPULATION * poly(TEMPERATURE, 3) + MASS_CENTERED +
                        RESTING_CHAMBER + RESTING_SUMP + 
                        MAX_CHAMBER + MAX_SUMP + (1|FISH_ID), 
                      #control=glmmTMBControl(optimizer=optim,
                                             #optArgs = list(method='BFGS')),
                      family=gaussian(),
                      data = resp4,
                      REML = TRUE)

AICc(pop.sqrt.MgO2.hr_NET, pop.sqrt.MgO2.hr_NET.p2, pop.sqrt.MgO2.hr_NET.p3, k = 2, REML = TRUE) 
check_model(pop.sqrt.MgO2.hr_NET.p2) 

#--- save model ---# 
saveRDS(pop.sqrt.MgO2.hr_NET.p2, file = "glmmTMB_sqrt.MgO2.hr_NET_p2_population.RDS") 

#--- load model ---#
#pop.sqrt.MgO2.hr_NET.poly3 <- readRDS("glmmTMB_sqrt.MgO2.hr_NET_p3_population.RDS")


#--- investigate model ---#
pop.sqrt.MgO2.hr_NET.p2 %>% plot_model(type='eff',  terms=c('TEMPERATURE','POPULATION'), show.data=TRUE)
pop.sqrt.MgO2.hr_NET.p2 %>% ggemmeans(~TEMPERATURE|POPULATION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
pop.sqrt.MgO2.hr_NET.p2 %>% plot_model(type='est')

pop.sqrt.MgO2.hr_NET.p2 %>% summary()
pop.sqrt.MgO2.hr_NET.p2 %>% confint()
pop.sqrt.MgO2.hr_NET.p2  %>% r.squaredGLMM()

pop.sqrt.MgO2.hr_NET.p3 %>% emtrends("POPULATION", var="TEMPERATURE") %>% pairs() %>% summary(infer=TRUE)

#########################################################################################




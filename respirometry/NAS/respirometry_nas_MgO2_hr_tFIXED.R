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
  dplyr::rename(MASS = DRY_WEIGHT) %>% 
  mutate(MASS_CENTERED = scale(MASS, scale = FALSE, center = TRUE)) %>% 
  drop_na(MASS)

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
           EXP_FISH_ID != "LCHA124_27"& 
           EXP_FISH_ID != "LCHA135_27" & 
           #EXP_FISH_ID != "LCKM162_27" & 
           EXP_FISH_ID != "CSUD026_31.5" & 
           EXP_FISH_ID != "CTON067_28.5" & 
           EXP_FISH_ID != "CVLA054_28.5" & 
           EXP_FISH_ID != "LCHA114_28.5"& 
           EXP_FISH_ID != "LCKM163_28.5" & 
           EXP_FISH_ID != "CVLA045_28.5")
#LCHA113_30 should be removed from resting and maximum. Data trend for resting looks very strange. 

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

#--- model formula ---# 
#--- base model ---# 
nas.1 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
                              family=gaussian(),
                              data = resp4,
                              REML = FALSE) 
 
#--- experimental equipment hypothesis ---#
nas.2 <- glmmTMB(MgO2.hr_Net ~ 1+ MASS_CENTERED + RESTING_SUMP + RESTING_CHAMBER + 
                                RESTING_AM_PM + RESTING_DATE + RESTING_DATE, 
                              family=gaussian(),
                              data = resp4,
                              REML = FALSE) 

#--- base model and equip ---# 
nas.3 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + 
                   MASS_CENTERED + RESTING_CHAMBER + 
                   RESTING_AM_PM, 
                 family=gaussian(),
                 data = resp4,
                 REML = FALSE) 

AICc(nas.1, nas.2, nas.3, k=2)
#followed by random effects
nas.1 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4,
                 REML = TRUE) 

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
#sqrt.MgO2.hr_NET.p3 <- readRDS("glmmTMB_sqrt.MgO2.hr_NET_p3.RDS") 

#--- investigate model ---#
#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")
check_model(nas.1a)
pha.resid <-  nas.1a %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
nas.1a %>% DHARMa::testResiduals()

nas.1a %>% plot_model(type='eff',  terms=c('TEMPERATURE','REGION'), show.data=TRUE)
nas.1a %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
nas.1a %>% plot_model(type='est')

nas.1a %>% summary()
nas.1a %>% confint()
nas.1a  %>% r.squaredGLMM()
nas.1a  %>% performance::r2_nakagawa()

nas.1a %>% emmeans(~ TEMPERATURE*REGION, type = "response")  %>% summary(infer=TRUE)
nas.1a %>% emmeans(~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(infer=TRUE)

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
  scale_y_continuous(limits = c(2,8), breaks = seq(2, 8, by = 2)) + 
  #scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+
  theme_classic() + ylab("NET AEROBIC SCOPE (NAS: MgO2/hr)") +
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Cairns (north)","Mackay (south)"),
                     name = "Regions")+ 
  geom_signif(
    y_position = c(4.11+1.5, 5.18+1.5,5.15+1.5,4.66+1.5), xmin = c(0.8, 1.8,2.8,3.8), xmax = c(1.2,2.2,3.2,4.2),
    annotation = c("ns", "ns", "**\np =0.046", "ns"), tip_length = 0.025, color = "black"); g2

pdf("nas_1a.pdf", width = 7, height = 5)
print(g2)
dev.off()

jpeg("MgO2.hr_NET_tfixedb2.jpeg", units="in", width=7, height=5, res=300) 
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
                   REML = FALSE)


pop.sqrt.MgO2.hr_NET.p2 <- glmmTMB(sqrt.MgO2.hr_NET ~ 1+ POPULATION * poly(TEMPERATURE, 2) + MASS_CENTERED +
                        RESTING_CHAMBER + RESTING_SUMP + 
                        MAX_CHAMBER + MAX_SUMP + (1|FISH_ID), 
                      #control=glmmTMBControl(optimizer=optim,
                                             #optArgs = list(method='BFGS')),
                      family=gaussian(),
                      data = resp4,
                      REML = FALSE)

pop.sqrt.MgO2.hr_NET.p3 <- glmmTMB(sqrt.MgO2.hr_NET ~ 1+ POPULATION * poly(TEMPERATURE, 3) + MASS_CENTERED +
                        RESTING_CHAMBER + RESTING_SUMP + 
                        MAX_CHAMBER + MAX_SUMP + (1|FISH_ID), 
                      #control=glmmTMBControl(optimizer=optim,
                                             #optArgs = list(method='BFGS')),
                      family=gaussian(),
                      data = resp4,
                      REML = FALSE)

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
tran <- make.tran("power", 1/2)
MgO2.hr_NET_tfixed <- with(tran, glmmTMB(linkfun(MgO2.hr_Net) ~ 1+ REGION * TEMPERATURE + MASS_CENTERED +
                                           RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP +
                                           (1|REGION:POPULATION) + (1|FISH_ID), 
                                         family=gaussian(),
                                         data = resp4,
                                         REML = TRUE))

MgO2.hr_NET_tfixedb <- with(tran, glmmTMB(linkfun(MgO2.hr_Net) ~ 1+ REGION * TEMPERATURE + MASS_CENTERED +
                                            RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP + (1|REGION) + (1|FISH_ID), 
                                          family=gaussian(),
                                          data = resp4,
                                          REML = TRUE))

MgO2.hr_NET_tfixedc <- with(tran, glmmTMB(linkfun(MgO2.hr_Net) ~ 1+ REGION * TEMPERATURE + MASS_CENTERED +
                                            RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP + (1|POPULATION) + (1|FISH_ID), 
                                          family=gaussian(),
                                          data = resp4,
                                          REML = TRUE))

MgO2.hr_NET_tfixedd <- with(tran, glmmTMB(linkfun(MgO2.hr_Net) ~ 1+ REGION * TEMPERATURE + MASS_CENTERED +
                                            RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP + (1|FISH_ID), 
                                          family=gaussian(),
                                          data = resp4,
                                          REML = TRUE)) 

Mg



#--- model compairson ---#
AICc(MgO2.hr_NET_tfixed, MgO2.hr_NET_tfixedb, MgO2.hr_NET_tfixedc, MgO2.hr_NET_tfixedd, k=2)



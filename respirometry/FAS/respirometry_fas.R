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
library(effects)
#--- set working directory ---#
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/")

#--- load data ---# 
resp <- read.delim("./SummaryData_2022_resp.txt")
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/FAS/")

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
           EXP_FISH_ID != "LCHA125_30")

#--- exploratory data analysis ---# 
hist(resp4$FAS); shapiro.test(resp4$FAS) 

resp4 %>% 
  group_by(REGION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(NAS), 
                   Max. = max(NAS), 
                   Mean = mean(NAS))

#--- model formula ---# 
#factorial aerobic scope
fas <- glmmTMB(FAS ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + 
                 RESTING_CHAMBER + RESTING_SUMP + MAX_CHAMBER + MAX_SUMP +
                 (1|REGION:POPULATION) + (1|FISH_ID), 
               family=gaussian('inverse'),
               data = resp4,
               REML = TRUE)


fas.p2 <- glmmTMB(FAS ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED + 
                    RESTING_CHAMBER + RESTING_SUMP + 
                    MAX_CHAMBER + MAX_SUMP +
                    (1|REGION:POPULATION) + (1|FISH_ID), 
                  family=gaussian('inverse'),
                  data = resp4,
                  REML = TRUE)

fas.p3 <- glmmTMB(FAS ~ 1+ REGION * poly(TEMPERATURE, 3) + MASS_CENTERED + 
                    RESTING_CHAMBER + RESTING_SUMP + 
                    MAX_CHAMBER + MAX_SUMP +
                    (1|REGION:POPULATION) + (1|FISH_ID), 
                  family=gaussian('inverse'),
                  data = resp4,
                  REML = TRUE)

#--- model compairson ---#
AICc(fas, fas.p2, fas.p3, k = 2, REML = TRUE) 

#--- saving model ---#
saveRDS(fas.p2, file = "glmmTMB_fas_p3.RDS") 

#--- load model ---# 
#nas.p3 <- readRDS("glmmTMB_nas_p3.RDS") 

#--- investigate model ---#
#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")
check_model(fas.p2)

fas.p2 %>% plot_model(type='eff',  terms=c('TEMPERATURE','REGION'), show.data=TRUE)
fas.p2 %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
fas.p2 %>% plot_model(type='est')

fas.p2 %>% summary()
fas.p2 %>% confint()
fas.p2  %>% r.squaredGLMM()
fas.p2  %>% performance::r2_nakagawa()

fas.p2 %>% emtrends("REGION", var="TEMPERATURE") %>% pairs() %>% summary(infer=TRUE)

#--- plot ---#
newdata <- fas.p2 %>% ggemmeans(~TEMPERATURE|REGION) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); g1

predict(fas.p2, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(fas.p2, type='response') 
#data points based on month/situation/random effects - to get the data points
obs <-  resp4 %>% 
  mutate(Pred=predict(fas.p2, re.form=NA),
         Resid = residuals(fas.p2, type='response'),
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
  scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+
  theme_classic() + ylab("NET AEROBIC SCOPE (MMR)")+
  scale_fill_manual(values=c("#2f3544", "#4c8494"))+ 
  scale_color_manual(values=c("#2f3544", "#4c8494")); g2

pdf("FAS.pdf")
print(g2)
dev.off()

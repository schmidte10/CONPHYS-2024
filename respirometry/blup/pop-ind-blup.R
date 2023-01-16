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

# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/blup")


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

#--- preparing data ---# 
resp6 <-  resp4 %>% 
  mutate(fTEMPERATURE = factor(TEMPERATURE), 
         cMASS_CENTERED = scale(MASS), 
         cTEMPERATURE = scale(as.numeric(TEMPERATURE)))

ggplot(resp6, aes(x= cTEMPERATURE, y = MgO2.hr_Net, group = POPULATION)) + 
  geom_line(aes(color = POPULATION)) + ylab("net aerobic scope") + 
  xlab("Mean-centered temperature") + theme_classic()

#--- model formula ---# 
#net aerobic scope looking at differences between populations
nas.pop <- glmmTMB(MgO2.hr_Net ~ 1+ POPULATION * cTEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                   family=gaussian(),
                   data = resp6,
                   REML = TRUE)

#--- linear ---#
model1.1 <- lmer(MgO2.hr_Net ~ TEMPERATURE*POPULATION + MASS_CENTERED + (1|POPULATION), 
                 REML = FALSE, 
                 data = resp6)

summary(model1.1)
r.squaredGLMM(model1.1)

model1.1 %>% emmeans(~ TEMPERATURE*POPULATION, type = "response") %>% pairs()
model1.1 %>% emmeans(consec ~ TEMPERATURE*POPULATION, type = "response") 

temperature_pred <- data.frame(cTEMPERATURE = seq(from = min(df.pop2$cTEMPERATURE, 
                                                             to = max(df.pop2$cTEMPERATURE), 
                                                             50)))
temperature_pred$fit1.1 <- predict(model1.1, newdata = temperature_pred, re.form = NA)

ggplot(temperature_pred, aes( x = cTEMPERATURE, y = fit1.1)) + 
  geom_line(data = df.pop2, aes(y = emmean, color = POPULATION), size =1.5) + 
  geom_line(size = 2)+ 
  theme_classic()

#--- quadratic ---#
model1.2 <- lmer(MgO2.hr_Net ~ poly(cTEMPERATURE, 2, raw = TRUE)*POPULATION + (1|POPULATION), 
                 REML = FALSE, 
                 data = resp6)

summary(model1.2)
r.squaredGLMM(model1.2)

temperature_pred$fit1.2 <- predict(model1.2, newdata = temperature_pred, re.form = NA)

ggplot(temperature_pred, aes(x = cTEMPERATURE, y = fit1.2)) +
  geom_line(data = df.pop2, aes(y = emmean, colour = POPULATION)) +
  geom_line(size = 2) +
  theme_classic()

summary(model1.1)$logLik
summary(model1.2)$logLik
chi2 <- 2*(summary(model1.2)$logLik - summary(model1.1)$logLik)
1-pchisq(chi2,1)
AIC(model1.1, model1.2) 

model1.2 %>% plot_model(type='eff',  terms=c('cTEMPERATURE'), show.data=TRUE)

#--- polynomial ---#
model1.3 <- lmer(MgO2.hr_Net ~ poly(cTEMPERATURE, 3, raw = TRUE)*POPULATION + (1|POPULATION), 
                 REML = FALSE, 
                 data = resp6)

summary(model1.3)
r.squaredGLMM(model1.3)

temperature_pred$fit1.2 <- predict(model1.3, newdata = temperature_pred, re.form = NA)

ggplot(temperature_pred, aes(x = cTEMPERATURE, y = fit1.2)) +
  geom_line(data = df.pop2, aes(y = emmean, colour = POPULATION)) +
  geom_line(size = 2) +
  theme_classic()

chi2 <- 2*(summary(model1.3)$logLik - summary(model1.1)$logLik)
1-pchisq(chi2,1)



AIC(model1.1, model1.2, model1.3)

#--- adding in the population random variable does improve the model 
#--- not by much but by P = 0.049 
#--- not we will analyses the changing of random slopes 

model1.5 <- lmer(MgO2.hr_Net ~ cTEMPERATURE*POPULATION + (1+cTEMPERATURE|POPULATION), 
                 REML = FALSE, 
                 data = resp6)

summary(model1.5)
r.squaredGLMM(model1.5)

temperature_pred$fit1.5 <- predict(model1.5, newdata = temperature_pred, re.form = NA)
df.pop2$pred_pop1.5 <- predict(model1.5, re.form = NA)
df.pop2$pred_pop1.5 <- predict(model1.5, re.form = ~1(+cTEMPERATURE|POPULATION))

ggplot(temperature_pred, aes(x = cTEMPERATURE, y = fit1.5)) +
  geom_line(data = df.pop2, aes(y = pred_pop1.5, group = POPULATION, colour = POPULATION), lty = 2) +
  geom_line(data = df.pop2, aes(y = emmean, group = POPULATION, color = POPULATION), size = 1) + 
  geom_line(size = 2) +
  theme_classic()


summary(model1.3)$logLik
summary(model1.5)$logLik
chi2 <- 2*(summary(model1.5)$logLik - summary(model1.4)$logLik)
1-pchisq(chi2, 2) 

AIC(model1.1, model1.2, model1.3, model1.5)

# no significant difference between model1.4 and model1.5 
# will try to improve model by allowing population slopes to vary in curvature 

model1.6 <- lmer(MgO2.hr_Net ~ cTEMPERATURE*POPULATION + 
                   (1 + cTEMPERATURE + I(cTEMPERATURE^2)|POPULATION), 
                 REML = FALSE, data = resp6)

summary(model1.6)
r.squaredGLMM(model1.6)

summary(model1.4)$logLik
summary(model1.6)$logLik
chi2 <- 2*(summary(model1.6)$logLik - summary(model1.4)$logLik)
1-pchisq(chi2, 6) 

AIC(model1.1, model1.2, model1.3, model1.5, model1.6, nas.1a)


nas.pop <- glmmTMB(MgO2.hr_Net ~ 1+ TEMPERATURE*POPULATION + MASS_CENTERED + (1|FISH_ID), 
                   family=gaussian(),
                   data = resp4,
                   REML = FALSE)

model1.1 <- lmer(MgO2.hr_Net ~ TEMPERATURE*POPULATION + MASS_CENTERED + (1|POPULATION), 
                 REML = FALSE, 
                 data = resp6)
model1.2 <- lmer(MgO2.hr_Net ~ poly(cTEMPERATURE, 2, raw = TRUE)*POPULATION + (1|POPULATION), 
                 REML = FALSE, 
                 data = resp6)
model1.3 <- lmer(MgO2.hr_Net ~ poly(cTEMPERATURE, 3, raw = TRUE)*POPULATION + (1|POPULATION), 
                 REML = FALSE, 
                 data = resp6)
model1.5 <- lmer(MgO2.hr_Net ~ cTEMPERATURE*POPULATION + (1+cTEMPERATURE|POPULATION), 
                 REML = FALSE, 
                 data = resp6)
model1.6 <- lmer(MgO2.hr_Net ~ cTEMPERATURE*POPULATION + 
                   (1 + cTEMPERATURE + I(cTEMPERATURE^2)|POPULATION), 
                 REML = FALSE, data = resp6)

AIC(model1.1, model1.2, model1.3, model1.5, model1.6, nas.pop)

chi2 <- 2*(summary(nas.pop)$logLik - summary(model1.1)$logLik)
1-pchisq(chi2, 0) 
# best model is --- nas.pop 
####################################################################################################################
#####################################################################################################################
##########################################################################################
##########################################################################################
####################                                                    ##################
####################                 individual analysis                ##################
###################                                                     ##################
##########################################################################################

#--- preparing data ---# 
resp6 <-  resp4 %>% 
  mutate(fTEMPERATURE = factor(TEMPERATURE), 
         cMASS_CENTERED = scale(MASS), 
         cTEMPERATURE = scale(as.numeric(TEMPERATURE)))
#--- model formula ---# 
#net aerobic scope looking at differences between populations
nas.ind <- lmer(MgO2.hr_Net ~ cTEMPERATURE + (1|MASS_CENTERED), 
                data = resp6,
                REML = FALSE)
model1.1 <- nas.ind
#--- linear ---#
summary(model1.1)
r.squaredGLMM(model1.1)

model1.1 %>% emmeans(~ cTEMPERATURE, type = "response") %>% summary(infer = TRUE)

temperature_pred <- data.frame(cTEMPERATURE = seq(from = -1, to = 1, length.out = 50))
temperature_pred$fit1.1 <- predict(model1.1, newdata = temperature_pred, re.form = NA)

ggplot(temperature_pred, aes( x = cTEMPERATURE, y = fit1.1)) + 
  geom_line(data = resp6, aes(y = MgO2.hr_Net, color = FISH_ID), size =1.5) + 
  geom_line(size = 2)+ 
  theme_classic()

#--- quadratic ---#
model1.2 <- lmer(MgO2.hr_Net ~ poly(cTEMPERATURE, 2, raw = TRUE) + (1|MASS_CENTERED), 
                 REML = FALSE, 
                 data = resp6)

summary(model1.2)
r.squaredGLMM(model1.2)

temperature_pred$fit1.2 <- predict(model1.2, newdata = temperature_pred, re.form = NA)

ggplot(temperature_pred, aes(x = cTEMPERATURE, y = fit1.2)) +
  geom_line(data = resp6, aes(y = MgO2.hr_Net, colour = FISH_ID)) +
  geom_line(size = 2) +
  theme_classic()

summary(model1.1)$logLik
summary(model1.2)$logLik
chi2 <- 2*(summary(model1.2)$logLik - summary(model1.1)$logLik)
1-pchisq(chi2,1)
AIC(model1.1, model1.2) 

model1.2 %>% plot_model(type='eff',  terms=c('cTEMPERATURE'), show.data=TRUE)

#--- polynomial ---#
model1.3 <- lmer(MgO2.hr_Net ~ poly(cTEMPERATURE, 3, raw = TRUE)  + (1|MASS_CENTERED), 
                 REML = FALSE, 
                 data = resp6)

summary(model1.3)
r.squaredGLMM(model1.3)

temperature_pred$fit1.3 <- predict(model1.3, newdata = temperature_pred, re.form = NA)

ggplot(temperature_pred, aes(x = cTEMPERATURE, y = fit1.3)) +
  geom_line(data = resp6, aes(y = MgO2.hr_Net, colour = FISH_ID)) +
  geom_line(size = 2) +
  theme_classic()

summary(model1.2)$logLik
summary(model1.3)$logLik
chi2 <- 2*(summary(model1.3)$logLik - summary(model1.2)$logLik)
1-pchisq(chi2,1)



AIC(model1.1, model1.2, model1.3)

#--- adding in the population random variable does improve the model 
#--- not by much but by P = 0.049 
#--- not we will analyses the changing of random slopes 

model1.4 <- lmer(MgO2.hr_Net ~ poly(cTEMPERATURE, 2, raw = TRUE) + (1|MASS_CENTERED) + (1+cTEMPERATURE|FISH_ID), 
                 REML = FALSE, 
                 data = resp6)

summary(model1.4)
r.squaredGLMM(model1.4)

temperature_pred$fit1.4 <- predict(model1.4, newdata = temperature_pred, re.form = NA)
resp6$pred_pop1.4 <- predict(model1.4, re.form = NA)
resp6$pred_pop1.4 <- predict(model1.4, re.form = ~1(+cTEMPERATURE|FISH_ID))

ggplot(temperature_pred, aes(x = cTEMPERATURE, y = fit1.4)) +
  geom_line(data = resp6, aes(y = pred_pop1.4, group = FISH_ID, colour = FISH_ID), lty = 2) +
  #geom_line(data = resp6, aes(y = MgO2.hr_Net, group = FISH_ID, color = FISH_ID), size = 1) + 
  geom_line(size = 2) +
  theme_classic()


summary(model1.2)$logLik
summary(model1.4)$logLik
chi2 <- 2*(summary(model1.4)$logLik - summary(model1.2)$logLik)
1-pchisq(chi2, 3) 

AIC(model1.1, model1.2, model1.3, model1.4)

# no significant difference between model1.4 and model1.5 
# will try to improve model by allowing population slopes to vary in curvature 

model1.5 <- lmer(MgO2.hr_Net ~ poly(cTEMPERATURE, 2, raw = TRUE) + (1|MASS_CENTERED) +
                   (1 + cTEMPERATURE + I(cTEMPERATURE^2)|FISH_ID), 
                 REML = FALSE, data = resp6)

summary(model1.5)
r.squaredGLMM(model1.5)

summary(model1.2)$logLik
summary(model1.5)$logLik
chi2 <- 2*(summary(model1.5)$logLik - summary(model1.2)$logLik)
1-pchisq(chi2, 6) 

AIC(model1.1, model1.2, model1.3, model1.4, model1.5)

####################

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
        EXP_FISH_ID !="LCHA127_27" # deceased during experiment
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

##########################################################################################
##########################################################################################
####################                                                    ##################
####################                 population analysis                ##################
###################                                                     ##################
##########################################################################################

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
# thus indicating that there is little difference in variation explained by populations

#--- custome contrast ---# 

# may be worth investigating individual performance curves. 


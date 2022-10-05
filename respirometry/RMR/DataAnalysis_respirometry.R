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

#--- load data ---# 
resp <- read.delim("./SummaryData_2022_resp.txt")
setwd("C:/Users/Elliott/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/respirometry/RMR/")

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

#--- remove individuals where data is irregular ---# 
resp3 <- resp2 %>% 
  subset(EXP_FISH_ID !="LCHA132_27" & 
           EXP_FISH_ID != "LKES168_27" & 
           EXP_FISH_ID != "LCHA113_30" & 
           EXP_FISH_ID != "CSUD088_27" & 
           EXP_FISH_ID != "CTON062_27")

#--- exploratory data analysis: covariates ---# 
ggplot(resp3, aes(MASS, RESTING)) + 
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
ggplot(resp3, aes(REGION, RESTING)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE)

ggplot(resp3, aes(REGION, MAX)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE)

ggplot(resp3, aes(MASS, RESTING, color = REGION)) + 
  geom_point() +
  theme_classic() + 
  geom_smooth(method = "lm")

ggplot(resp3, aes(REGION, RESTING_MgO2.hr)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE)

ggplot(resp3, aes(REGION, MAX_MgO2.hr)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE) 

ggplot(resp3, aes(REGION, NAS)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE) 

ggplot(resp3, aes(REGION, FAS)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_grid(~TEMPERATURE) 

#--- exploratory data analysis ---# 
hist(resp3$RESTING); shapiro.test(resp3$RESTING)
hist(resp3$MAX); shapiro.test(resp3$MAX) 
hist(resp3$NAS); shapiro.test(resp3$NAS) 
hist(resp3$FAS); shapiro.test(resp3$FAS)

#--- model formula ---# 
#resting metablic rate
rest <- glmmTMB(RESTING ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + RESTING_CHAMBER + RESTING_SUMP + (1|REGION:POPULATION) + (1|FISH_ID), 
                family=gaussian(),
                data = resp3,
                REML = TRUE)


rest.poly2 <- glmmTMB(RESTING ~ 1+ REGION * poly(TEMPERATURE, 2) + MASS_CENTERED + RESTING_CHAMBER + RESTING_SUMP + (1|REGION:POPULATION) + (1|FISH_ID), 
                      family=gaussian(),
                      data = resp3,
                      REML = TRUE)

rest.poly3 <- glmmTMB(RESTING ~ 1+ REGION * poly(TEMPERATURE, 3) + MASS_CENTERED + RESTING_CHAMBER + RESTING_SUMP + (1|REGION:POPULATION) + (1|FISH_ID), 
                     family=gaussian(),
                     data = resp3,
                     REML = TRUE)
 
AICc(rest, rest.poly2, rest.poly3, k = 2, REML = TRUE)
check_model(rest.poly3) 
#saveRDS(rest.poly3, file = "glmmTMB_restpoly3.RDS")

#--- investigate model ---#
#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")

rest.poly3 %>% plot_model(type='eff',  terms=c('TEMPERATURE','REGION'), show.data=TRUE)
rest.poly3 %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
rest.poly3 %>% plot_model(type='est')

rest.poly3 %>% summary()
rest.poly3 %>% confint()
rest.poly3  %>% r.squaredGLMM()
rest.poly3  %>% performance::r2_nakagawa()

rest.poly3 %>% emtrends("REGION", var="TEMPERATURE") %>% pairs() %>% summary(infer=TRUE)

#--- plot ---#
newdata <- rest.poly3 %>% ggemmeans(~TEMPERATURE|REGION) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); g1

predict(rest.poly3, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(rest.poly3, type='response') 
#data points based on month/situation/random effects - to get the data points
obs <-  resp3 %>% 
  mutate(Pred=predict(rest.poly3, re.form=NA),
         Resid = residuals(rest.poly3, type='response'),
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
  theme_classic() + ylab("RESTING METABOLIC RATE (RMR)")+
  scale_fill_manual(values=c("#2f3544", "#4c8494"))+ 
  scale_color_manual(values=c("#2f3544", "#4c8494")); g2

pdf("RMR.pdf")
print(g2)
dev.off()










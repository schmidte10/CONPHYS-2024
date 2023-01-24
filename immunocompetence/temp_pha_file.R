#--- load libraries ---# 
library(tidyverse) 
library(plyr)
library(dplyr)
library(lubridate) 
library(ggplot2)
library(glmmTMB)
library(sjPlot)
library(gridExtra)
library(performance)
library(car)
library(DHARMa)
library(MuMIn)
library(rcompanion)
library(emmeans) 
library(glmmTMB) 
library(ggeffects) 
library(ggsignif)
library(rcompanion)
#--- data import ---# 
# uni computer
pha <- read.delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/import_files/pha_data.txt")
#personal computer


#--- set working directory ---# 
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/immunocompetence")

#--- preparing data ---# 
pha2 <- pha %>% 
  dplyr::rename(PHA_28.5 = PHA_285) %>%
  mutate(FISH_ID = factor(FISH_ID), 
         POPULATION = factor(POPULATION), 
         REGION = factor(REGION), 
         TANK = factor(TANK), 
         PHA_28.5 = as.numeric(PHA_28.5), 
         MASS_CENTERED = scale(MASS, center = TRUE, scale = FALSE)) %>% 
  pivot_longer(cols = c(PHA_27, 
                        PHA_28.5, 
                        PHA_30, 
                        PHA_31.5), 
               names_to = 'PHA', 
               values_to = 'IMMUNE_RESPONSE') %>% 
  separate(col = PHA, 
           into = c('TEST','TEMPERATURE'), sep = '_') %>% 
  filter(IMMUNE_RESPONSE >= 0.01) %>% # removing negative values greater than -0.05
  mutate(TEMPERATURE = factor(TEMPERATURE))

#--- initial look at data ---# 

plotNormalHistogram(pha2$IMMUNE_RESPONSE, prob = FALSE); shapiro.test(pha2$MASS) 
#plotNormalHistogram(pha2$IMMUNE_RESPONSE, prob = FALSE) # SWp = 0.002341
#plotNormalHistogram(pha2$IMMUNE_RESPONSE^(0.33), prob = FALSE); shapiro.test(pha2$IMMUNE_RESPONSE^0.33) # SWp = 0.2159

ggplot(pha2, aes(x=TEMPERATURE, y=IMMUNE_RESPONSE)) + 
  geom_violin(alpha = 0.5) +  # four potential outliers but will keep for now 
  geom_point()                # looking at histogram they don't appear to be way out of order   

ggplot(pha2, aes(x=TEMPERATURE, y=IMMUNE_RESPONSE, fill = REGION, color = REGION)) + 
  geom_violin(alpha = 0.5) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0), color = "black")

temp_leading_hist <- pha2 %>% 
  filter(REGION == "Leading")  
temp_core_hist <- pha2 %>% 
  filter(REGION == "Core")  
plotNormalHistogram(temp_core_hist$IMMUNE_RESPONSE); shapiro.test(temp_core_hist$IMMUNE_RESPONSE)
plotNormalHistogram(temp_leading_hist$IMMUNE_RESPONSE); shapiro.test(temp_leading_hist$IMMUNE_RESPONSE)

temp_leading_hist %>% ggplot(aes(x=IMMUNE_RESPONSE, fill = TEMPERATURE)) + 
  geom_density(alpha=0.5) + 
  facet_wrap(~TEMPERATURE)

temp_core_hist %>% ggplot(aes(x=IMMUNE_RESPONSE, fill = TEMPERATURE)) + 
  geom_density(alpha=0.5) + 
  facet_wrap(~TEMPERATURE)

#plotNormalHistogram((temp_core_hist$IMMUNE_RESPONSE)^(0.33)); shapiro.test((temp_core_hist$IMMUNE_RESPONSE)^(0.33))
#plotNormalHistogram((temp_leading_hist$IMMUNE_RESPONSE)^(0.33)); shapiro.test((temp_leading_hist$IMMUNE_RESPONSE)^(0.33))

pha2 %>% 
  group_by(REGION, TEMPERATURE)  %>% 
  dplyr::summarise(samples_size = n(), 
                   Min. = min(IMMUNE_RESPONSE), 
                   Max. = max(IMMUNE_RESPONSE), 
                   Mean = mean(IMMUNE_RESPONSE))

#--- model ---# 
#Gaussian model with transformed data
#tran <- make.tran("power", 1/3) 
#pha.model <- with(tran, glmmTMB(linkfun(IMMUNE_RESPONSE) ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                             #   family=gaussian(), 
                               # data = pha2, 
                               # REML = TRUE)) 

pha.model <- glmmTMB(IMMUNE_RESPONSE ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                     family=gaussian(), 
                     data = pha2, 
                     REML = TRUE) 


#--- model validation ---#
pha.model %>% check_model() 
pha.resid <-  pha.model %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
pha.model %>% DHARMa::testResiduals()

pha3 <- pha2 %>% 
  drop_na(IMMUNE_RESPONSE) 

plot(pha3$TEMPERATURE, resid(pha.model))
plot(pha3$REGION, resid(pha.model))

#--- model validation to elminate heterogenity ---# 
pha.modelb <- glmmTMB(IMMUNE_RESPONSE ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                      family=gaussian(), 
                      data = pha2,
                      dispformula = ~TEMPERATURE,
                      REML = TRUE) 

pha.modelc <- glmmTMB(IMMUNE_RESPONSE ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                      family=gaussian(), 
                      data = pha2,
                      dispformula = ~REGION,
                      REML = TRUE)

AIC(pha.modelb,pha.modelc)

pha.modelb %>% check_model() 

#--- model still not performing well ---# 
pha.model.gamma <- glmmTMB(IMMUNE_RESPONSE ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                     family=Gamma(link="log"), 
                     data = pha2, 
                     REML = TRUE) 

pha.model.gamma %>% check_model() 
pha.resid <-  pha.model.gamma %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
pha.model.gamma %>% DHARMa::testResiduals()

#--- looks better but not great ---# 

pha.model.gamma <- glmmTMB(IMMUNE_RESPONSE ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                           family=Gamma(link="inverse"), 
                           data = pha2, 
                           REML = TRUE) 

pha.model.gamma %>% check_model() 
pha.resid <-  pha.model.gamma %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
pha.model %>% DHARMa::testResiduals()


#--- looks better but not great ---# 

pha.model.tweedie <- glmmTMB(IMMUNE_RESPONSE ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                           family=tweedie(), 
                           data = pha2, 
                           REML = TRUE) 

pha.model.tweedie %>% check_model() 
pha.resid <-  pha.model.tweedie %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
pha.model.tweedie %>% DHARMa::testResiduals()
#--- partial plots ---# 
pha.model.gamma %>% ggemmeans(~TEMPERATURE*REGION) %>% plot()
pha.model.gamma %>% summary()
pha.model.gamma %>% confint()
pha.model.gamma %>% performance::r2()
pha.model.gamma %>% performance::r2_nakagawa()


#--- investigation on which random factors to include ---#
pha.glmm.gamma <- glmmTMB(IMMUNE_RESPONSE ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                                family=Gamma('log'), 
                                data = pha2, 
                                REML = TRUE) 

pha.glmm.gamma.b <- glmmTMB(IMMUNE_RESPONSE ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1|POPULATION/FISH_ID), 
                                 family=Gamma('log'), 
                                 data = pha2,
                                 REML = TRUE) 


pha.glmm.gamma.c <- glmmTMB(IMMUNE_RESPONSE ~ 1 + REGION * TEMPERATURE + MASS_CENTERED + (1 + REGION|POPULATION) + (1|FISH_ID), 
                                 family=Gamma('log'), 
                                 data = pha2,
                            control=glmmTMBControl(optimizer=optim,
                                                   optArgs = list(method='BFGS')),
                            
                                 REML = TRUE)


AIC(pha.glmm.gamma, pha.glmm.gamma.b, pha.glmm.gamma.c, k=2)

#--- save model ---# 
saveRDS(pha.glmm.gamma.b, "pha_gamma_b.RDS")

#--- model validation ---#
pha.glmm.gamma.b %>% check_model() 
pha.resid <-  pha.glmm.gamma.b %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
pha.glmm.gamma %>% DHARMa::testResiduals()

#--- partial plots ---# 
pha.glmm.gamma.b %>% ggemmeans(~TEMPERATURE*REGION) %>% plot()
pha.glmm.gamma.b %>% summary()
pha.glmm.gamma.b %>% confint()
pha.glmm.gamma.b %>% performance::r2_nakagawa()

#--- Results ---#
pha.glmm.gamma.b %>% emmeans(~ TEMPERATURE*REGION, type = "response") %>% regrid() %>% pairs(by = "TEMPERATURE") %>% summary(infer=TRUE) 
pha.glmm.gamma.b %>% emmeans(~ TEMPERATURE*REGION, type = "response") %>% regrid() %>% pairs(by = "REGION") %>% summary(infer=TRUE) 
#--- plot ---# 
newdata <- pha.glmm.gamma.b %>% ggemmeans(~TEMPERATURE|REGION) %>% 
  as.data.frame() %>% 
  dplyr::rename(TEMPERATURE = x) 

g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color = group)) + 
  geom_point() + 
  theme_classic(); g1

# predict(pha.glmm.gamma, re.form=NA)
# residuals(pha.glmm.gamma, type = "response") 

obs <- pha2 %>% 
  mutate(Pred = predict(pha.glmm.gamma.b, re.form=NA), 
         Resid = residuals(pha.glmm.gamma.b, type = 'response'), 
         Fit = Pred - Resid)

g2 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color = group)) + 
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape=19, 
                  size=1, 
                  position = position_dodge(0.2)) + 
  scale_y_continuous(limits = c(0,0.9), breaks = seq(0, 0.9, by =0.15)) + 
  theme_classic() + ylab("PHA Swelling response (mm)") + 
  scale_color_manual(values=c("#DA3A36","#0D47A1"), labels = c("Cairns (north)","Mackay (south)"),
                     name = "Regions"); g2
  #geom_signif( 
    #y_position = c(0.305+0.05, 0.55+0.05, 0.267+0.05, 0.125+0.05), 
    #xmin = c(0.8, 1.8, 2.8, 3.8), 
    #xmax = c(1.2, 2.2, 3.2, 4.2), 
    #annotation = c("ns", "ns", "ns", "ns"), 
    #tip_length = 0.025, 
    #color = "black")

pdf("pha_figure.pdf", width  = 7, height = 5)
print(g2)
dev.off()

jpeg("pha_figure.jpeg", units="in", width=7, height=5, res=300)
print(g2)
dev.off()

#--- set working directory ---# 
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/hematocrit-hemoglobin") 

#--- load libraries ---# 
library(tidyverse)
library(readr)
library(ggplot2)
library(glmmTMB)
library(MuMIn)
library(performance)
library(DHARMa)
library(emmeans) 
library(ggeffects)
library(sjPlot)

#--- load data ---# 
hema <- read_delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/import_files/HematocritHemoglobin.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE) %>% 
  select(c(1:9))

#--- Data Manipulation ---# 
hema <-  hema %>% 
  mutate(PERC_RBC = as.numeric(PERC_RBC), 
         MASS = as.numeric(MASS),
         MASS_CENTERED = scale(MASS, center = TRUE, scale = FALSE)) %>% 
  drop_na(PERC_RBC)
#--- exploratory analysis ---# 
ggplot(hema, aes(y=PERC_RBC, x=REGION)) + 
  geom_boxplot() + 
  theme_classic() 

hist(hema$PERC_RBC)

hema %>% ggplot(aes(x=PERC_RBC)) + 
  geom_density() +
  facet_wrap(~REGION)

#--- model exploration ---# 
#--- best model with fixed effects ---#
model.1 <- glmmTMB(PERC_RBC ~ REGION, 
                family = gaussian(), 
                REML = FALSE, 
                data = hema)

model.2 <- glmmTMB(PERC_RBC ~ REGION + MASS_CENTERED, 
                 family = gaussian(), 
                 REML = FALSE, 
                 data = hema) 

AICc(model.1, model.2, k=2) 

#--- best model with random effects ---# 
model.1 <- glmmTMB(PERC_RBC ~ REGION, 
                   family = gaussian(), 
                   REML = TRUE, 
                   data = hema) 

model.1a <- glmmTMB(PERC_RBC ~ REGION + (1|POPULATION), 
                    family = gaussian(), 
                    REML = TRUE, 
                    data = hema) 

model.1b <- glmmTMB(PERC_RBC ~ REGION + (REGION|POPULATION), 
                    family = gaussian(), 
                    REML = TRUE, 
                    data = hema) 

AIC(model.1, model.1a, model.1b, k=2)

#--- save model ---# 
saveRDS(model.1, "hema_model_1.RDS")
# model 1 is the best 
model.1 <- readRDS("./hema_model_1.RDS")
#--- checking model performance ---#
model.1 %>% check_model()
model.1.reside <- model.1 %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE)
model.1 %>% DHARMa::testResiduals()  

#---plotting model estimates ---# 
model.1 %>% summary()
model.1 %>% confint()
model.1 %>% r.squaredGLMM()

#--- results ---# 
model.1 %>% emmeans(~ REGION, type = "response") %>% pairs() %>% summary(infer=TRUE) 

#--- plot ---# 
hema.newdata <-  model.1 %>% ggemmeans(~REGION) %>% 
  as.data.frame() %>% 
  dplyr::rename(REGION = x)

hema.plot1 <- ggplot(hema.newdata, aes(y=predicted, x=REGION)) + 
  geom_point() + 
  theme_classic(); hema.plot1

obs <- hema %>% drop_na(MASS) %>%
  mutate(Pred = predict(model.1, re.form=NA), 
         Resid = residuals(model.1, type = "response"), 
         Fit = Pred + Resid)
 
hematocrit.plot <- ggplot(hema.newdata, aes(y=predicted, x=REGION, color=REGION))  + 
  geom_jitter(data=obs, aes(y=Pred, x=REGION, color =REGION), 
              width = 0.05, alpha=0.3)+
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape = 19, 
                  size = 1, 
                  position = position_dodge(0.2)) + 
  scale_color_manual(values=c("#B2182B", "#4393C3"), labels = c("Low","High"),
                     name = "Latitude") +
  ylab("HEMATOCRIT RATIO") +
  scale_x_discrete(name = "LATITUDE", 
                   labels = c("Low","High"))+
  theme_classic() + 
  theme(legend.position = 'none'); hematocrit.plot


pdf("hematocrit_plot.pdf", width = 7, height = 5)
print(hematocrit.plot)
dev.off()

hema.plot2 <- plot_grid(hematocrit.plot); hema.plot2

ggsave("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/supplemental_figures/Supplemental_figure4.pdf", width = 18, height = 13, units = 'cm', dpi = 360)
hema.plot2 
dev.off()

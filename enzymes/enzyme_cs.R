#--- libraries ---# 
library(plyr)
library(tidyverse)
library(chron) 
library(lubridate)
library(janitor) 
library(ggpubr) 
library(MuMIn)
library(glmmTMB)
library(performance)
library(DHARMa)
library(emmeans) 
library(ggeffects)
library(cowplot)
#--- set working directory ---# 
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/enzymes") 

#--- import data ---# 
cs <- read_delim("CS_LocalAdapt6.txt", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(...21 = col_skip(), 
                                              ...22 = col_skip()), trim_ws = TRUE) %>% 
  clean_names() %>% 
  mutate(creation_time = as.POSIXct(creation_time, format = "%d/%m/%Y %H:%M:%S"))
tissue.mass <- read.delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/enzymes/tissue_mass.txt") %>% 
  dplyr::rename(FISH_ID = fish_id)


cs2 <- cs %>%
  mutate(muscle_type = str_replace(muscle_type, " ", ".")) %>%
  unite("UNIQUE_SAMPLE_ID", c(fish_id,temperature,sample_index), sep="_", remove = FALSE) %>% 
  separate(creation_time, into=c('DATE','TIME'), sep = " ", remove = FALSE) %>% 
  arrange(sample_id_1, DATE, TIME) 

cs3 <- cs2 %>% 
  mutate(DATE = as.Date(creation_time), 
         TIME = format(creation_time, "%H:%M:%S")) %>%
  mutate(TIME = hms(cs2$TIME)) %>% 
  mutate(TIME = chron(times=cs2$TIME)) %>% 
  arrange(TIME) %>%
  group_by(UNIQUE_SAMPLE_ID, sample_id_1) %>% 
  mutate(TIME_DIFF = TIME - first(TIME)) %>% 
  filter(TIME != first(TIME)) %>%
  ungroup() %>% 
  mutate(TIME_DIFF_SECS = period_to_seconds(hms(TIME_DIFF))) %>% 
  mutate(MINUTES = TIME_DIFF_SECS/60) %>% 
  mutate(MINUTES = round(MINUTES, digits = 2)) %>% 
  dplyr::rename(CUVETTE = sample_id_1) %>% 
  mutate(REGION = substr(fish_id, 1, 1 ), 
         POPULATION = substr(fish_id, 2, 4), 
         SAMPLE_NO = substr(fish_id, 5, 7)) %>% 
  mutate(REGION = case_when( REGION =="L"~ "Leading", 
                             REGION == "C" ~ "Core", 
                             TRUE ~ "na")) 

#---- filter out samples ---# 
cs3.filtered <- cs3 %>% 
  dplyr::rename(TEMPERATURE = temperature, 
         FISH_ID = fish_id) %>%
  filter(!(TEMPERATURE == "50" & FISH_ID == "LCKM158")) %>% 
  filter(!(TEMPERATURE == "50" & FISH_ID == "CSUD010")) %>% 
  filter(!(TEMPERATURE == "40" & FISH_ID == "CSUD018")) %>% 
  filter(!(TEMPERATURE == "20" & FISH_ID == "CTON061")) %>% 
  filter(!(TEMPERATURE == "50" & FISH_ID == "CTON065")) %>%
  filter(!(FISH_ID == "CTON069")) %>% 
  filter(!(UNIQUE_SAMPLE_ID %in% c("CTON060_50_3", 
                                   "CTON061_30_3", 
                                   "CTON061_40_3", 
                                   "CTON061_50_2", 
                                   "CTON061_50_3")))%>% 
  mutate(UNIQUE_SAMPLE_ID = str_sub(UNIQUE_SAMPLE_ID, end = -3))

#--- Template - use for data analysis---#
CS_activity <- cs3.filtered %>% 
  group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
  do({
    mod = lm(result ~ MINUTES, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2], 
               r2 = summary(mod)$adj.r.squared)
  }) %>%
  ungroup() %>%
  filter(CUVETTE != ("6"))%>% 
  filter(CUVETTE != ("4"))%>% 
  filter(CUVETTE != ("5"))#%>% 
  #mutate(UNIQUE_SAMPLE_ID = str_sub(UNIQUE_SAMPLE_ID, end = -3))

CS_activity_means <- CS_activity %>% 
  group_by(UNIQUE_SAMPLE_ID) %>% 
  mutate(Mean = mean(Slope))

distinct(CS_activity_means[,c(1,5)]) 

CS_background <- cs3.filtered %>% 
  group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
  do({
    mod = lm(result ~ MINUTES, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2], 
               r2 = summary(mod)$adj.r.squared)
  }) %>%
  ungroup() %>%
  filter(CUVETTE == ("5")) %>% 
  dplyr::rename(Background = Slope) #%>% 
  #mutate(UNIQUE_SAMPLE_ID = str_sub(UNIQUE_SAMPLE_ID, end = -2))

final_table <- CS_activity %>% 
  full_join(distinct(CS_activity_means[,c(1,6)]), by = "UNIQUE_SAMPLE_ID") %>% 
  full_join(CS_background[,c(1,4)], by = "UNIQUE_SAMPLE_ID") 
final_table$Mean[duplicated(final_table$Mean)] <- ""
final_table$Background[duplicated(final_table$Background)] <- ""
final_table <- final_table %>% 
  mutate(Mean = as.numeric(Mean), 
         Background = as.numeric(Background), 
         Background_perc = Background/Mean) %>% 
  drop_na()


CS.data <- final_table %>% 
  select(c(UNIQUE_SAMPLE_ID, Mean, Background, Background_perc)) %>% 
  mutate(Mean = as.numeric(Mean), 
         Background = as.numeric(Background), 
         Background_perc = as.numeric(Background_perc)) %>% 
  mutate(Background2 = case_when(Background_perc <= 0.05 ~ 0, 
                                 TRUE ~ Background), 
         CS_ABSORBANCE = Mean - Background2) %>%
  inner_join(select(cs3.filtered, c(UNIQUE_SAMPLE_ID, REGION, POPULATION, TEMPERATURE, FISH_ID)), by ="UNIQUE_SAMPLE_ID") %>% 
  inner_join(tissue.mass, by = "FISH_ID") %>% 
  mutate(TISSUE_MASS_CENTERED = scale(TISSUE_MASS, center = TRUE, scale = FALSE)) %>%
  distinct(UNIQUE_SAMPLE_ID, REGION, POPULATION, .keep_all = TRUE) %>% 
  mutate(temperature = factor(TEMPERATURE), 
         PATH_LENGTH = 1, 
         EXTINCTION_COEFFICIENT = 13.6, 
         TISSUE_CONCENTRATION = 0.01, 
         ASSAY_VOL = 0.930, 
         SAMPLE_VOL = 0.020, 
         CS_ACTIVITY = ((CS_ABSORBANCE/(PATH_LENGTH*EXTINCTION_COEFFICIENT*TISSUE_CONCENTRATION))*(ASSAY_VOL/SAMPLE_VOL))*-1)  
#filter(LDH_ACTIVITY >= 0)

#--- modelling data ---# 
CS.summary.table <- CS.data %>% 
  group_by(FISH_ID) %>% 
  summarise(#sample_size = count(), 
    Min. = min(CS_ACTIVITY), 
    Max. = max(CS_ACTIVITY), 
    Mean = mean(CS_ACTIVITY)) 

#--- models - fixed factor ---# 
cs.model.1 <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE + TISSUE_MASS_CENTERED, 
                       family=gaussian(), 
                       data = CS.data, 
                       REML = TRUE)  

cs.model.2 <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*TEMPERATURE, 
                       family=gaussian(), 
                       data = CS.data, 
                       REML = TRUE)  

AIC(cs.model.1, cs.model.2, k=2)

#--- models ---# 
cs.model.1 <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                       family=gaussian(), 
                       data = CS.data, 
                       REML = TRUE) 

cs.model.2 <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|POPULATION/FISH_ID), 
                       family=gaussian(), 
                       data = CS.data,
                       control=glmmTMBControl(optimizer=optim, 
                                              optArgs = list(method='BFGS')), 
                       REML = TRUE) 

cs.model.3 <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|FISH_ID) + (1 + REGION|POPULATION), 
                       family=gaussian(), 
                       data = CS.data,
                       control=glmmTMBControl(optimizer=optim, 
                                              optArgs = list(method='BFGS')), 
                       REML = TRUE)

+ (1|POPULATION/fish_id)

#control=glmmTMBControl(optimizer=optim,
#optArgs = list(method='BFGS')),


#--- Model comparison ---# 
AIC(cs.model.1, cs.model.2, k=2)

#---final model ---# 
cs.model.2 <- glmmTMB(CS_ACTIVITY ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|POPULATION/FISH_ID), 
                      family=gaussian(), 
                      data = CS.data, 
                      REML = TRUE) 

#--- save model ---# 
saveRDS(cs.model.2, "cs_model_2.RDS")

#--- model investigation ---# 
cs.model.2 %>% check_model()
cs.model.2 %>% simulateResiduals(plot = TRUE, integerResponse = TRUE) 
cs.model.2 %>% testResiduals()

#--- partial plots ---# 
cs.model.2 %>% ggemmeans(~temperature*REGION) %>% plot()
cs.model.2 %>% summary()
cs.model.2 %>% confint()
cs.model.2 %>% performance::r2()

#--- results ---# 
cs.model.2 %>% emmeans(~ temperature*REGION, type = "response") %>% pairs(by = "temperature") %>% summary(infer = TRUE) 
cs.model.2 %>% emmeans(~ temperature*REGION, type = "response")  %>% summary(infer = TRUE) 
#--- plot ---# 
newdata <- cs.model.2 %>% ggemmeans(~temperature|REGION) %>% 
  as.data.frame() %>% 
  rename(TEMPERATURE = x)

g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color = group)) + 
  geom_point() + 
  theme_classic(); g1

obs <- CS.data %>% 
  mutate(Pred = predict(cs.model.2, re.form=NA), 
         Resid = residuals(cs.model.2, type = 'response'), 
         Fit = Pred - Resid)

cs.plot2 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color = group)) + 
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape=19, 
                  size=1, 
                  position = position_dodge(0.2)) + 
  #scale_y_continuous(limits = c(0,0.9), breaks = seq(0, 0.9, by =0.15)) + 
  theme_classic() + ylab("CS ACTIVITY SLOPE") + 
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Low","High"),
                     name = "Latitude") + 
  scale_y_continuous(limits=c(-450,0), breaks = seq(-450,0,100))+
  theme(legend.position = c(.8, .90)); cs.plot2

pdf("cs.pdf", width= 7, height = 5)
print(g2)
dev.off()

jpeg("cs.jpeg", units="in", width=7, height=5, res=300) 
print(g2)
dev.off()

#--- general relationship between LDH and temperature (all fish) ---#
cs.data2 <- CS.data %>% 
  mutate(temperature = as.numeric(temperature))
cs.cmb <- glmmTMB(CS_ACTIVITY ~ 1 + temperature + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                   family=gaussian(), 
                   data = cs.data2, 
                   REML = TRUE)

cs.cmb.2 <- glmmTMB(CS_ACTIVITY ~ 1 + poly(temperature, 2) + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                     family=gaussian(), 
                     data = cs.data2, 
                     REML = TRUE) 

cs.cmb.3 <- glmmTMB(CS_ACTIVITY ~ 1 + poly(temperature, 3) + TISSUE_MASS_CENTERED + (1|FISH_ID), 
                     family=gaussian(), 
                     data = cs.data2, 
                     REML = TRUE) 

AIC(cs.cmb, cs.cmb.2, cs.cmb.3, k=2)

saveRDS(cs.cmb.3, "cs_cmb_3.RDS")

summary(cs.cmb.3)
cs.cmb.3 %>% check_model()
cs.cmb.3 %>% emtrends(~ temperature, var = "temperature") %>% summary(infer = TRUE)
cs.cmb.3 %>% performance::r2_nakagawa()
####### 
load("./ldh.plot.RData")
enzyme.plot <- ggarrange(ldh2, cs.plot2, labels = c("A","B")) ; enzyme.plot
enzyme.plot2 <- plot_grid(ldh2, cs.plot2, align = "h", axis="bt", rel_widths = c(1,1.1), labels = c("A","B")); enzyme.plot2

ggsave("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/figures/figure4.pdf", width = 20, height = 12, units = 'cm', dpi = 360)
enzyme.plot2 
dev.off()

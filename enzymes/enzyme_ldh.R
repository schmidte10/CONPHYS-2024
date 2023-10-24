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
#--- set working directory ---# 
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Local_adaptation/Chapter1_LocalAdaptation/enzymes") 

#--- import data ---# 
ldh <- read_delim("LDH_LocalAdapt.txt", delim = "\t", 
                  escape_double = FALSE, col_types = cols(`Creation time` = col_datetime(format = "%d/%m/%Y %H:%M")), 
                  trim_ws = TRUE)
tissue.mass <- read.delim("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/enzymes/tissue_mass.txt")
#--- data preparation/manipulation ---# 
ldh2 <- ldh %>%
  clean_names() %>%
  mutate(muscle_type = str_replace(muscle_type, " ", ".")) %>%
  unite("UNIQUE_SAMPLE_ID", c(fish_id,temperature,sample_index), sep="_", remove = FALSE) %>% 
  separate(creation_time, into=c('DATE','TIME'), sep = " ", remove = FALSE) %>% 
  arrange(sample_id_1, DATE, TIME) 

ldh3 <- ldh2 %>% 
  mutate(TIME = hms(ldh2$TIME)) %>% 
  mutate(TIME = chron(times=ldh2$TIME)) %>% 
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

#--- filter out samples that need do not pass quailty check ---# 
grp1 <- c("CSUD008_20_1","CSUD008_20_2","CSUD008_20_3","CSUD008_20_4","CSUD008_20_5","CSUD008_20_6", 
          "CVLA047_50_1","CVLA047_50_2","CVLA047_50_3","CVLA047_50_4","CVLA047_50_5","CVLA047_50_6", 
          "CVLA046_50_1","CVLA046_50_2","CVLA046_50_3","CVLA046_50_4","CVLA046_50_5","CVLA046_50_6") 
grp2 <- c("LCKM180_30_1","LCKM180_30_2","LCKM180_30_3","LCKM180_30_4","LCKM180_30_5","LCKM180_30_6", 
          "LKES172_50_1","LKES172_50_2","CLKES172_50_3","LKES172_50_4","LKES172_50_5","LKES172_50_6", 
          "LCHA114_50_1","LCHA114_50_2","LCHA114_50_3","LCHA114_50_4","LCHA114_50_5","LCHA114_50_6", 
          "CSUD074_50_1","CSUD074_50_2","CSUD074_50_3","CSUD074_50_4","CSUD074_50_5","CSUD074_50_6")
grp3 <- c("LCKM165_50_1","LCKM165_50_2","LCKM165_50_3","LCKM165_50_4","LCKM165_50_5","LCKM165_50_6", 
          "LCKM163_50_1","LCKM163_50_2","CLCKM163_50_3","LCKM163_50_4","LCKM163_50_5","LCKM163_50_6", 
          "CTON068_50_1","CTON068_50_2","CTON068_50_3","CTON068_50_4","CTON068_50_5","CTON068_50_6", 
          "CVLA104_50_1","CVLA104_50_2","CVLA104_50_3","CVLA104_50_4","CVLA104_50_5","CVLA104_50_6") 
grp4 <- c("LCHA135_50_1","LCHA135_50_2","LCHA135_50_3","LCHA135_50_4","LCHA135_50_5","LCHA135_50_6", 
          "CTON069_50_1","CTON069_50_2","CCTON069_50_3","CTON069_50_4","CTON069_50_5","CTON069_50_6", 
          "CVLA045_50_1","CVLA045_50_2","CVLA045_50_3","CVLA045_50_4","CVLA045_50_5","CVLA045_50_6") 
grp5 <- c("CSUD014_50_1","CSUD014_50_2","CSUD014_50_3","CSUD014_50_4","CSUD014_50_5","CSUD014_50_6", 
          "CTON110_50_1","CTON110_50_2","CCTON110_50_3","CTON110_50_4","CTON110_50_5","CTON110_50_6")

ldh3.filtered <- ldh3 %>% 
  filter(!(UNIQUE_SAMPLE_ID %in% c("LCKM154_20_1", 
                                   "LKES143_30_3", 
                                   "LKES143_20_2", 
                                   "CSUD010_40_2"))) %>% 
  group_by(UNIQUE_SAMPLE_ID) %>% 
  arrange(UNIQUE_SAMPLE_ID, TIME) %>% 
  filter(!(UNIQUE_SAMPLE_ID %in% grp1 & row_number() > (n() - 1))) %>% 
  filter(!(UNIQUE_SAMPLE_ID %in% grp2 & row_number() > (n() - 2))) %>% 
  filter(!(UNIQUE_SAMPLE_ID %in% grp3 & row_number() > (n() - 3))) %>% 
  filter(!(UNIQUE_SAMPLE_ID %in% grp4 & row_number() > (n() - 4))) %>% 
  filter(!(UNIQUE_SAMPLE_ID %in% grp5 & row_number() > (n() - 5))) %>% 
  ungroup() %>% 
  mutate(UNIQUE_SAMPLE_ID = str_sub(UNIQUE_SAMPLE_ID, end = -3))

#--- Template - use for data analysis---#
LDH_activity <- ldh3.filtered %>% 
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
  filter(CUVETTE != ("5"))

LDH_activity_means <- LDH_activity %>% 
  group_by(UNIQUE_SAMPLE_ID) %>% 
  mutate(Mean = mean(Slope)) 

distinct(LDH_activity_means[,c(1,5)]) 

LDH_background <- ldh3 %>% 
  group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
  do({
    mod = lm(result ~ MINUTES, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  }) %>%
  ungroup() %>%
  filter(CUVETTE == ("5")) %>% 
  dplyr::rename(Background = Slope) %>% 
  mutate(UNIQUE_SAMPLE_ID = str_sub(UNIQUE_SAMPLE_ID, end = -3))

final_table <- LDH_activity %>% 
  full_join(distinct(LDH_activity_means[,c(1,6)]), by = "UNIQUE_SAMPLE_ID") %>% 
  full_join(LDH_background[,c(1,4)], by = "UNIQUE_SAMPLE_ID") 
final_table$Mean[duplicated(final_table$Mean)] <- ""
final_table$Background[duplicated(final_table$Background)] <- ""
final_table <- final_table %>% 
  mutate(Mean = as.numeric(Mean), 
         Background = as.numeric(Background), 
         Background_perc = Background/Mean) 


ldh.data <- final_table %>% 
  select(c(UNIQUE_SAMPLE_ID, Mean, Background, Background_perc)) %>% 
  mutate(Mean = as.numeric(Mean), 
         Background = as.numeric(Background), 
         Background_perc = as.numeric(Background_perc)) %>% 
  mutate(Background2 = case_when(Background_perc <= 0.05 ~ 0, 
                                    TRUE ~ Background), 
         LDH_ABSORBANCE = Mean - Background2) %>%
  drop_na() %>% 
  inner_join(select(ldh3.filtered, c(UNIQUE_SAMPLE_ID, REGION, POPULATION, temperature, fish_id)), by ="UNIQUE_SAMPLE_ID") %>% 
  inner_join(tissue.mass, by = "fish_id") %>% 
  mutate(TISSUE_MASS_CENTERED = scale(TISSUE_MASS, center = TRUE, scale = FALSE)) %>%
  distinct(UNIQUE_SAMPLE_ID, REGION, POPULATION, .keep_all = TRUE) %>% 
  mutate(temperature = factor(temperature), 
         PATH_LENGTH = 1, 
         EXTINCTION_COEFFICIENT = 6.22, 
         TISSUE_CONCENTRATION = 0.2, 
         ASSAY_VOL = 2.975, 
         SAMPLE_VOL = 0.025, 
         LDH_ACTIVITY = ((LDH_ABSORBANCE/(PATH_LENGTH*EXTINCTION_COEFFICIENT*TISSUE_CONCENTRATION))*(ASSAY_VOL/SAMPLE_VOL))*-1)  
  #filter(LDH_ACTIVITY >= 0)

#--- quailty check ---# 
# use app to compelte quailty check 


  
ggplot(ldh.data, aes(x =as.numeric(temperature), y= LDH_ACTIVITY, color = POPULATION)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)
#--- begin data analysis ---# 
ggplot(ldh.data, aes(x = LDH_ACTIVITY, fill = temperature, color = temperature)) + 
  geom_density(alpha =0.5, position = "identity") 

# data is not normal - perhaps more samples will end up helping 

ldh.data %>% 
  group_by(fish_id) %>% 
  summarise(#sample_size = count(), 
            Min. = min(LDH_ACTIVITY), 
            Max. = max(LDH_ACTIVITY), 
            Mean = mean(LDH_ACTIVITY)) 

#--- models - fixed factor ---# 
ldh.model.1 <- glmmTMB(LDH_ACTIVITY ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED, 
                       family=gaussian(), 
                       data = ldh.data, 
                       REML = TRUE)  

ldh.model.2 <- glmmTMB(LDH_ACTIVITY ~ 1 + REGION*temperature, 
                       family=gaussian(), 
                       data = ldh.data, 
                       REML = TRUE)  

AIC(ldh.model.1, ldh.model.2, k=2)

#--- models ---# 
ldh.model.1 <- glmmTMB(LDH_ACTIVITY ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|fish_id), 
                       family=gaussian(), 
                       data = ldh.data, 
                       REML = TRUE) 

ldh.model.2 <- glmmTMB(LDH_ACTIVITY ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|POPULATION/fish_id), 
                  family=gaussian(), 
                  data = ldh.data,
                  control=glmmTMBControl(optimizer=optim, 
                                         optArgs = list(method='BFGS')), 
                  REML = TRUE) 

ldh.model.3 <- glmmTMB(LDH_ACTIVITY ~ 1 + REGION*temperature + TISSUE_MASS_CENTERED + (1|fish_id) + (1 + REGION|POPULATION), 
                       family=gaussian(), 
                       data = ldh.data,
                       control=glmmTMBControl(optimizer=optim, 
                                              optArgs = list(method='BFGS')), 
                       REML = TRUE)

#control=glmmTMBControl(optimizer=optim,
#optArgs = list(method='BFGS')),


#--- Model comparison ---# 
AIC(ldh.model.1, ldh.model.2, ldh.model.3, k=2)

#--- save model ---# 
saveRDS(ldh.model.2, "ldh_model_2.RDS")

#--- model investigation ---# 
ldh.model.2 %>% check_model()
ldh.model.2 %>% simulateResiduals(plot = TRUE, integerResponse = TRUE) 
ldh.model.2 %>% testResiduals()

#--- partial plots ---# 
ldh.model.2 %>% ggemmeans(~temperature*REGION) %>% plot()
ldh.model.2 %>% summary()
ldh.model.2 %>% confint()
ldh.model.2 %>% performance::r2()

#--- results ---# 
ldh.model.2 %>% emmeans(~ temperature*REGION, type = "response") %>% pairs(by = "temperature") %>% summary(infer = TRUE) 
ldh.model.2 %>% emmeans(~ temperature*REGION, type = "response")  %>% summary(infer = TRUE) 
#--- plot ---# 
ldh.newdata <- ldh.model.2 %>% ggemmeans(~temperature|REGION) %>% 
  as.data.frame() %>% 
  rename(TEMPERATURE = x)
  
ldh.obs <- ldh.data %>% 
  mutate(Pred = predict(ldh.model.2, re.form=NA), 
         Resid = residuals(ldh.model.2, type = 'response'), 
         Fit = Pred - Resid)

ldh2 <- ggplot(ldh.newdata, aes(y=predicted, x=TEMPERATURE, color = group)) + 
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape=19, 
                  size=1, 
                  position = position_dodge(0.2)) + 
  #scale_y_continuous(limits = c(0,0.9), breaks = seq(0, 0.9, by =0.15)) + 
  theme_classic() + ylab("LDH activity slope") + 
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Low-latitude","High-latitude"),
                     name = "Regions") +
  scale_y_continuous(limits=c(0,7), breaks = seq(0,6,1.5))+
  theme(legend.position = 'none'); ldh2

pdf("LDH.pdf", width= 7, height = 5)
print(ldh2)
dev.off()

jpeg("LDH.jpeg", units="in", width=7, height=5, res=300) 
print(ldh2)
dev.off()

save(ldh2, file="lda.plot.RData")

#--- general relationship between LDH and temperature (all fish) ---#
ldh.data2 <- ldh.data %>% 
  mutate(temperature = as.numeric(temperature))
ldh.cmb <- glmmTMB(LDH_ACTIVITY ~ 1 + temperature + TISSUE_MASS_CENTERED + (1|fish_id), 
                   family=gaussian(), 
                   data = ldh.data2, 
                   REML = TRUE)

ldh.cmb.2 <- glmmTMB(LDH_ACTIVITY ~ 1 + poly(temperature, 2) + TISSUE_MASS_CENTERED + (1|fish_id), 
                   family=gaussian(), 
                   data = ldh.data2, 
                   REML = TRUE) 

ldh.cmb.3 <- glmmTMB(LDH_ACTIVITY ~ 1 + poly(temperature, 3) + TISSUE_MASS_CENTERED + (1|fish_id), 
                   family=gaussian(), 
                   data = ldh.data2, 
                   REML = TRUE) 

AIC(ldh.cmb, ldh.cmb.2, ldh.cmb.3, k=2)

saveRDS(ldh.cmb.3, "ldh_cmb_3.RDS")

summary(ldh.cmb.3)
ldh.cmb.3 %>% check_model()
ldh.cmb.3 %>% emtrends(~ temperature, var = "temperature") %>% summary(infer = TRUE)
ldh.cmb.3 %>% performance::r2_nakagawa()

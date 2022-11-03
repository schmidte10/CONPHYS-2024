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
#--- data preparation/manipulation ---# 
ldh2 <- ldh %>%
  clean_names() %>%
  mutate(muscle_type = str_replace(muscle_type, " ", ".")) %>%
  unite("UNIQUE_SAMPLE_ID", c(muscle_type,fish_id,temperature), sep="_", remove = FALSE) %>% 
  separate(creation_time, into=c('DATE','TIME'), sep = " ", remove = FALSE) %>% 
  arrange(sample_id_1, DATE, TIME) 

ldh3 <- ldh2 %>% 
  mutate(TIME = hms(ldh2$TIME)) %>% 
  mutate(TIME = chron(times=ldh2$TIME)) %>% 
  group_by(UNIQUE_SAMPLE_ID, sample_id_1) %>% 
  mutate(TIME_DIFF = TIME - first(TIME)) %>% 
  ungroup() %>% 
  mutate(TIME_DIFF_SECS = period_to_seconds(hms(TIME_DIFF))) %>% 
  mutate(MINUTES = TIME_DIFF_SECS/60) %>% 
  mutate(MINUTES = round(MINUTES, digits = 2)) %>% 
  rename(CUVETTE = sample_id_1) %>% 
  mutate(REGION = substr(fish_id, 1, 1 ), 
         POPULATION = substr(fish_id, 2, 4), 
         SAMPLE_NO = substr(fish_id, 5, 7)) %>% 
  mutate(REGION = case_when( REGION =="L"~ "Leading", 
                             REGION == "C" ~ "Core", 
                             TRUE ~ "na"))

  

#--- making plot ---# 
sample_names_list <- unique(ldh3$UNIQUE_SAMPLE_ID)
Plots <- list()
for (i in sample_names_list) {
  Plots[[i]] <- ldh3 %>% 
    filter(UNIQUE_SAMPLE_ID == i) %>%
    ggplot(aes(MINUTES, result)) + 
    geom_point() +
    facet_wrap(~CUVETTE) + 
    geom_smooth(method = "lm") + 
    theme_bw() + 
    ggtitle(i) + 
    stat_regline_equation(label.y = 0.7) + 
    stat_cor(label.y = 0.6)
  Plots[[i]] = Plots[[i]]
}


#--- making table - quick processing for app---# 
tables <- list()

for (i in sample_names_list) {
  
  LDH_activity <- ldh3 %>% 
    filter(UNIQUE_SAMPLE_ID == i) %>%
    group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
    do({
      mod = lm(result ~ MINUTES, data = .)
      data.frame(Intercept = coef(mod)[1],
                 Slope = coef(mod)[2])
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
    filter(UNIQUE_SAMPLE_ID == i) %>%
    group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
    do({
      mod = lm(result ~ MINUTES, data = .)
      data.frame(Intercept = coef(mod)[1],
                 Slope = coef(mod)[2])
    }) %>%
    ungroup() %>%
    filter(CUVETTE == ("5")) %>% 
    rename(Background = Slope)
  
  final_table <- LDH_activity %>% 
    full_join(distinct(LDH_activity_means[,c(1,5)]), by = "UNIQUE_SAMPLE_ID") %>% 
    full_join(LDH_background[,c(1,4)], by = "UNIQUE_SAMPLE_ID")
  final_table$Mean[duplicated(final_table$Mean)] <- ""
  final_table$Background[duplicated(final_table$Background)] <- ""
  tables[[i]] <- final_table 
  
  
}


tables[["white.muscle_LCKM163_40"]]



#--- Template - use for data analysis---#
LDH_activity <- ldh3 %>% 
  group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
  do({
    mod = lm(result ~ MINUTES, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
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
  rename(Background = Slope)

final_table <- LDH_activity %>% 
  full_join(distinct(LDH_activity_means[,c(1,5)]), by = "UNIQUE_SAMPLE_ID") %>% 
  full_join(LDH_background[,c(1,4)], by = "UNIQUE_SAMPLE_ID") 
final_table$Mean[duplicated(final_table$Mean)] <- ""
final_table$Background[duplicated(final_table$Background)] <- ""

ldh.data <- final_table %>% 
  select(c(UNIQUE_SAMPLE_ID, Mean, Background)) %>% 
  mutate(Mean = as.numeric(Mean), 
         Background = as.numeric(Background),
         LDH_ACTIVITY = Mean - Background) %>% 
  drop_na() %>% 
  inner_join(select(ldh3, c(UNIQUE_SAMPLE_ID, REGION, POPULATION, temperature, fish_id)), by ="UNIQUE_SAMPLE_ID") %>% 
  distinct(UNIQUE_SAMPLE_ID, REGION, POPULATION, .keep_all = TRUE) %>% 
  mutate(LDH_ACTIVITY_POSITIVE = LDH_ACTIVITY * -1, 
         temperature = factor(temperature))

#--- begin data analysis ---# 
hist(ldh.data$LDH_ACTIVITY); hist(scale(ldh.data$LDH_ACTIVITY, center=TRUE, scale=FALSE)*-1) 
# data is not normal - perhaps more samples will end up helping 

ldh.data %>% 
  group_by(fish_id) %>% 
  summarise(sample_size = n(), 
            Min. = min(LDH_ACTIVITY), 
            Max. = max(LDH_ACTIVITY), 
            Mean = mean(LDH_ACTIVITY)) 

#--- models ---# 

ldh.model.1 <- glmmTMB(LDH_ACTIVITY ~ 1 + REGION*temperature + (1|fish_id), 
                       family=gaussian(), 
                       data = ldh.data, 
                       REML = TRUE) 


ldh.model.2 <- glmmTMB(LDH_ACTIVITY_POSITIVE ~ 1 + REGION*temperature + (1|fish_id), 
                       family=Gamma('inverse'), 
                       data = ldh.data, 
                       REML = TRUE)  

ldh.model.3 <- glmmTMB(LDH_ACTIVITY_POSITIVE ~ 1 + REGION*temperature + (1|fish_id) + (1|POPULATION), 
                       family=Gamma(), 
                       data = ldh.data, 
                       REML = TRUE) 

ldh.model.1 %>% check_model()
ldh.model.1 %>% simulateResiduals(plot = TRUE, integerResponse = TRUE) 
ldh.model.1 %>% testResiduals()

#--- partial plots ---# 
ldh.model.1 %>% ggemmeans(~temperature*REGION) %>% plot()
ldh.model.1 %>% summary()
ldh.model.1 %>% confint()
ldh.model.1 %>% performance::r2_nakagawa()

#--- results ---# 
ldh.model.1 %>% emmeans(~ temperature*REGION, type = "response") %>% pairs(by = "temperature") %>% summary(infer = TRUE) 

#--- plot ---# 
newdata <- ldh.model.1 %>% ggemmeans(~temperature|REGION) %>% 
  as.data.frame() %>% 
  rename(TEMPERATURE = x)
  
g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color = group)) + 
  geom_point() + 
  theme_classic(); g1

obs <- ldh.model.1 %>% 
  mutate(Pred = predict(ldh.model.1, re.form=NA), 
         Resid = residuals(ldh.model.1, type = 'response'), 
         Fit = Pred - Resid)

g2 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color = group)) + 
  geom_pointrange(aes(ymin=conf.low, 
                      ymax=conf.high), 
                  shape=19, 
                  size=1, 
                  position = position_dodge(0.2)) + 
  #scale_y_continuous(limits = c(0,0.9), breaks = seq(0, 0.9, by =0.15)) + 
  theme_classic() + ylab("PHA Swelling response (mm)") + 
  scale_color_manual(values=c("#DA3A36","#0D47A1"), name = "Regions"); g2

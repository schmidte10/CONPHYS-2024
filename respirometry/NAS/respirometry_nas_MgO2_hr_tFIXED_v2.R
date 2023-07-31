#--- load libraries ---#  
library(tidyverse) 
library(tidyr)
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
         RESTING_END_TIME = hms(RESTING_ENDTIME),
         MAX_DATE = factor(MAX_DATE), 
         MAX_CHAMBER = factor(MAX_CHAMBER), 
         MAX_SYSTEM = factor(MAX_SYSTEM), 
         MAX_SUMP = factor(MAX_SUMP), 
         MAX_AM_PM = factor(MAX_AM_PM), 
         MAX_START_TIME = hms(MAX_START_TIME), 
         Swim.performance = factor(Swim.performance), 
         MgO2.hr_Net = as.numeric(MgO2.hr_Net), 
         RESTING_RUNTIME_SECONDS = as.numeric(hms(RESTING_RUNTIME))) %>% 
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
      EXP_FISH_ID !="CSUD026_30" & # max. value low 
      EXP_FISH_ID !="CSUD074_28.5" & # fas value low 
      EXP_FISH_ID !="CSUD079_30" &
      EXP_FISH_ID !="CVLA052_27" & #nas value low 
      EXP_FISH_ID !="CVLA054_28.5" & # low max value? 
      EXP_FISH_ID !="CVLA104_27" &
      EXP_FISH_ID !="LCHA113_27" & # poor data quality 
      EXP_FISH_ID !="LCHA113_30" & # poor swim 
      EXP_FISH_ID !="LCHA127_27" & # deceased during experiment
      EXP_FISH_ID !="LCHA114_28.5"  # poor swim 
  )   


#--- exploratory data analysis ---# 

resp4 %>% ggplot(aes(x=TEMPERATURE, y=MgO2.hr_Net, fill = REGION)) + geom_boxplot() 
resp4 %>% ggplot(aes(x=MgO2.hr_Net, fill=TEMPERATURE)) + geom_density(alpha = 0.5) + 
  facet_wrap(~TEMPERATURE)


#--- make sqrt transformed column ---# 
resp4 %>% 
  group_by(REGION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(MgO2.hr_Net), 
                   Max. = max(MgO2.hr_Net), 
                   Mean = mean(MgO2.hr_Net))

pop.sample.size <- resp4 %>% 
  group_by(POPULATION, TEMPERATURE)  %>%    
  dplyr::summarise(sample_size = n(), 
                   Min. = min(MgO2.hr_Net), 
                   Max. = max(MgO2.hr_Net), 
                   Mean = mean(MgO2.hr_Net)); pop.sample.size

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
nas.3 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE +  MAX_SUMP + MAX_CHAMBER + 
                   MAX_AM_PM, 
                 family=gaussian(),
                 data = resp4,
                 REML = FALSE) 

AIC(nas.1, nas.2, nas.3, k=2)

#--- followed by inclusion of random variables
nas.1a <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE) 

nas.1b <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|POPULATION/FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

nas.1c <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|FISH_ID) + (REGION|POPULATION), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

AIC(nas.1, nas.1a, nas.1b, nas.1c, k=2)

#--- Final model ---# 
nas.1c <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED + (1|POPULATIONFISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

#--- saving model ---#
saveRDS(nas.1b, file = "nas_1b.RDS") 

#--- load model ---# 
 nas.1b <- readRDS("nas_1b.RDS") 

#--- investigate model ---#
#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")
check_model(nas.1b)
pha.resid <-  nas.1b %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
nas.1b %>% DHARMa::testResiduals() 

#sim <- simulateResiduals(nas.1a)
#which(residuals(sim) == 1 | residuals(sim) == 0)

nas.1b %>% plot_model(type='eff',  terms=c('TEMPERATURE','REGION'), show.data=TRUE)
nas.1b %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
nas.1b %>% plot_model(type='est')

nas.1b %>% summary()
nas.1b %>% confint()
nas.1b  %>% performance::r2_nakagawa()

nas.1b %>% emmeans(~ TEMPERATURE*REGION, type = "response")  %>% summary(infer=TRUE)
nas.1b %>% emmeans(~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(infer=TRUE)
nas.1b %>% emmeans(~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "REGION") %>% summary(infer=TRUE)
#--- plot ---#
newdata <- nas.1b %>% ggemmeans(~TEMPERATURE|REGION) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

g1 <- ggplot(newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); g1

predict(nas.1b, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(nas.1b, type='response')
#data points based on month/situation/random effects - to get the data points
obs <-  resp4 %>% 
  mutate(Pred=predict(nas.1b, re.form=NA),
         Resid = residuals(nas.1b, type='response'),
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

pdf("nas_1b.pdf", width = 7, height = 5)
print(g2)
dev.off()

jpeg("nas_1b.jpeg", units="in", width=7, height=5, res=300) 
print(g2)
dev.off()

#--- figure for 3MT ---# 
g2 <- ggplot(newdata, aes(y=predicted, x=as.numeric(TEMPERATURE), color=group)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se=F, fill=NA, 
              formula=y ~ poly(x, 3, raw=TRUE), linewidth = 2)+
  xlab("TEMPERATURE") +
  scale_y_continuous(limits = c(6,13), breaks = seq(6, 13, by = 2)) + 
  scale_x_discrete(limits=c("27","28.5","30","31.5"))+
  #scale_x_continuous(limits = c(26,32), breaks = seq(27,31.5, by =1.5))+
  #scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+
  theme_classic() + 
  theme(axis.title.x=element_text(colour="white"), 
        axis.title.y=element_text(colour="white"), 
        axis.text = element_text(color="white"), 
        axis.line = element_line(color="white"), 
        axis.ticks = element_line(color="white"), 
        panel.background = element_rect(fill = "transparent", 
                                        color = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), 
        panel.grid.major = element_line(colour = "NA")) +
  ylab("NET AEROBIC SCOPE (NAS: MgO2/hr)") +
  scale_color_manual(values=c("#DA3A36", "#00FFFF"), labels = c("Cairns (north)","Mackay (south)"),
                     name = "Regions");g2

ggsave('mt_figure2.png',g2,bg='transparent', 
       width = 6, height = 3.2, units = "in", dpi=300, device = 'png')

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
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/import_files/")
#--- load data ---# 
resp <- read.delim("./SummaryData_2022_resp_updated.txt")
# uni computer
setwd("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/respirometry/NAS")


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
      EXP_FISH_ID !="CSUD008_28.5" & # poor swim
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
save(resp4, file="./resp4.RData")

#--- exploratory data analysis ---# 


mass.distr <- resp4 %>% distinct(FISH_ID, .keep_all = TRUE) %>% 
  ggplot(aes(x=MASS, y=REGION, fill=REGION)) + 
  scale_fill_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Low","High"),
                    name = "Latitude") +
  ylab("")+
  geom_density_ridges(scale = 2, jittered_points=TRUE, position = position_points_jitter(height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) + 
  theme_classic() + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) 

pdf(file = "C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/supplemental_figures/Supplemental_figure3.pdf", 
    width = 10, height=6)
mass.distr
dev.off()



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

#---linear ---# 
nas.1 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * TEMPERATURE + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4,
                 REML = FALSE) 

#--- second order polynomial ---# 
nas.1.p2 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE,2) + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4,
                 REML = FALSE) 

#--- third order polynomial ---#
nas.1.p3 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE, 3) + MASS_CENTERED, 
                 family=gaussian(),
                 data = resp4,
                 REML = FALSE) 

AICc(nas.1, nas.1.p2, nas.1.p3, k=2)
#--- followed by inclusion of random variables
nas.1a <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE,2) + MASS_CENTERED + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE) 

nas.1b <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE,2) + MASS_CENTERED + (1|POPULATION/FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

nas.1c <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE,2) + MASS_CENTERED + (1|FISH_ID) + (REGION|POPULATION), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

AICc( nas.1a, nas.1b, nas.1c, k=2)

#--- Final model ---# 
nas.1.p2 <- glmmTMB(MgO2.hr_Net ~ 1+ REGION * poly(TEMPERATURE,2) + MASS_CENTERED + (1|FISH_ID), 
                  family=gaussian(),
                  data = resp4,
                  REML = TRUE)

#--- saving model ---#
saveRDS(nas.1.p2, file = "nas_1_p2.RDS") 

#--- load model ---# 
nas.1.p2 <- readRDS("./nas_1_p2.RDS") 

#--- investigate model ---#
#rest.poly3 <- readRDS("glmmTMB_restpoly3.RDS")
check_model(nas.1.p2)
pha.resid <-  nas.1.p2 %>% 
  DHARMa::simulateResiduals(plot = TRUE, integerResponse = TRUE) 
nas.1.p2 %>% DHARMa::testResiduals() 

#sim <- simulateResiduals(nas.1a)
#which(residuals(sim) == 1 | residuals(sim) == 0)

nas.1.p2 %>% plot_model(type='eff',  terms=c('TEMPERATURE','REGION'), show.data=TRUE)
nas.1.p2 %>% ggemmeans(~TEMPERATURE|REGION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
nas.1.p2 %>% plot_model(type='est')

nas.1.p2 %>% summary()
nas.1.p2 %>% Anova()
nas.1.p2 %>% confint()
nas.1.p2  %>% performance::r2_nakagawa()

nas.1.p2 %>% emmeans(~ TEMPERATURE*REGION, type = "response")  %>% summary(infer=TRUE)
nas.1.p2 %>% emmeans(pairwise ~ TEMPERATURE*REGION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
nas.1.p2 %>% emtrends(var = "TEMPERATURE", type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
#--- plot ---#
nas.emm <- emmeans(nas.1.p2, ~ TEMPERATURE*REGION, 
                  at = list(TEMPERATURE = seq(from=27, to = 31.5, by=.1)))
nas.emm.df=as.data.frame(nas.emm)

nas.newdata <- nas.1.p2 %>% ggemmeans(terms = c("TEMPERATURE[all]","REGION")) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

nas.g1 <- ggplot(nas.newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high), shape=21,
                  position=position_dodge(0.2))+
  theme_classic(); nas.g1

predict(nas.1b, re.form=NA) 
#data points based on month and situation - to get the group means
residuals(nas.1b, type='response')
#data points based on month/situation/random effects - to get the data points

nas.obs <-  resp4 %>% 
  mutate(Pred=predict(nas.1.p2, re.form=NA),
         Resid = residuals(nas.1.p2, type='response'),
         Fit = Pred + Resid)


nas.1.p2.g2 <- ggplot(nas.newdata, aes(y=predicted, x=TEMPERATURE, color=group)) + 
  geom_jitter(data=nas.obs, aes(y=Fit, color=REGION), width=0.05, alpha = 0.3) +
  stat_smooth(method = "lm", 
              formula =y ~ poly(x, 2, raw=TRUE)) + 
  geom_ribbon(aes(x=TEMPERATURE, ymin= conf.low, ymax= conf.high, fill = group), 
              alpha = 0.2, color=NA)+
  scale_y_continuous(limits = c(4,20), breaks = seq(4, 20, by = 2)) + 
  scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+
  theme_classic() + ylab("ABSOLUTE AEROBIC SCOPE (AAS: MgO2/hr)") +
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Cairns (north)","Mackay (south)"),
                     name = "Regions") +
  scale_fill_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Cairns (north)","Mackay (south)"),
                     name = "Regions") + 
  theme(legend.position = "none") + 
  annotate("text", x=31, y= 19, label="P =0.0039", fontface="italic", size=4); nas.1.p2.g2

nas.1.p2.g2 <- ggplot(nas.emm.df, aes(y=emmean, x=TEMPERATURE, color=REGION)) + 
  geom_jitter(data=nas.obs, aes(y=Fit, color=REGION), width=0.05, alpha = 0.3) +
  stat_smooth(method = "lm", 
              formula =y ~ poly(x, 2, raw=TRUE)) + 
  geom_ribbon(aes(x=TEMPERATURE, ymin= lower.CL, ymax= upper.CL, fill = REGION), 
              alpha = 0.2, color=NA)+
  scale_y_continuous(limits = c(4,20), breaks = seq(4, 20, by = 2)) + 
  scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+
  theme_classic() + ylab("ABSOLUTE AEROBIC SCOPE (AAS: MgO2/hr)") +
  scale_color_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Cairns (north)","Mackay (south)"),
                     name = "Regions") +
  scale_fill_manual(values=c("#DA3A36", "#0D47A1"), labels = c("Cairns (north)","Mackay (south)"),
                    name = "Regions") + 
  theme(legend.position = "none") + 
  annotate("text", x=31, y= 19, label="P =0.0039", fontface="italic", size=4); nas.1.p2.g2

#geom_signif(
#y_position = c(4.11+1.5, 5.18+1.5,5.15+1.5,4.66+1.5), xmin = c(0.8, 1.8,2.8,3.8), xmax = c(1.2,2.2,3.2,4.2),
# annotation = c("ns", "ns", "**\np =0.046", "ns"), tip_length = 0.025, color = "black"); g2

pdf("nas_1_p2.pdf", width = 7, height = 5)
print(nas.1.p2.g2)
dev.off()

jpeg("nas_1_p2.jpeg", units="in", width=7, height=5, res=300) 
print(nas.1.p2.g2)
dev.off()


#--- figure 2 ---# 
 resp.plot <- ggarrange(rmr.g2, mmr.g2, nas.1.p2.g2, labels =c("A", "B", "C")); resp.plot 
 resp.plot2 <- plot_grid(rmr.g2, mmr.g2, nas.1.p2.g2, nrow=3, align = "h", axis="bt", rel_widths = c(1,1,1), labels=c("A","B","C")); resp.plot2

 ggsave("C:/Users/jc527762/OneDrive - James Cook University/PhD dissertation/Data/Chapter1_LocalAdaptation/figures/figure2.pdf", width = 12, height = 29, units = 'cm', dpi = 360)
 resp.plot2
 dev.off()















#---mackay, core, cahuvel ---# 
resp7 <- resp4 |> 
  mutate(LOCATION = case_when(POPULATION == "Sudbury Reef" ~ "Core", 
                              POPULATION == "Vlassof Cay" ~ "Core",
                              POPULATION == "Tongue Reef" ~ "Core",
                              POPULATION == "Cockermouth Island" ~ "Mackay (inshore)",
                              POPULATION == "Keswick Island" ~ "Mackay (inshore)",
                              POPULATION == "Chauvel Reef" ~ "Chauvel")) 

nas.pop2 <- glmmTMB(MgO2.hr_Net ~ 1+ poly(TEMPERATURE, 2)*LOCATION + MASS_CENTERED + (1|FISH_ID), 
                    family=gaussian(),
                    data = resp7,
                    REML = FALSE)

nas.pop2 %>% summary()
nas.pop2 %>% Anova()
nas.pop2 %>% confint()
nas.pop2  %>% performance::r2_nakagawa()

nas.pop2 %>% emmeans(~ TEMPERATURE*LOCATION, type = "response")  %>% summary(infer=TRUE)
nas.pop2 %>% emmeans(pairwise ~ TEMPERATURE*LOCATION, type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)
nas.pop2 %>% emtrends(var = "TEMPERATURE", type = "response") %>% pairs(by = "TEMPERATURE") %>% summary(by = NULL, adjust = "tukey", infer=TRUE)


nas.pop2 %>% plot_model(type='eff',  terms=c('TEMPERATURE','LOCATION'), show.data=TRUE)
nas.pop2 %>% ggemmeans(~TEMPERATURE|LOCATION) %>% plot(add.data=TRUE, jitter=c(0.05,0))
nas.pop2 %>% plot_model(type='est')

nas.newdata <- nas.pop2 %>% ggemmeans(~TEMPERATURE|LOCATION) %>%
  as.data.frame %>% 
  dplyr::rename(TEMPERATURE = x)

nas.g1 <- ggplot(nas.newdata, aes(y=predicted, x=TEMPERATURE, color=group)) +
  geom_point()+
  theme_classic(); nas.g1

nas.1.p2.g2_location <- ggplot(nas.newdata, aes(y=predicted, x=TEMPERATURE, color=group)) + 
  stat_smooth(method = "lm", SE=FALSE,
              formula =y ~ poly(x, 2, raw=TRUE)) + 
  #geom_ribbon(aes(x=TEMPERATURE, ymin= conf.low, ymax= conf.high, fill = group), 
              #alpha = 0.4, color=NA)+
  #scale_y_continuous(limits = c(6,13), breaks = seq(6, 13, by = 2)) + 
  #scale_x_continuous(limits = c(26.9, 31.6), breaks = seq(27, 31.5, by = 1.5))+
  theme_classic() + ylab("ABSOLUTE AEROBIC SCOPE (AAS: MgO2/hr)") +
  scale_color_manual(values=c("#DA3A36",  "orange","#0D47A1"), labels = c("Cairns (north)","Mackay (Inshore)", "Mackay (Outshore)"),
                     name = "Regions") +
  scale_fill_manual(values=c("#DA3A36",  "orange", "#0D47A1"))+
  theme(legend.position = "none"); nas.1.p2.g2_location
